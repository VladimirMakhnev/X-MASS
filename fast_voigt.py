"""Vectorized Voigt evaluation for the X-MASS fast path (Stage 2).

The line list is split once per worker into:
  - "slow" lines: any beyond-Voigt parameter (speed dependence, HT, or
    profile-specific line mixing) is present in the table for the chosen
    profile. These keep going through HAPI, evaluated on the exact same
    wavenumber grid via a subset table (XMASS_SLOW).
  - "fast" lines: plain Voigt (optionally with first-order Rosenkranz line
    mixing from y_air/y_self), evaluated here in line-blocks with
    scipy.special.wofz.

The classification mirrors HAPI's `ladder` fallback semantics: a parameter
lookup falls through to the Voigt parametrization exactly when the profile-
specific column is absent or masked, so "no such column value" == "HAPI
would have used the Voigt route for this line". The temperature/pressure
scaling formulas below replicate hapi's PRESSURE_INDUCED_ENVDEP presets
(PowerLaw/LinearLaw with the same defaults and the n_self -> n_air
fallback), and the Doppler width uses hapi's own CGS constants.
"""
import numpy as np
from scipy.special import wofz

import hapi as hapi1

T_REF = 296.0
P_REF = 1.0
SQRT_LN2 = np.sqrt(np.log(2.0))
SQRT_PI = np.sqrt(np.pi)

# ~50 bytes per element of the (lines x window) complex block; 2M elements
# keeps a worker near 100 MB of scratch
DEFAULT_BLOCK_ELEMS = 2_000_000

# profile-specific column prefixes whose presence sends a line to the slow
# (HAPI) path; lowercase, matched against the table column names. Prefix
# matching also covers the HT "multitemp" reference-temperature variants.
SLOW_PREFIXES = {
    'SDVoigt': ['gamma_sdv_', 'delta_sdv_', 'deltap_sdv_', 'sd_'],
    'HT':      ['gamma_ht_', 'delta_ht_', 'deltap_ht_', 'nu_ht_',
                'kappa_ht_', 'eta_ht_'],
    'Voigt':   [],
}
# Line-mixing columns do NOT make a line slow for Voigt/SDVoigt: with no
# speed-dependence widths present, HAPI evaluates such lines as Voigt plus
# the first-order Y term (pcqsdhc with Gamma2=Delta2=0), which the fast path
# reproduces analytically -- only the Y source columns differ per profile
# (y_air/y_self for Voigt, Y_SDV_*_296 for SDVoigt).
SLOW_PREFIXES_LM = {
    'SDVoigt': [],
    'HT':      ['y_ht_'],
    'Voigt':   [],
}

# per-profile database columns feeding the analytic Y term
Y_COLUMNS = {
    'Voigt':   ('y_air', 'n_y_air', 'y_self', 'n_y_self'),
    'SDVoigt': ('Y_SDV_air_296', 'n_Y_SDV_air_296',
                'Y_SDV_self_296', 'n_Y_SDV_self_296'),
    'HT':      ('Y_HT_air_296', 'n_Y_HT_air_296',
                'Y_HT_self_296', 'n_Y_HT_self_296'),
}


def _missing_mask(col):
    """True where the column value is missing (masked or NaN)."""
    if isinstance(col, np.ma.MaskedArray):
        filled = np.ma.filled(col.astype(float), np.nan)
        return np.ma.getmaskarray(col) | ~np.isfinite(filled)
    arr = np.asarray(col, dtype=float)
    return ~np.isfinite(arr)


def _filled_column(data, name, n, fill=0.0):
    """Numeric column with missing entries replaced by `fill`;
    all-`fill` when the column is not in the table at all."""
    if name not in data:
        return np.full(n, fill, dtype=float), np.ones(n, dtype=bool)
    col = data[name]
    miss = _missing_mask(col)
    if isinstance(col, np.ma.MaskedArray):
        arr = np.ma.filled(col.astype(float), fill)
    else:
        arr = np.asarray(col, dtype=float).copy()
    arr[miss] = fill
    return arr, miss


def classify_lines(data, profile_name, flag_lm):
    """Boolean mask, True = slow (HAPI) path, for the chosen profile."""
    n = len(data['nu'])
    slow = np.zeros(n, dtype=bool)
    prefixes = list(SLOW_PREFIXES.get(profile_name, []))
    if flag_lm:
        prefixes += SLOW_PREFIXES_LM.get(profile_name, [])
    if not prefixes:
        return slow
    for key in list(data.keys()):
        kl = str(key).lower()
        if any(kl.startswith(p) for p in prefixes):
            slow |= ~_missing_mask(data[key])
    return slow


def make_slow_table(src, dst, mask_slow):
    """Register an in-memory HAPI table holding only the slow lines.
    Masked columns keep their masks so HAPI's missing-value handling
    (np.ma.core.MaskedConstant checks) behaves identically."""
    tbl = hapi1.LOCAL_TABLE_CACHE[src]
    header = dict(tbl['header'])
    header['table_name'] = dst
    header['number_of_rows'] = int(np.count_nonzero(mask_slow))
    data_dst = hapi1.CaselessDict()
    for key in list(tbl['data'].keys()):
        col = tbl['data'][key]
        if isinstance(col, np.ma.MaskedArray):
            data_dst[key] = col[mask_slow]
        else:
            data_dst[key] = np.asarray(col)[mask_slow]
    hapi1.LOCAL_TABLE_CACHE[dst] = {'header': header, 'data': data_dst}


def build_fast_context(data, mask_fast, flag_lm, profile_name='Voigt'):
    """Contiguous per-line arrays for the fast lines."""
    n_all = len(data['nu'])
    idx = np.flatnonzero(np.asarray(mask_fast))

    def take(name, fill=0.0):
        arr, _ = _filled_column(data, name, n_all, fill)
        return np.ascontiguousarray(arr[idx])

    nu = take('nu')
    ctx = {
        'nu': nu,
        'sw': take('sw'),
        'elower': take('elower'),
        'gamma_air': take('gamma_air'),
        'gamma_self': take('gamma_self'),
        'delta_air': take('delta_air'),
        'deltap_air': take('deltap_air'),
        'delta_self': take('delta_self'),
        'deltap_self': take('deltap_self'),
    }

    # temperature exponent: n_self where present, n_air otherwise
    # (hapi's 'Lorentz 1'/'Lorentz 2' lookup cases)
    n_air, _ = _filled_column(data, 'n_air', n_all, 0.0)
    ctx['n_air'] = np.ascontiguousarray(n_air[idx])
    n_self, miss_self = _filled_column(data, 'n_self', n_all, 0.0)
    n_self_eff = np.where(miss_self, n_air, n_self)
    ctx['n_self_eff'] = np.ascontiguousarray(n_self_eff[idx])

    if flag_lm:
        src_names = Y_COLUMNS[profile_name]
        for name, src_name in zip(('y_air', 'n_y_air', 'y_self', 'n_y_self'),
                                  src_names):
            ctx[name] = take(src_name)   # data dict lookups are caseless
    ctx['flag_lm_columns'] = bool(flag_lm)

    # Doppler HWHM = doppler_coef * sqrt(T); replicates hapi's
    # calculate_parameter_GammaD (CGS constants, mass in grams)
    mols = np.asarray(data['molec_id'], dtype=int)[idx]
    isos = np.asarray(data['local_iso_id'], dtype=int)[idx]
    cMassMol = 1.66053873e-27
    coef = np.empty(len(idx), dtype=float)
    iso_groups = []
    for (m, i) in sorted(set(zip(mols.tolist(), isos.tolist()))):
        sel = (mols == m) & (isos == i)
        mass_g = hapi1.molecularMass(m, i) * cMassMol * 1000.
        coef[sel] = np.sqrt(2*hapi1.cBolts*np.log(2.)/mass_g/hapi1.cc**2)
        iso_groups.append((m, i, sel))
    ctx['doppler_coef'] = coef * nu
    ctx['iso_groups'] = iso_groups
    return ctx


def xsec_fast_voigt(ctx, wngrid, T, p, diluent, flag_lm,
                    wing=25.0, block_elems=DEFAULT_BLOCK_ELEMS):
    """Cross-section (cm2/molecule) of the fast lines on `wngrid`."""
    nu = ctx['nu']
    n = len(nu)
    sigma = np.zeros(len(wngrid))
    if n == 0:
        return sigma

    ab = float(diluent.get('air', 0.))
    sf = float(diluent.get('self', 0.))
    tr = T_REF / T

    # line intensities at T (TIPS partition sums per isotopologue)
    S = np.empty(n)
    for (m, i, sel) in ctx['iso_groups']:
        SigmaT = hapi1.partitionSum(m, i, T)
        SigmaTref = hapi1.partitionSum(m, i, T_REF)
        S[sel] = hapi1.EnvironmentDependency_Intensity(
            ctx['sw'][sel], T, T_REF, SigmaT, SigmaTref,
            ctx['elower'][sel], nu[sel])

    # hapi PowerLaw / LinearLaw environment dependences, summed over diluent
    gamma0 = (p/P_REF)*(ab*ctx['gamma_air']*tr**ctx['n_air']
                        + sf*ctx['gamma_self']*tr**ctx['n_self_eff'])
    delta0 = (p/P_REF)*(ab*(ctx['delta_air'] + ctx['deltap_air']*(T - T_REF))
                        + sf*(ctx['delta_self'] + ctx['deltap_self']*(T - T_REF)))
    if flag_lm and ctx['flag_lm_columns']:
        y = (p/P_REF)*(ab*ctx['y_air']*tr**ctx['n_y_air']
                       + sf*ctx['y_self']*tr**ctx['n_y_self'])
    else:
        y = None

    gammaD = ctx['doppler_coef']*np.sqrt(T)
    cte = SQRT_LN2/gammaD

    # per-line hard windows, matching hapi's bisect(Omegas, nu +/- wing)
    li = np.searchsorted(wngrid, nu - wing, side='right')
    ri = np.searchsorted(wngrid, nu + wing, side='right')

    b0 = 0
    while b0 < n:
        if ri[b0] <= li[b0]:      # window entirely outside the grid
            b0 += 1
            continue
        b1 = b0 + 1
        while b1 < n and (b1 + 1 - b0)*(ri[b1] - li[b0]) <= block_elems:
            b1 += 1
        lo, hi = li[b0], ri[b1-1]
        blk = slice(b0, b1)

        x = wngrid[lo:hi][None, :]
        z = ((x - (nu[blk] + delta0[blk])[:, None])
             + 1j*gamma0[blk][:, None]) * cte[blk][:, None]
        w = wofz(z)
        if y is not None:
            prof = (w.real + y[blk][:, None]*w.imag)
        else:
            prof = w.real
        prof *= (cte[blk]/SQRT_PI)[:, None]

        cols = np.arange(lo, hi)
        inwin = (cols[None, :] >= li[blk][:, None]) & (cols[None, :] < ri[blk][:, None])
        sigma[lo:hi] += S[blk] @ np.where(inwin, prof, 0.0)
        b0 = b1
    return sigma
