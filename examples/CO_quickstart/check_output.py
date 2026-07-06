"""Smoke test for the CO quickstart example.

Run the example first:

    cd examples/CO_quickstart
    python ../../main.py

then check the result:

    python check_output.py

Structural checks are exact; the physical reference values were computed
with HITRAN2020 (fetched 2026-07) and HAPI 1.3.0.0 and are compared with
a 5% tolerance, so future HITRAN line-list updates should not break them.
"""
import sys
import numpy as np
import h5py

FNAME = 'CO_quickstart.hdf5'

# reference values (HITRAN2020, HAPI 1.3.0.0, SDVoigt + line mixing,
# 25 cm-1 wing with fetch margin)
REF_GLOBAL_MAX = 2.0141e-17     # cm2/molecule, p=0.1 atm, T=250 K, air
REF_GLOBAL_MAX_WN = 2165.6      # cm-1
REF_PROBE_PEAK = 2.3267e-18     # cm2/molecule at p=1 atm, T=300 K, air
REF_PROBE_PEAK_WN = 2169.2      # cm-1
RTOL = 0.05

def close(a, b, rtol=RTOL):
    return abs(a - b) <= rtol*abs(b)

def main():
    with h5py.File(FNAME, 'r') as f:
        for name in ('Gas_05_Absorption', 'Pressure', 'Temperature',
                     'Broadener_00_VMR', 'Wavenumber', 'Gas_Index', 'Broadener_Index'):
            assert name in f, 'missing dataset: %s'%name
        x = f['Gas_05_Absorption'][()]
        wn = f['Wavenumber'][()]
        assert f['Gas_05_Absorption'].attrs['units'] == 'cm2/molecule'

    assert x.shape == (2, 2, 2, 2001), 'unexpected shape: %s'%(x.shape,)
    assert np.isfinite(x).all(), 'NaN rows present: some grid points failed'
    assert (x >= 0).all(), 'negative cross-sections'

    i = np.unravel_index(x.argmax(), x.shape)
    ok_max = close(x.max(), REF_GLOBAL_MAX) and abs(wn[i[3]] - REF_GLOBAL_MAX_WN) < 0.5
    print('global max: %.4e at %.1f cm-1 (reference %.4e at %.1f)  %s'
          % (x.max(), wn[i[3]], REF_GLOBAL_MAX, REF_GLOBAL_MAX_WN,
             'OK' if ok_max else 'MISMATCH'))

    probe = x[1, 1, 1, :]
    k = probe.argmax()
    ok_probe = close(probe[k], REF_PROBE_PEAK) and abs(wn[k] - REF_PROBE_PEAK_WN) < 0.5
    print('probe peak: %.4e at %.1f cm-1 (reference %.4e at %.1f)  %s'
          % (probe[k], wn[k], REF_PROBE_PEAK, REF_PROBE_PEAK_WN,
             'OK' if ok_probe else 'MISMATCH'))

    if not (ok_max and ok_probe):
        print('CHECK FAILED')
        return 1
    print('ALL CHECKS PASSED')
    return 0

if __name__ == '__main__':
    sys.exit(main())
