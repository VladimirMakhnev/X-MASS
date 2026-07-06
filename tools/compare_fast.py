"""Validate the Stage-2 fast Voigt path against the pure-HAPI reference.

Run from a directory that contains a fetched HITRAN2020 table and a
params.inp (e.g. after any X-MASS run):

    python ../tools/compare_fast.py [--points "p,T,vms;p,T,vms;..."]

For every requested (p, T, VMR) point the cross-section is computed twice
through the actual X-MASS worker code path (core_calcs.CalculateXsec) —
once with the fast path, once without — and compared:

    peak-normalized max error   <= 1e-5
    pointwise |rel| where sigma > 1e-6*peak   <= 1e-4
    integrated-intensity rel difference        <= 1e-5

Tolerance rationale (measured, 2026-07): hapi's working CPF is hum1_wei
(Humlicek region 1 + Weideman), whose relative error reaches ~4e-6 for
|z| ~ 15-25 (near-core of pressure-broadened lines), while scipy's wofz
is ~1e-13-accurate. With per-line parameters verified bitwise between
the two paths, residuals at this level are the REFERENCE's CPF error --
the fast path is the more accurate of the two.
"""
import argparse
import os
import sys
import time

import numpy as np

sys.path.insert(0, os.path.join(os.path.dirname(os.path.abspath(__file__)), '..'))
import initial
import core_calcs

TOL_PEAK = 1e-5
TOL_POINT = 1e-4
TOL_INTEGRAL = 1e-5
FLOOR = 1e-6

DEFAULT_POINTS = "1e-4,200,0.0;1e-4,340,0.5;1.0,200,0.5;1.0,340,0.0;0.1,296,1.0"


def compute(param, Nwn, fast, p, T, v):
    flags = {
        'line_mixing':     initial.readSwitchByName(param, 'Line_mixing'),
        'remove_pedestal': initial.readSwitchByName(param, 'Remove_pedestal'),
        'keep_dat':        False,
        'fast':            fast,
        'quiet_workers':   False,
    }
    core_calcs.WORKER.clear()
    core_calcs._init_worker('HITRAN2020', param, Nwn, flags)
    if core_calcs.WORKER.get('init_error'):
        sys.exit('worker init failed:\n%s' % core_calcs.WORKER['init_error'])
    t0 = time.time()
    ip, it, iv, coef, err = core_calcs.CalculateXsec((0, 0, 0, p, T, v))
    dt = time.time() - t0
    if err is not None:
        sys.exit('point (%g,%g,%g) failed:\n%s' % (p, T, v, err))
    return coef, dt


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument('--points', default=DEFAULT_POINTS,
                    help='semicolon-separated p,T,vms triples (default: "%(default)s")')
    ap.add_argument('--params', default='params.inp')
    args = ap.parse_args()

    param = initial.openParametersFile(args.params)
    Nwn = int(open('wn.inp').readline()) if os.path.exists('wn.inp') \
        else int(initial.readParamByName(param, 'Wavenumber'))
    print('profile=%s  LM=%s  pedestal=%s  Nwn=%d'
          % (param[18][1], initial.readSwitchByName(param, 'Line_mixing'),
             initial.readSwitchByName(param, 'Remove_pedestal'), Nwn))

    all_ok = True
    for triple in args.points.split(';'):
        p, T, v = (float(x) for x in triple.split(','))
        ref, t_ref = compute(param, Nwn, False, p, T, v)
        n_slow = None
        fast, t_fast = compute(param, Nwn, True, p, T, v)
        n_slow = core_calcs.WORKER.get('N_SLOW')
        n_fast = len(core_calcs.WORKER['FAST']['nu']) if core_calcs.WORKER.get('FAST') else 0

        peak = ref.max()
        d = np.abs(fast - ref)
        err_peak = d.max()/peak
        m = ref > FLOOR*peak
        err_point = np.max(d[m]/ref[m]) if m.any() else 0.0
        err_int = abs(fast.sum() - ref.sum())/ref.sum()
        ok = (err_peak <= TOL_PEAK) and (err_point <= TOL_POINT) and (err_int <= TOL_INTEGRAL)
        all_ok &= ok
        print('p=%-8g T=%-6g vmr=%-4g | fast/slow lines %d/%d | '
              'peak-norm %.2e  pointwise %.2e  integral %.2e | ref %.2fs fast %.2fs | %s'
              % (p, T, v, n_fast, n_slow, err_peak, err_point, err_int,
                 t_ref, t_fast, 'OK' if ok else 'FAIL'))

    print('GATE:', 'PASSED' if all_ok else 'FAILED')
    return 0 if all_ok else 1


if __name__ == '__main__':
    sys.exit(main())
