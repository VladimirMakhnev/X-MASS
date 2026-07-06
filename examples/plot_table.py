"""Inspect and plot an X-MASS / ABSCO HDF5 cross-section table.

Usage:
    python plot_table.py TABLE.hdf5 [--ip 0] [--it 0] [--ivmr 0] [-o out.png]

Prints the table structure (datasets, shapes, attributes) and saves a
log-scale plot of the selected cross-section slice.
"""
import argparse
import sys

import numpy as np
import h5py
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt


def find_absorption_dataset(f):
    for name in f.keys():
        if name.startswith('Gas_') and name.endswith('_Absorption'):
            return name
    sys.exit('No Gas_XX_Absorption dataset found in the file.')


def find_broadener_dataset(f):
    for name in f.keys():
        # ABSCO-conformant name first, legacy X-MASS name second
        if name.startswith('Broadener_') and (name.endswith('_VMR') or name.endswith('_VMS')):
            return name
    return None


def main():
    parser = argparse.ArgumentParser(description=__doc__.splitlines()[0])
    parser.add_argument('table', help='path to the HDF5 table')
    parser.add_argument('--ip', type=int, default=0, help='pressure index (default 0)')
    parser.add_argument('--it', type=int, default=0, help='temperature index (default 0)')
    parser.add_argument('--ivmr', type=int, default=0, help='broadener VMR index (default 0)')
    parser.add_argument('-o', '--output', default=None,
                        help='output image name (default: <table>.png)')
    args = parser.parse_args()

    with h5py.File(args.table, 'r') as f:
        print('Datasets in %s:'%args.table)
        for name in f.keys():
            ds = f[name]
            attrs = {k: ds.attrs[k] for k in ds.attrs.keys()}
            print('  %-28s %-18s %s'%(name, ds.shape, attrs if attrs else ''))

        ds_name = find_absorption_dataset(f)
        wn = f['Wavenumber'][()]
        pres = f['Pressure'][()]
        temp = f['Temperature'][()]
        broad_name = find_broadener_dataset(f)
        vmr = f[broad_name][()] if broad_name else [0.0]
        xsec = f[ds_name][args.ip, args.it, args.ivmr, :]

    if np.isnan(xsec).all():
        sys.exit('Selected slice is all NaN (failed grid point).')

    label = 'p=%.3g atm, T=%.5g K, VMR=%.2f'%(
        pres[args.ip], np.atleast_2d(temp)[args.ip][args.it], vmr[args.ivmr])
    fig, ax = plt.subplots(figsize=(12, 5))
    ax.plot(wn, xsec, lw=0.6, label=label)
    ax.set_yscale('log')
    ax.set_xlabel(r'Wavenumber (cm$^{-1}$)')
    ax.set_ylabel(r'Cross-section (cm$^2$/molecule)')
    ax.set_title('%s : %s'%(args.table, ds_name))
    ax.legend()

    out = args.output or (args.table.rsplit('.', 1)[0] + '.png')
    fig.savefig(out, dpi=200, bbox_inches='tight')
    print('Plot saved to %s'%out)


if __name__ == '__main__':
    main()
