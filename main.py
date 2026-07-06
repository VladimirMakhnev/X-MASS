import numpy as np
import time
import sys
import os
import argparse
import logging
import cProfile
import pstats
import initial
import hdf5_io
import core_calcs

try:
    import winsound
    def _beep(freq, dur):
        winsound.Beep(freq, dur)
except ImportError:          # winsound exists on Windows only
    def _beep(freq, dur):
        pass

LOG = logging.getLogger('xmass')

XMASSSEC_VERSION = '0.9.1'; __version__ = XMASSSEC_VERSION
XMASSSEC_HISTORY = [
'INITIATION OF INPUT FILE WITH PARAMETERS 31.01.23 (ver. 0.1)',
'CREATION OF HDF5 FILE + SOME EXCEPTIONS HANDLING (ver. 0.2)',
'CLOSING HDF5 FILE AND SATURATION OF ATTRIBUTES (ver. 0.2.1)',
'SATURATION OF ATTRIBUTE (ROOT) (ver. 0.2.2)',
'INPUTS FOR P,T AND SATURATION OF ATTRUBUTES (ver. 0.2.3)',
'INPUTS FOR VMS, WN AND SATURATION OF ATTRUBUTES + OUTPUT LOG FILE (ver. 0.2.4)',
'CALCULATING X-SEC (ver. 0.3)',
'CALCULATING X-SEC, DEBUG INFO (ver.0.3.1)',
'FIX: MISSING WN RANGE IN HDF5 FILE (ver. 0.3.2)',
'FIX: NEW EXCEPTION IN X-SEC RELATED TO NAN VALUES OF P,T,VMS (ver. 0.3.3)',
'FIRST ITERATION OF PARALLEL ATTEMPT (ver. 0.4)',
'FIX: FLOAT64 FOR DTYPE (ver.0.4.1)',
'A NEW STRUCTURE FOR HYDRA USES (ver. 0.5)',
'MULTIPLE PARALLEL OPTIONS (ver. 0.6)',
'COMMENTS AND STYLE (ver. 0.7)',
'NEW INPUT FILE, PROFILE PICK (ver. 0.8)',
'TIME PROFILER (ver. 0.8.1)',
'LINE MIXING (ROSENKRANZ), MT_CKD PEDESTAL REMOVAL (ver. 0.8.2)',
'PER-WORKER TABLE CACHE, DIRECT HDF5 WRITES, NAN FAILURE ROWS (ver. 0.9)',
'CLI, LOGGING, VALIDATE-ONLY MODE (ver. 0.9.1)'
]


def parse_cli():
    parser = argparse.ArgumentParser(
        prog='x-mass',
        description='X-MASS: batch calculation of ABSCO-format absorption '
                    'cross-section tables from HITRAN via HAPI.')
    parser.add_argument('config', nargs='?', default='filenames.inp',
                        help='master input file listing the other input files '
                             '(default: filenames.inp)')
    parser.add_argument('--method', choices=['PLAIN', 'PC', 'MULTITHREADING'],
                        help='override the calculation method from the config')
    parser.add_argument('--cores', type=int, metavar='N',
                        help='override Number_cores from params file')
    parser.add_argument('--keep-dat', action='store_true',
                        help='also write the legacy per-point .dat text files')
    parser.add_argument('--validate-only', action='store_true',
                        help='read and check all inputs, print a run summary '
                             'and exit without calculating')
    parser.add_argument('--quiet', action='store_true',
                        help='no console output (everything goes to output.log only)')
    parser.add_argument('--version', action='version',
                        version='X-MASS %s'%XMASSSEC_VERSION)
    return parser.parse_args()


def setup_logging(quiet):
    handlers = [logging.FileHandler('output.log', mode='w')]
    if (not quiet):
        handlers.append(logging.StreamHandler(sys.stdout))
    logging.basicConfig(level=logging.INFO, format='%(message)s', handlers=handlers)


def print_banner():
    LOG.info('**************************************************************************')
    LOG.info('XMASSSECTION: program to calculate and store cross-sections in HDF5 files.')
    LOG.info('**************************************************************************')
    LOG.info('X-MASS-SEC version: %s'%(XMASSSEC_VERSION))
    LOG.info('')
    LOG.info('           MIT license: Copyright 2024 HITRAN team, see more at http://hitran.org. ')
    LOG.info('')
    LOG.info('           If you use X-MASS in your research or software development,')
    LOG.info('           please cite it using the following reference:')
    LOG.info('           V.Yu. Makhnev, I.E. Gordon, L.S. Rothman, R.J. Hargreaves')
    LOG.info('           DOI: ')
    LOG.info('')


def run_summary(args, METHOD, ParametersCalculation, Np, Ntt, Nvms, Nwn, HDF5FileName, VMSs):
    wn_begin = float(ParametersCalculation[3][1])
    wn_end = float(ParametersCalculation[4][1])
    n_points = Np*Ntt*Nvms
    est_bytes = n_points*Nwn*8
    LOG.info('*** RUN SUMMARY *******************************')
    LOG.info('Config file:        %s'%args.config)
    LOG.info('Molecule id:        %s (%s)'%(ParametersCalculation[10][1], ParametersCalculation[8][1]))
    LOG.info('Spectral range:     %.2f - %.2f cm-1, %d points (step %.4g cm-1)'%(
             wn_begin, wn_end, Nwn, (wn_end-wn_begin)/(Nwn-1)))
    LOG.info('Grid:               %d pressures x %d temperatures x %d VMR values = %d points'%(
             Np, Ntt, Nvms, n_points))
    LOG.info('Profile:            %s'%ParametersCalculation[18][1])
    LOG.info('Line mixing:        %s'%('ON' if initial.readSwitchByName(ParametersCalculation, 'Line_mixing') else 'OFF'))
    LOG.info('Pedestal removal:   %s'%('ON' if initial.readSwitchByName(ParametersCalculation, 'Remove_pedestal') else 'OFF'))
    LOG.info('Method:             %s (cores: %s of %d available)'%(
             METHOD, ParametersCalculation[17][1], os.cpu_count()))
    LOG.info('Output:             %s (uncompressed data size: %.2f GB)'%(
             HDF5FileName, est_bytes/1024.**3))
    if (np.min(VMSs) < 0.0) or (np.max(VMSs) > 1.0):
        LOG.warning('WARNING: VMR values outside [0, 1]: %s'%np.array2string(np.asarray(VMSs)))
    if (METHOD == 'MULTITHREADING') and (int(ParametersCalculation[17][1]) > os.cpu_count()):
        LOG.warning('WARNING: Number_cores exceeds available CPUs')
    LOG.info('***********************************************')


if __name__ == "__main__":

    args = parse_cli()
    setup_logging(args.quiet)
    print_banner()

    LOG.info("Timer started")
    t_begin = time.time()
    profiler = cProfile.Profile()
    profiler.enable()

    LOG.info("***********************************************")
    LOG.info("*** PARAMETERS HANDLING ***********************")
    LOG.info("***********************************************")

    argvlen = ((open(args.config, mode='r')).read()).split()
    LOG.info('%s'%argvlen)
    (PARAM_FILENAME := argvlen[1]) if (len(argvlen)>1) else (PARAM_FILENAME := 'params.inp')
    (PRES_FILENAME  := argvlen[2]) if (len(argvlen)>2) else (PRES_FILENAME := 'pres_pRT.inp')
    (TEMP_FILENAME  := argvlen[3]) if (len(argvlen)>3) else (TEMP_FILENAME := 'temps_pRT.inp')
    (VMS_FILENAME   := argvlen[4]) if (len(argvlen)>4) else (VMS_FILENAME := 'vms.inp')
    (WN_FILENAME    := argvlen[5]) if (len(argvlen)>5) else (WN_FILENAME := 'wn.inp')
    (HDF5FileName   := argvlen[6]) if (len(argvlen)>6) else (HDF5FileName := '01.H2O.SDV.HITRAN2020.25wing.hdf5')
    (METHOD         := argvlen[7]) if (len(argvlen)>7) else (METHOD := 'MULTITHREADING')

    # command-line overrides
    if (args.method):
        METHOD = args.method
    ParametersCalculation = initial.openParametersFile(PARAM_FILENAME)
    if (args.cores is not None):
        initial.setParamByName(ParametersCalculation, 'Number_cores', str(args.cores))
    if (args.keep_dat):
        initial.setParamByName(ParametersCalculation, 'Keep_dat', 'ON')
    LOG.info('%s'%ParametersCalculation)

    # opening pressure file
    Pressures, Np = initial.openPressure(PRES_FILENAME)

    # Opening temperature file
    (Temps, Npp, Ntt) = initial.openTemp(TEMP_FILENAME, Np)

    # Opening volume mixing ratio file
    (VMSs, Nvms) = initial.openVMS(VMS_FILENAME)

    # Opening wavenumber range file
    (WNs, Nwn) = initial.openXgenetareWn(WN_FILENAME, ParametersCalculation)

    run_summary(args, METHOD, ParametersCalculation, Npp, Ntt, Nvms, Nwn, HDF5FileName, VMSs)
    if (args.validate_only):
        LOG.info('Validation finished: inputs are consistent. Exiting (--validate-only).')
        sys.exit(0)

    LOG.info("***********************************************")
    LOG.info("*** OPENING HDF5 FILE *************************")
    LOG.info("***********************************************")

    # Initialazing the core HDF5 file
    co_hdf5 = hdf5_io.OpenHDF5(HDF5FileName, ParametersCalculation, Pressures, Temps, VMSs, WNs, Npp, Ntt, Nvms, Nwn)
    dataset_name = 'Gas_%02d_Absorption'%(int(ParametersCalculation[10][1]))

    # Constructing the flat task list: one (ip,it,iv,p,T,vms) per grid point
    tasks = initial.mergeParamsIndexed(Pressures, Temps, VMSs)

    LOG.info("***********************************************")
    LOG.info("*** CALCULATIONS PART *************************")
    LOG.info("***********************************************")

    # Calculation part: cross-sections are written into the HDF5 dataset
    # slice by slice as they are computed
    failures = core_calcs.ParallelPart(tasks, ParametersCalculation, Nwn, co_hdf5, dataset_name, METHOD, quiet=args.quiet)

    if (failures):
        LOG.warning('WARNING: %d of %d grid points FAILED; their rows are left as NaN in %s.'%(len(failures), len(tasks), HDF5FileName))
        for (ip, it, iv, err) in failures:
            LOG.warning('    (ip=%d, it=%d, iv=%d): %s'%(ip, it, iv, err.strip().splitlines()[-1]))

    # Closing the HDF5 file
    hdf5_io.CloseHDF5(co_hdf5)

    profiler.disable()
    # profiling results go to their own files, not into the run log
    profiler.dump_stats("profiling_results.prof")
    with open("profiling_results.txt", 'w') as pf:
        stats = pstats.Stats("profiling_results.prof", stream=pf)
        pf.write("Sorted by cumulative time:\n")
        stats.sort_stats("cumulative").print_stats()
    LOG.info('Profiler output: profiling_results.prof / profiling_results.txt')

    t_end = time.time()
    LOG.info('Time taken: %d seconds'%(t_end-t_begin))
    _beep(261, 400)
    _beep(329, 400)
    _beep(392, 400)
    _beep(523, 700)
    LOG.info('\nDone.')

    # nonzero exit when some grid points failed (their rows are NaN)
    sys.exit(1 if failures else 0)
