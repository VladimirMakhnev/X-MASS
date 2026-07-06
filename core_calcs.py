from initial import FLAG_DEBUG_PRINT, readSwitchByName

import numpy as np

#from hapi2 import *
#from hapi import *
import hapi as hapi1
# import hapi2

import sys
import os
import logging
import contextlib
import traceback

LOG = logging.getLogger('xmass')


from multiprocessing import Pool



import asyncio
import nest_asyncio
nest_asyncio.apply()
def background(f):
    def wrapped(*args, **kwargs):
        return asyncio.get_event_loop().run_in_executor(None, f, *args, **kwargs)

    return wrapped

# Absolute line-wing cutoff (cm-1). The value of 25 cm-1 follows the Earth
# remote-sensing convention and matches the "local line shape" definition of
# the MT_CKD continuum model (Mlawer et al., JQSRT 306, 108645, 2023).
WING_WN = 25.

# Extra margin (cm-1) added to the HITRAN fetch range so that lines centered
# just outside [wn_begin, wn_end] still contribute their wings to the edge
# points of the output grid.
FETCH_MARGIN = WING_WN

# per-process worker state, populated once by _init_worker
WORKER = {}

def _init_worker(tab_name, param, Nwn, flags):
    """Pool initializer: make the line list available in this process ONCE.

    On fork platforms (Linux) the table is usually inherited from the parent
    and the disk re-parse is skipped; on spawn platforms (Windows) each worker
    parses the storage files a single time. Never raises: a raising
    initializer would make multiprocessing respawn workers in a loop, so
    errors are stored and reported by every task instead."""
    try:
        if (flags.get('quiet_workers') and not FLAG_DEBUG_PRINT):
            sys.stdout = open(os.devnull, 'w')   # silence per-point HAPI chatter
        if (tab_name not in hapi1.LOCAL_TABLE_CACHE):
            hapi1.storage2cache(tab_name)
        wn_begin = float(param[3][1])
        wn_end = float(param[4][1])
        WORKER.update(tab_name=tab_name, param=param, Nwn=Nwn, flags=flags,
                      wngrid=np.linspace(wn_begin, wn_end, Nwn),
                      init_error=None)
    except Exception:
        WORKER['init_error'] = traceback.format_exc()

# single source of truth for the per-point text-file name (the historical
# pattern is kept verbatim so existing tooling can still glob for it)
def dat_filename(param, pres, Temp, VMS):
    IndexMol = int(param[10][1])
    IndexBroad = int(param[15][1])
    return './datafiles/%06.2fT_Id%02d_%06.4eatm_IdBroad%02d_%06.4fVMS_H2O_SDV_hitran2020.dat'%(Temp,IndexMol,pres,IndexBroad,VMS)

def _consume_results(results_iter, ds_coef, hdf5file, Nwn, n_total):
    """Write each finished (ip,it,iv) slice into the HDF5 dataset as results
    arrive, so parent memory stays at one slice regardless of grid size.
    Failed points are left at the dataset fill value (NaN).
    Returns the list of failures [(ip, it, iv, error_text), ...]."""
    failures = []
    n_done = 0
    for (ip, it, iv, coef, err) in results_iter:
        n_done += 1
        if (err is None) and (coef is not None) and (coef.shape == (Nwn,)):
            ds_coef[ip, it, iv, :] = coef
        else:
            if (err is None):
                err = 'bad result shape: %s'%(str(getattr(coef, 'shape', None)))
            failures.append((ip, it, iv, err))
            LOG.warning('FAILED point (ip=%d, it=%d, iv=%d):\n%s'%(ip, it, iv, err))
        if (n_done % 16 == 0):
            hdf5file.flush()
        LOG.info('Progress: %d/%d points done'%(n_done, n_total))
    hdf5file.flush()
    return failures

# Lorentz pedestal ("plinth") spectrum: for every line, the value of its
# Lorentz profile at +/-wing from the (pressure-shifted) line center, spread
# as a constant over the truncation window. Subtracting it from the truncated
# cross-section yields the local-line contribution required for consistency
# with the MT_CKD continuum without double-counting of absorption.
def PedestalSpectrum(wn_grid, tab_name, Temp, pres, diluent, wing):
    data = hapi1.LOCAL_TABLE_CACHE[tab_name]['data']
    nu = np.array(data['nu'], dtype=float)
    sw = np.array(data['sw'], dtype=float)
    elower = np.array(data['elower'], dtype=float)
    gamma_air = np.array(data['gamma_air'], dtype=float)
    gamma_self = np.array(data['gamma_self'], dtype=float)
    n_air = np.array(data['n_air'], dtype=float)
    delta_air = np.array(data['delta_air'], dtype=float)
    mols = np.array(data['molec_id'], dtype=int)
    isos = np.array(data['local_iso_id'], dtype=int)

    Tref = 296.
    # line intensities at Temp (TIPS partition sums per isotopologue)
    S = np.zeros(len(nu))
    for (tmol, tiso) in set(zip(mols.tolist(), isos.tolist())):
        mask = (mols == tmol) & (isos == tiso)
        SigmaT = hapi1.partitionSum(tmol, tiso, Temp)
        SigmaTref = hapi1.partitionSum(tmol, tiso, Tref)
        S[mask] = hapi1.EnvironmentDependency_Intensity(sw[mask], Temp, Tref,
                                                        SigmaT, SigmaTref,
                                                        elower[mask], nu[mask])

    # Lorentz HWHM of the diluent mixture (n_air temperature exponent for both
    # air and self, consistently with the 160-char HITRAN record)
    gammaL = pres*((Tref/Temp)**n_air)*(diluent.get('air', 0.)*gamma_air
                                        + diluent.get('self', 0.)*gamma_self)
    nu0 = nu + delta_air*pres
    height = S*gammaL/(np.pi*(gammaL*gammaL + wing*wing))

    # accumulate the constant pedestals via a difference array
    pedestal = np.zeros(len(wn_grid) + 1)
    i_lo = np.searchsorted(wn_grid, nu0 - wing, side='left')
    i_hi = np.searchsorted(wn_grid, nu0 + wing, side='right')
    np.add.at(pedestal, i_lo, height)
    np.add.at(pedestal, i_hi, -height)
    return np.cumsum(pedestal[:-1])


    
def ParallelPart(tasks, ParametersCalculation, Nwn, co_hdf5, dataset_name, METHOD, quiet=False):
    """Computes all grid points in `tasks` = [(ip,it,iv,p,T,vms), ...] and
    writes them into co_hdf5[dataset_name] as they finish.
    Returns the list of failed points (empty on full success)."""

    wn_begin = float(ParametersCalculation[3][1])
    wn_end = float(ParametersCalculation[4][1])

    molec_id = int(ParametersCalculation[10][1])

    par_group = ParametersCalculation[19][1]
    par_group = [par_group]  
    
    iso_array = [[1,2,3,4,5,6,129],                           # H2O
                 [7,8,9,10,11,12,13,14,15],#,120,121,122],       # CO2
                 [16,17,18,19,20],                            # O3
                 [21,22,23,24,25],                            # N2O
                 [26,27,28,29,30,31],                         # CO
                 [32,33,34,35],                               # CH4
                 [36,37,38],                                  # O2
                 [39,40,41],                                  # NO
                 [42,43,137,138],                             # SO2
                 [44,130],                                    # NO2
                 [45,46],                                     # NH3
                 [47,117],                                    # HNO3
                 [48,49,50],                                  # OH
                 [51,110],                                    # HF
                 [52,53,107,108],                             # HCl
                 [54,55,111,112],                             # HBr
                 [56,113],                                    # HI
                 [57,58],                                     # ClO
                 [59,60,61,62,63,135],                        # OCS
                 [64,65,66],                                  # H2CO
                 [67,68],                                     # HOCl
                 [69,118],                                    # N2
                 [70,71,72],                                  # HCN
                 [73,74],                                     # CH3Cl
                 [75],                                        # H2O2
                 [76,77,105],                                 # C2H2
                 [78,106],                                    # C2H6
                 [79],                                        # PH3
                 []]
                 # To be updated!                          

    iso_list = iso_array[molec_id-1]




    FLAG_LINE_MIXING = readSwitchByName(ParametersCalculation, 'Line_mixing')
    flags = {
        'line_mixing':     FLAG_LINE_MIXING,
        'remove_pedestal': readSwitchByName(ParametersCalculation, 'Remove_pedestal'),
        'keep_dat':        readSwitchByName(ParametersCalculation, 'Keep_dat'),
        'quiet_workers':   quiet,
    }

    LOG.info('*** CORE_CALCS ***')
    LOG.info('%s %s %s'%(molec_id, par_group, iso_list))
    LOG.info('Line mixing (1st-order Rosenkranz): %s'%('ON' if flags['line_mixing'] else 'OFF'))
    LOG.info('Pedestal removal: %s'%('ON' if flags['remove_pedestal'] else 'OFF'))
    LOG.info('Keep per-point .dat files: %s'%('ON' if flags['keep_dat'] else 'OFF'))

    tab_name = 'HITRAN2020'

    fetch_groups = list(par_group)
    if (flags['line_mixing']):
        # line-mixing parameters, where HITRAN provides them
        fetch_groups += ['voigt_linemixing', 'sdvoigt_linemixing']
    def _fetch(groups):
        if (quiet):
            with open(os.devnull, 'w') as devnull, contextlib.redirect_stdout(devnull):
                hapi1.fetch_by_ids(  tab_name, iso_list,
                                     max(0., wn_begin-FETCH_MARGIN),  wn_end+FETCH_MARGIN,
                                     ParameterGroups=groups)
        else:
            hapi1.fetch_by_ids(  tab_name, iso_list,
                                 max(0., wn_begin-FETCH_MARGIN),  wn_end+FETCH_MARGIN,
                                 ParameterGroups=groups)
    try:
        _fetch(fetch_groups)
    except Exception as err:
        if (flags['line_mixing']):
            LOG.warning('WARNING: fetch with line-mixing groups failed (%s), refetching without them'%(err))
            _fetch(par_group)
        else:
            raise

    # NOTE: do not call hapi1.cache2storage() here: fetch_by_ids() has already
    # persisted the table, and cache2storage() in HAPI <= 1.3.0.0 rewrites the
    # header dropping the 'extra' column info (e.g. the line-mixing parameters),
    # which breaks storage2cache() in the worker processes.

    if (flags['keep_dat'] and ('datafiles' not in os.listdir('./'))):
        os.mkdir('datafiles')

    ds_coef = co_hdf5[dataset_name]
    n_total = len(tasks)

    if (METHOD=='PC'):
        LOG.info('METHOD IS ASYNCIO (compatibility mode: GIL-bound, no speedup for CPU-bound work)')
        LOG.info('Number of CPUs in the system: %d'%os.cpu_count())
        _init_worker(tab_name, ParametersCalculation, Nwn, flags)
        loop = asyncio.get_event_loop()                                              # Have a new event loop
        looper = asyncio.gather(*[CalculateXsecAS(task) for task in tasks])          # Run the loop
        results = loop.run_until_complete(looper)                                    # Wait until finish
        failures = _consume_results(iter(results), ds_coef, co_hdf5, Nwn, n_total)
    elif (METHOD=='PLAIN'):
        LOG.info('METHOD IS PLAIN')
        _init_worker(tab_name, ParametersCalculation, Nwn, flags)
        failures = _consume_results(map(CalculateXsec, tasks), ds_coef, co_hdf5, Nwn, n_total)
    elif (METHOD=='MULTITHREADING'):
        class InvalidCoreCount(Exception):
            "Raised when the number of cores requested more than existed"
            pass

        try:
            LOG.info('Number of CPUs in the system: %d'%os.cpu_count())
            LOG.info('METHOD IS MULTIPROCESSING')
            N_threads = int(ParametersCalculation[17][1])
            if (N_threads>os.cpu_count()):
                raise InvalidCoreCount
            LOG.info('Number of CPUs used: %d'%N_threads)
            # each worker parses the line list once in the initializer; tasks
            # carry only plain scalars, results stream back one slice at a time
            with Pool(N_threads, initializer=_init_worker,
                      initargs=(tab_name, ParametersCalculation, Nwn, flags)) as pool:
                failures = _consume_results(
                    pool.imap_unordered(CalculateXsec, tasks, chunksize=1),
                    ds_coef, co_hdf5, Nwn, n_total)
        except InvalidCoreCount:
            LOG.error("Exception occurred: requested too many cores!")
            sys.exit()

    else:
        raise NameError('ERROR: Unknown method!')

    return failures


def save_xsc( filename, vals_nu, vals_abs  , vals_unc_l, vals_unc_u):
    xsc_w   = open( filename, 'w')    
    for i_wrt in range(0,len(vals_nu)):
        xsc_w.write( '{:>10.3f}'.format( vals_nu[i_wrt] ) )   
        xsc_w.write( '{:>15.6E}'.format( vals_abs[i_wrt] ) )
        xsc_w.write( '{:>15.6E}'.format( vals_unc_l[i_wrt] ) )
        xsc_w.write( '{:>15.6E}'.format( vals_unc_u[i_wrt] ) )
        xsc_w.write( '\n'  ) 
    xsc_w.close()


# calculate x-sec for exact P, T, VMS of exact molecule
def CalculateXsec(task):
    """task = (ip, it, iv, pres, Temp, VMS) with plain ints/floats.
    Returns (ip, it, iv, coef | None, error_text | None).
    Never raises: failures are returned, so one bad grid point cannot
    kill the whole run. Worker state comes from _init_worker()."""
    ip, it, iv, pres, Temp, VMS = task
    try:
        if (WORKER.get('init_error')):
            return (ip, it, iv, None, 'worker init failed:\n%s'%WORKER['init_error'])
        if ((pres!=pres) or (Temp!=Temp) or (VMS!=VMS)):
            return (ip, it, iv, None, 'NaN in p/T/VMS value')

        param    = WORKER['param']
        tab_name = WORKER['tab_name']
        wngrid   = WORKER['wngrid']

        wn_begin = float(param[3][1])
        wn_end = float(param[4][1])
        profile_name = param[18][1]

        FLAG_LINE_MIXING = WORKER['flags']['line_mixing']
        FLAG_REMOVE_PEDESTAL = WORKER['flags']['remove_pedestal']
        diluent = {'self':1.00-VMS, 'air':VMS}

        if (FLAG_DEBUG_PRINT):
            print('*** DEBUG: X-sec ***')
            print('VMS=%4.2f, type='%VMS, type(VMS))
            print('Range from %8.2f to %8.2f, %d points'%(wn_begin, wn_end, WORKER['Nwn']))
            print('Pressure=%6.2e, temperature=%7.2f'%(pres,Temp))
            print('*** END: X-sec ***\n')

        # WavenumberGrid pins the calculation to the exact grid stored in the
        # 'Wavenumber' dataset (HAPI's internal grid differs in the last ulp)
        kwargs = dict(SourceTables=tab_name, HITRAN_units=True,
                      OmegaRange=[wn_begin,wn_end],
                      WavenumberGrid=wngrid,
                      WavenumberWing=WING_WN, OmegaWingHW=0.0,
                      Environment={'T':Temp,'p':pres},
                      Diluent=diluent)
        if (profile_name == 'HT'):
            nu_co,coef_co = hapi1.absorptionCoefficient_HT(LineMixingRosen=False, **kwargs)
        elif (profile_name == 'SDVoigt'):
            nu_co,coef_co = hapi1.absorptionCoefficient_SDVoigt(LineMixingRosen=FLAG_LINE_MIXING, **kwargs)
        elif (profile_name == 'Voigt'):
            nu_co,coef_co = hapi1.absorptionCoefficient_Voigt(LineMixingRosen=FLAG_LINE_MIXING, **kwargs)
        elif (profile_name == 'Lorentz'):
            nu_co,coef_co = hapi1.absorptionCoefficient_Lorentz(LineMixingRosen=False, **kwargs)
        else:
            return (ip, it, iv, None, 'unknown profile name: %r'%profile_name)

        if (FLAG_REMOVE_PEDESTAL):
            coef_co = coef_co - PedestalSpectrum(nu_co, tab_name, Temp, pres, diluent, WING_WN)

        coef_co = np.asarray(coef_co, dtype=np.float64)
        if (WORKER['flags']['keep_dat']):
            zeros = np.zeros(len(nu_co))
            save_xsc( dat_filename(param, pres, Temp, VMS), nu_co, coef_co, zeros, zeros )
        return (ip, it, iv, coef_co, None)
    except Exception:
        return (ip, it, iv, None, traceback.format_exc())

@background
def CalculateXsecAS(task):
    # compatibility wrapper for the asyncio method: threads share the module
    # state set by _init_worker() in the main thread
    return CalculateXsec(task)
