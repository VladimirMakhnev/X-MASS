from initial import FLAG_DEBUG_PRINT

import json
from getpass import getpass

import numpy as np
import pylab as pl
import matplotlib as mpl

#from hapi2 import *
#from hapi import *
import hapi as hapi1
# import hapi2

import sys
import os


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

# reads an optional ON/OFF switch from the params.inp lines (OFF when absent)
def readSwitch(param, index):
    return (len(param) > index) and (len(param[index]) > 1) and \
           (param[index][1].strip().upper() in ('ON', 'TRUE', 'YES', '1'))

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


    
def ParallelPart(pTVMS,WNs,ParametersCalculation,Nwn,Npp,Ntt,Nvms,co_hdf5,METHOD):
    
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




    FLAG_LINE_MIXING = readSwitch(ParametersCalculation, 20)

    print('*** CORE_CALCS ***')
    print(molec_id, par_group,iso_list)
    print('Line mixing (1st-order Rosenkranz): %s'%('ON' if FLAG_LINE_MIXING else 'OFF'))

    tab_name = 'HITRAN2020'

    fetch_groups = list(par_group)
    if (FLAG_LINE_MIXING):
        # line-mixing parameters, where HITRAN provides them
        fetch_groups += ['voigt_linemixing', 'sdvoigt_linemixing']
    try:
        hapi1.fetch_by_ids(  tab_name, iso_list,   wn_begin,  wn_end,   ParameterGroups=fetch_groups)
    except Exception as err:
        if (FLAG_LINE_MIXING):
            print('WARNING: fetch with line-mixing groups failed (%s), refetching without them'%(err))
            hapi1.fetch_by_ids(  tab_name, iso_list,   wn_begin,  wn_end,   ParameterGroups=par_group)
        else:
            raise

    hapitable = hapi1.LOCAL_TABLE_CACHE
    
    hapi1.cache2storage(tab_name)
    
    if ('datafiles' not in os.listdir('./')):
        os.mkdir('datafiles')
    
    
    if (METHOD=='PC'):
        print('METHOD IS ASYNCIO')
        print('Number of CPUs in the system: {}'.format(os.cpu_count()))

        loop = asyncio.get_event_loop()                                              # Have a new event loop
        looper = asyncio.gather(*[CalculateXsecAS(p, T,VMS, WNs,ParametersCalculation,Nwn,hapitable) for (p, T, VMS) in pTVMS])         # Run the loop
        results = loop.run_until_complete(looper)                                    # Wait until finish
    elif (METHOD=='PLAIN'):
        print('METHOD IS PLAIN')
        for (p, T, VMS) in pTVMS:
            t_myarg = (p, T, VMS, WNs, ParametersCalculation, Nwn, tab_name)
            CalculateXsec(t_myarg) 
    elif (METHOD=='MULTITHREADING'):
        class InvalidCoreCount(Exception):
            "Raised when the number of cores requested more than existed"
            pass

        try:
            print('Number of CPUs in the system: {}'.format(os.cpu_count()))
            print('METHOD IS MULTIPROCESSING')
            N_threads = int(ParametersCalculation[17][1])
            if (N_threads>os.cpu_count()):
                raise InvalidCoreCount
            print('Number of CPUs used: %d'%N_threads)
            Nptvms = len(pTVMS)
            myargs = []
            for iptvms in np.arange(Nptvms):
                p, T, VMS = pTVMS[iptvms]
                # t_myarg = (p, T, VMS, WNs, ParametersCalculation, Nwn, hapitable)
                t_myarg = (p, T, VMS, WNs, ParametersCalculation, Nwn, tab_name)

                myargs.append(t_myarg)
            # for item in myargs: 
            #     print('   ',item)
            pool = Pool(N_threads)
            results = pool.map(CalculateXsec, myargs)
        except InvalidCoreCount:
            print("Exception occurred: requested too many cores!")
            sys.exit()

        
    else:
        raise NameError('ERROR: Unknown method!') 
     

    return co_hdf5


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
def CalculateXsec(args):
    pres, Temp, VMS, WN_range, param, Nwn, tab_name = args
    # print('Calculate xsec, hapitable')
    class NaNError(Exception):
        'Corrupted p, T or VMS value'
        pass
    class ProfileNameError(Exception):
        'Corrupted profile name'
        pass
    
    try:
        if ((pres!=pres) or (Temp!=Temp) or (VMS!=VMS)):
            raise NaNError
        
        wn_begin = float(param[3][1])
        wn_end = float(param[4][1])
    
        wn_step = (wn_end-wn_begin)/(Nwn-1)
        
        profile_name = param[18][1]
        # print(profile_name)

        FLAG_LINE_MIXING = readSwitch(param, 20)
        FLAG_REMOVE_PEDESTAL = readSwitch(param, 21)
        diluent = {'self':1.00-VMS, 'air':VMS}

        # print('%25.22f'%wn_step)
        # print('core_calcs: wn_len', param[3][1],param[4][1])
        # print('core_calcs: wn_step',wn_step)

        wngrid =  np.linspace(wn_begin,wn_end,Nwn)

        if (FLAG_DEBUG_PRINT):
            print('*** DEBUG: X-sec ***')
            print('VMS=%4.2f, type='%VMS, type(VMS))
    #        print(tableList())
            print('Range from %8.2f to %8.2f, step %6.2f'%(wn_begin, wn_end,wn_step))
            print('Pressure=%6.2e, temperature=%7.2f'%(pres,Temp))
            print('*** END: X-sec ***\n')
        IndexMol = int(param[10][1])
        IndexBroad = int(param[15][1])
        
        
        # PROBLEM: FIX THE NAME AT UPD_HDF5 TOO!
        CoefFileName = './datafiles/%06.2fT_Id%02d_%06.4eatm_IdBroad%02d_%06.4fVMS_H2O_SDV_hitran2020.dat'%(Temp,IndexMol,pres,IndexBroad,VMS)
        hapi1.storage2cache(tab_name)
        # print(hapi1.LOCAL_TABLE_CACHE.keys())
        if (profile_name == 'HT'):
            nu_co,coef_co = hapi1.absorptionCoefficient_HT(SourceTables='HITRAN2020', HITRAN_units=True,
                                                                OmegaRange=[wn_begin,wn_end],WavenumberStep=wn_step,
                                                                WavenumberWing=WING_WN,OmegaWingHW=0.0,LineMixingRosen=False,
                                                                Environment={'T':Temp,'p':pres},
                                                                Diluent=diluent,
                                                                File = CoefFileName)

        elif (profile_name == 'SDVoigt'):
            nu_co,coef_co = hapi1.absorptionCoefficient_SDVoigt(SourceTables='HITRAN2020', HITRAN_units=True,
                                                                OmegaRange=[wn_begin,wn_end],WavenumberStep=wn_step,
                                                                WavenumberWing=WING_WN,OmegaWingHW=0.0,LineMixingRosen=FLAG_LINE_MIXING,
                                                                Environment={'T':Temp,'p':pres},
                                                                Diluent=diluent,
                                                                File = CoefFileName)
        elif (profile_name == 'Voigt'):
            nu_co,coef_co = hapi1.absorptionCoefficient_Voigt(SourceTables='HITRAN2020', HITRAN_units=True,
                                                                OmegaRange=[wn_begin,wn_end],WavenumberStep=wn_step,
                                                                WavenumberWing=WING_WN,OmegaWingHW=0.0,LineMixingRosen=FLAG_LINE_MIXING,
                                                                Environment={'T':Temp,'p':pres},
                                                                Diluent=diluent,
                                                                File = CoefFileName)
        elif (profile_name == 'Lorentz'):
            nu_co,coef_co = hapi1.absorptionCoefficient_Lorentz(SourceTables='HITRAN2020', HITRAN_units=True,
                                                                OmegaRange=[wn_begin,wn_end],WavenumberStep=wn_step,
                                                                WavenumberWing=WING_WN,OmegaWingHW=0.0,LineMixingRosen=False,
                                                                Environment={'T':Temp,'p':pres},
                                                                Diluent=diluent,
                                                                File = CoefFileName)
        else:
            raise ProfileNameError

        if (FLAG_REMOVE_PEDESTAL):
            coef_co = coef_co - PedestalSpectrum(nu_co, tab_name, Temp, pres, diluent, WING_WN)

        xunc_l, xunc_u = np.zeros(len(nu_co)),np.zeros(len(nu_co))
        save_xsc( CoefFileName, nu_co, coef_co, xunc_l, xunc_u  )
        # print('saved?')
        return
    except NaNError:
        print('Corrupted p, T or VMS value')
        return np.linspace(-1.0,-1.0,Nwn)
    else:
        err = Exception
        print("Unexpected %s"%(err))
        sys.exit()
    
@background
def CalculateXsecAS(pres, Temp, VMS,WN_range, param, Nwn, hapitable):
    
    # print('Calculate xsec, hapitable')
    class NaNError(Exception):
        'Corrupted p, T or VMS value'
        pass
    class ProfileNameError(Exception):
        'Corrupted profile name'
        pass
    
    try:
        if ((pres!=pres) or (Temp!=Temp) or (VMS!=VMS)):
            raise NaNError
        
        wn_begin = float(param[3][1])
        wn_end = float(param[4][1])
    
        wn_step = (wn_end-wn_begin)/(Nwn-1)
        print('%25.22f'%wn_step)
        print('core_calcs: wn_len', param[3][1],param[4][1])
        print('core_calcs: wn_step',wn_step)
        wngrid =  np.linspace(wn_begin,wn_end,Nwn)
        if (FLAG_DEBUG_PRINT):
            print('*** DEBUG: X-sec ***')
            print('VMS=%4.2f, type='%VMS, type(VMS))
    #        print(tableList())
            print('Range from %8.2f to %8.2f, step %6.2f'%(wn_begin, wn_end,wn_step))
            print('Pressure=%6.2f, temperature=%7.2f'%(pres,Temp))
            print('*** END: X-sec ***\n')
        IndexMol = int(param[10][1])
        IndexBroad = int(param[15][1])
        
        
        FLAG_LINE_MIXING = readSwitch(param, 20)
        FLAG_REMOVE_PEDESTAL = readSwitch(param, 21)
        diluent = {'self':1.00-VMS, 'air':VMS}

        # PROBLEM: FIX THE NAME AT UPD_HDF5 TOO!
        CoefFileName = './datafiles/%06.2fT_Id%02d_%06.4eatm_IdBroad%02d_%06.4fVMS_H2O_SDV_hitran2020.dat'%(Temp,IndexMol,pres,IndexBroad,VMS)
        hapi1.LOCAL_TABLE_CACHE = hapitable
        # print(hapitable.keys())
        nu_co,coef_co = hapi1.absorptionCoefficient_SDVoigt(SourceTables='HITRAN2020', HITRAN_units=True,
                                                            OmegaRange=[wn_begin,wn_end],WavenumberStep=wn_step,
                                                            WavenumberWing=WING_WN,OmegaWingHW=0.0,LineMixingRosen=FLAG_LINE_MIXING,
                                                            Environment={'T':Temp,'p':pres},
                                                            Diluent=diluent,
                                                            File = CoefFileName)
        if (FLAG_REMOVE_PEDESTAL):
            coef_co = coef_co - PedestalSpectrum(nu_co, 'HITRAN2020', Temp, pres, diluent, WING_WN)
        xunc_l, xunc_u = np.zeros(len(nu_co)),np.zeros(len(nu_co))
        save_xsc( CoefFileName, nu_co, coef_co, xunc_l, xunc_u  )
        # print('saved?')
        return
    except NaNError:
        print('Corrupted p, T or VMS value')
        return np.linspace(-1.0,-1.0,Nwn)
    else:
        err = Exception
        print("Unexpected %s"%(err))
        sys.exit()
        
