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


    
def ParallelPart(pTVMS,WNs,ParametersCalculation,Nwn,Npp,Ntt,Nvms,co_hdf5,METHOD):
    
    wn_begin = float(ParametersCalculation[3][1])
    wn_end = float(ParametersCalculation[4][1])

    molec_id = int(ParametersCalculation[10][1])

    par_group = ParametersCalculation[19][1]
    par_group = [par_group]  
    
    iso_array = [[1,2,3,4,5,6,129],                           # H2O
                 [7,8,9,10,11,12,13,14,15,120,121,122],       # CO2
                 [16,17,18,19,20],                            # O3
                 [21,22,23,24,25],                            # N2O
                 [26,27,28,29,30,31],                         # CO
                 [32,33,34,35],                               # CH4
                 [36,37,38],                                  # O2
                 [39,40,41],                                  # NO
                 [42,43,137,138],                             # SO2
                 [44,130],                                    # NO2
                 [45,46],                                     # NH3
                 [47,117] ]                                   # HNO3
                 # To be updated!                          

    iso_list = iso_array[molec_id-1]




    print('*** CORE_CALCS ***')
    print(molec_id, par_group,iso_list)

    tab_name = 'HITRAN2020'
    
    
    hapi1.fetch_by_ids(  tab_name, iso_list,   wn_begin,  wn_end,   ParameterGroups=par_group)

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
    #        print(myargs)
            pool = Pool(N_threads)
            results = pool.map(CalculateXsec, myargs)
        except InvalidCoreCount:
            print("Exception occurred: requested too many cores!")
            sys.exit()

        
    else:
        raise NameError('ERROR: Unknown method!') 
    # for (p, T, VMS) in pTVMS:
    #     CalculateXsec(p, T,VMS, WNs,ParametersCalculation,Nwn)
     

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
# @background
def CalculateXsec(args):
    pres, Temp, VMS, WN_range, param, Nwn, tab_name = args
    print('Calculate xsec, hapitable')
    class NaNError(Exception):
        'Corrupted p, T or VMS value'
        pass
    
    try:
        if ((pres!=pres) or (Temp!=Temp) or (VMS!=VMS)):
            raise NaNError
        # print('entered?')
        
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
            print('Pressure=%6.2e, temperature=%7.2f'%(pres,Temp))
            print('*** END: X-sec ***\n')
        IndexMol = int(param[10][1])
        IndexBroad = int(param[15][1])
        
        
        # PROBLEM: FIX THE NAME AT UPD_HDF5 TOO!
        CoefFileName = './datafiles/%06.2fT_Id%02d_%06.4eatm_IdBroad%02d_%06.4fVMS_H2O_SDV_hitran2020.dat'%(Temp,IndexMol,pres,IndexBroad,VMS)
        # print("named?")
        # nu_co,coef_co = absorptionCoefficient_Voigt(SourceTables='05_HITEMP2019-full',#'07_all_iso_hit20_0k-35k',#'05_HITEMP2019',#'CO_all_HITRAN',#
        #                                               HITRAN_units=True, OmegaRange=[wn_begin,wn_end],
        #                                               WavenumberStep=wn_step,
        #                                               WavenumberWing=25.0,
        #                                               OmegaWingHW=0.0,
        #                                               Diluent={'self':1.00-VMS, 'air':VMS},
        #                                               Environment={'T':Temp,'p':pres},
        #                                               File = CoefFileName)
        # tab_name = '03_HITRAN2020'
        # hapi1.fetch_by_ids(  tab_name, [16,17,18,19,20],   wn_begin,  wn_end,   ParameterGroups=['160-char'])
        # hapi1.LOCAL_TABLE_CACHE = hapi1.storage2cache(tab_name)
        hapi1.storage2cache(tab_name)
        # print(hapi1.LOCAL_TABLE_CACHE.keys())
        # nu_co,coef_co,xunc_l, xunc_u = hapi2.opacity.lbl.numba.absorptionCoefficient_Voigt(SourceTables='05_SDV_HITRAN2020',
        #                                                                     OmegaStep = wn_step,
        #                                                                     OmegaWing=25.0,
        #                                                                     OmegaWingHW=0.0,
        #                                                                     WavenumberGrid=wngrid,
        #                                                                     Diluent={'self':1.00-VMS, 'air':VMS},
        #                                                                     Environment={'T':Temp,'p':pres},table_itself = hapitable)
        nu_co,coef_co = hapi1.absorptionCoefficient_SDVoigt(SourceTables='HITRAN2020', HITRAN_units=True,
                                                            OmegaRange=[wn_begin,wn_end],WavenumberStep=wn_step,  
                                                            WavenumberWing=25.,OmegaWingHW=0.0,LineMixingRosen=False,
                                                            Environment={'T':Temp,'p':pres},
                                                            Diluent={'self':1.00-VMS, 'air':VMS},
                                                            File = CoefFileName)
                                                            # File = CoefFileName,hapitab=hapitable)
                                                            # Components=[(5,1,(1-VMS)),(5,2,(1-VMS)),(5,3,(1-VMS)),(5,4,(1-VMS)),(5,5,(1-VMS)),(5,6,(1-VMS)),],


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
    
    print('Calculate xsec, hapitable')
    class NaNError(Exception):
        'Corrupted p, T or VMS value'
        pass
    
    try:
        if ((pres!=pres) or (Temp!=Temp) or (VMS!=VMS)):
            raise NaNError
        # print('entered?')
        
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
        
        
        # PROBLEM: FIX THE NAME AT UPD_HDF5 TOO!
        CoefFileName = './datafiles/%06.2fT_Id%02d_%06.4eatm_IdBroad%02d_%06.4fVMS_H2O_SDV_hitran2020.dat'%(Temp,IndexMol,pres,IndexBroad,VMS)
        # print("named?")
        # nu_co,coef_co = absorptionCoefficient_Voigt(SourceTables='05_HITEMP2019-full',#'07_all_iso_hit20_0k-35k',#'05_HITEMP2019',#'CO_all_HITRAN',#
        #                                               HITRAN_units=True, OmegaRange=[wn_begin,wn_end],
        #                                               WavenumberStep=wn_step,
        #                                               WavenumberWing=25.0,
        #                                               OmegaWingHW=0.0,
        #                                               Diluent={'self':1.00-VMS, 'air':VMS},
        #                                               Environment={'T':Temp,'p':pres},
        #                                               File = CoefFileName)
        # tab_name = '03_HITRAN2020'
        # hapi1.fetch_by_ids(  tab_name, [16,17,18,19,20],   wn_begin,  wn_end,   ParameterGroups=['160-char'])
        hapi1.LOCAL_TABLE_CACHE = hapitable
        print(hapitable.keys())
        # nu_co,coef_co,xunc_l, xunc_u = hapi2.opacity.lbl.numba.absorptionCoefficient_Voigt(SourceTables='05_SDV_HITRAN2020',
        #                                                                     OmegaStep = wn_step,
        #                                                                     OmegaWing=25.0,
        #                                                                     OmegaWingHW=0.0,
        #                                                                     WavenumberGrid=wngrid,
        #                                                                     Diluent={'self':1.00-VMS, 'air':VMS},
        #                                                                     Environment={'T':Temp,'p':pres},table_itself = hapitable)
        nu_co,coef_co = hapi1.absorptionCoefficient_SDVoigt(SourceTables='HITRAN2020', HITRAN_units=True,
                                                            OmegaRange=[wn_begin,wn_end],WavenumberStep=wn_step,  
                                                            WavenumberWing=25.,OmegaWingHW=0.0,LineMixingRosen=False,
                                                            Environment={'T':Temp,'p':pres},
                                                            Diluent={'self':1.00-VMS, 'air':VMS},
                                                            File = CoefFileName)
                                                            # Components=[(5,1,(1-VMS)),(5,2,(1-VMS)),(5,3,(1-VMS)),(5,4,(1-VMS)),(5,5,(1-VMS)),(5,6,(1-VMS)),],

        

                                                      #'05_HITEMP2019-full',#'07_all_iso_hit20_0k-35k',#'05_HITEMP2019',#'CO_all_HITRAN',#
                                                      #HITRAN_units=True, OmegaRange=[wn_begin,wn_end],
                                                      #WavenumberStep=wn_step,
                                                      #WavenumberWing=25.0,
                                                      #OmegaWingHW=0.0,
                                                      #Diluent={'self':1.00-VMS, 'air':VMS},
                                                      #Environment={'T':Temp,'p':pres})
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
        
