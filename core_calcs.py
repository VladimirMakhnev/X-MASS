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
    # tab_name = '05_HITEMP2019-full'
    # hapi1.fetch('05_HITEMP2019-full', 05, I, numin, numax)
    # hapi1.db_begin('HITRAN')




    tab_name = 'H2O_SDV_HITRAN2020'
#    hapi1.fetch_by_ids(  tab_name, [26],   wn_begin,  wn_end,   ParameterGroups=['160-char',"SDvoigt"],
#                        Parameters=['n_SDV_self_296','n_SDV_air_296','n_gamma_SDV_2_self_296','n_gamma_SDV_2_air_296'])
    hapi1.fetch_by_ids(  tab_name, [1],   wn_begin,  wn_end,   ParameterGroups=['160-char'])

    # tab_name = '51_SDV_HITRAN2020'
    # hapi1.db_begin('HITRAN')
    # hapi1.getHelp(hapi1.fetch)

    hapitable = hapi1.LOCAL_TABLE_CACHE

    if (METHOD=='PC'):
        print('METHOD IS ASYNCIO')
        print('Number of CPUs in the system: {}'.format(os.cpu_count()))

        loop = asyncio.get_event_loop()                                              # Have a new event loop
        looper = asyncio.gather(*[CalculateXsecAS(p, T,VMS, WNs,ParametersCalculation,Nwn,hapitable) for (p, T, VMS) in pTVMS])         # Run the loop
        results = loop.run_until_complete(looper)                                    # Wait until finish
    elif (METHOD=='PLAIN'):
        print('METHOD IS PLAIN')
        for (p, T, VMS) in pTVMS:
            t_myarg = (p, T, VMS, WNs, ParametersCalculation, Nwn, hapitable)
            CalculateXsec(t_myarg) 
    elif (METHOD=='MULTITHREADING'):
        print('Number of CPUs in the system: {}'.format(os.cpu_count()))
        print('METHOD IS MULTIPROCESSING')
        N_threads = 7
        Nptvms = len(pTVMS)
        myargs = []
        for iptvms in np.arange(Nptvms):
            p, T, VMS = pTVMS[iptvms]
            t_myarg = (p, T, VMS, WNs, ParametersCalculation, Nwn, hapitable)
            myargs.append(t_myarg)
#        print(myargs)
        pool = Pool(N_threads)
        results = pool.map(CalculateXsec, myargs)
### WORKING PART
#        pros = np.empty(len(pTVMS),dtype=type(Process))
#        for iptvms in np.arange(len(pTVMS)):
#            (p,T,VMS) = pTVMS[iptvms]
#            print(p,T,VMS)
#            print(iptvms)
#            pros[iptvms] = Process(target=CalculateXsec, args=(p, T,VMS, WNs,ParametersCalculation,Nwn,hapitable))
#            pros[iptvms].start()
#            
#        for i in np.arange(len(pTVMS)):
#            pros[i].join() 
###
#        Nptvms = len(pTVMS)
#        Nucore = 15
#        pros = np.empty(Nucore,dtype=type(Process))
#        for n_current in np.arange(Nptvms//Nucore):
#            pros = np.empty(Nucore,dtype=type(Process))
#
#            for i_itern in np.arange(Nucore):
#                (p,T,VMS) = pTVMS[i_itern + n_current*Nucore]
#                print(p,T,VMS)
#                print(i_itern + n_current*Nucore)
#                pros[i_itern] = Process(target=CalculateXsec, args=(p, T,VMS, WNs,ParametersCalculation,Nwn,hapitable))
#                pros[i_itern].start()
#            for i in np.arange(Nucore):
#                pros[i].join() 
#        pros = np.empty(Nucore,dtype=type(Process))
#        for i_rest in np.arange(Nucore*(Nptvms//Nucore),Nptvms):
#            (p,T,VMS) = pTVMS[rest]
#            print(p,T,VMS)
#            print(i_rest)
#            i_resty = i_rest - Nucore*Nptvms//Nucore
#            pros[i_resty] = Process(target=CalculateXsec, args=(p, T,VMS, WNs,ParametersCalculation,Nwn,hapitable))
#            pros[i_resty].start()
#        for i_rest in np.arange(-Nucore*(Nptvms//Nucore)+Nptvms):
#            pros[i_rest].join()

        
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
    pres, Temp, VMS, WN_range, param, Nwn, hapitable = args
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
        nu_co,coef_co = hapi1.absorptionCoefficient_SDVoigt(SourceTables='H2O_SDV_HITRAN2020', HITRAN_units=True,
                                                            OmegaRange=[wn_begin,wn_end],WavenumberStep=wn_step,  
                                                            WavenumberWing=25.,OmegaWingHW=0.0,LineMixingRosen=False,
                                                            Environment={'T':Temp,'p':pres},
                                                            Diluent={'self':1.00-VMS, 'air':VMS},
                                                            File = CoefFileName,hapitab=hapitable)
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
        nu_co,coef_co = hapi1.absorptionCoefficient_SDVoigt(SourceTables='H2O_SDV_HITRAN2020', HITRAN_units=True,
                                                            OmegaRange=[wn_begin,wn_end],WavenumberStep=wn_step,  
                                                            WavenumberWing=25.,OmegaWingHW=0.0,LineMixingRosen=False,
                                                            Environment={'T':Temp,'p':pres},
                                                            Diluent={'self':1.00-VMS, 'air':VMS},
                                                            File = CoefFileName,hapitab=hapitable)
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
        
