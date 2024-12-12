import numpy as np
import time
import sys
import initial
import hdf5_io
import core_calcs

import winsound

# TO RECORD OUTPUT INTO THE FILE -- TRUE
# TO KEEP IT IN THE CONSOLE      -- FALSE
FLAG_LOG_FILE = True
FLAG_PROFILER = True

if __name__ == "__main__":
    
    
    
    fLog = open('output.log', 'w')
    fLog.close()
    if (FLAG_LOG_FILE):
        import cProfile
        import pstats

        orig_stdout = sys.stdout
        fLog = open('output.log', 'a')
        sys.stdout = fLog

    print('**************************************************************************\nXMASSSECTION: program to calculate and store cross-sections in HDF5 files.\n**************************************************************************')
    
    #############################################################
    ### BEGIN OF MAIN PART ######################################
    #############################################################
    XMASSSEC_VERSION = '0.8.1'; __version__ = XMASSSEC_VERSION
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
    'TIME PROFILER (ver. 0.8.1)'
    ]
    
    # version header
    print('X-MASS-SEC version: %s'%(XMASSSEC_VERSION))
#    print('To get the most up-to-date version please check http://hitran.org/')
    # print('ATTENTION: Python versions of partition sums from TIPS-2021 are now available in HAPI code')
    #print('ATTENTION: Python versions of partition sums from TIPS-2017 are available at http://hitran.org/suppl/TIPS/')
    #print('           To use them in HAPI ver. 1.1.0.7, use partitionFunction parameter of the absorptionCoefficient_ routine.')
    print('')
    print('           MIT license: Copyright 2024 HITRAN team, see more at http://hitran.org. ')
    print('')
    print('           If you use X-MASS in your research or software development,')
    print('           please cite it using the following reference:')
    print('           V.Yu. Makhnev, I.E. Gordon, L.S. Rothman, R.J. Hargreaves')
    print('           ')
    print('           ')
    print('           DOI: ')
    print('')
    # print('           ATTENTION: This is the core version of the HITRAN Application Programming Interface.')
    # print('                      For more efficient implementation of the absorption coefficient routine, ')
    # print('                      as well as for new profiles, parameters and other functional,')
    # print('                      please consider using HAPI2 extension library.')
    # print('                      HAPI2 package is available at http://github.com/hitranonline/hapi2')
    print('')


    #####################################################
    
    
    #####################################################
    
    print("Timer started")
    t_begin = time.time()
    profiler = cProfile.Profile()
    if (FLAG_PROFILER):
        profiler.enable()
    
    print("***********************************************")
    print("*** PARAMETERS HANDLING ***********************")
    print("***********************************************")
    
    
    argvlen = ((open('filenames.inp',mode='r')).read()).split()
    print(argvlen)
    (PARAM_FILENAME := argvlen[1]) if (len(argvlen)>1) else (PARAM_FILENAME := 'params.inp')
    (PRES_FILENAME  := argvlen[2]) if (len(argvlen)>2) else (PRES_FILENAME := 'pres_pRT.inp')
    (TEMP_FILENAME  := argvlen[3]) if (len(argvlen)>3) else (TEMP_FILENAME := 'temps_pRT.inp')
    (VMS_FILENAME   := argvlen[4]) if (len(argvlen)>4) else (VMS_FILENAME := 'vms.inp')
    (WN_FILENAME    := argvlen[5]) if (len(argvlen)>5) else (WN_FILENAME := 'wn.inp')
    (HDF5FileName   := argvlen[6]) if (len(argvlen)>6) else (HDF5FileName := '01.H2O.SDV.HITRAN2020.25wing.hdf5')
    (METHOD         := argvlen[7]) if (len(argvlen)>7) else (METHOD := 'MULTITHREADING')
    # (METHOD         := sys.argv[7]) if (argvlen>7) else (METHOD := 'PLAIN')
    # (METHOD         := sys.argv[7]) if (argvlen>7) else (METHOD := 'PC')

    # Opening parameter file
    ParametersCalculation = initial.openParametersFile(PARAM_FILENAME)
    print(ParametersCalculation)
    
    # opening pressure file
    Pressures, Np = initial.openPressure(PRES_FILENAME)
    print(Pressures,'\n', Np)
    
    # Opening temperature file 
    (Temps, Npp,Ntt) = initial.openTemp(TEMP_FILENAME,Np)
    print(Temps, Npp,Ntt)
    
    # Opening volume mixing ratio file
    (VMSs, Nvms) = initial.openVMS(VMS_FILENAME)
    print(VMSs, Nvms)
    
    # Opening wavenumber range file
    (WNs, Nwn) = initial.openXgenetareWn(WN_FILENAME,ParametersCalculation)
    print(WNs, Nwn)
    
    print("***********************************************")
    print("*** OPENING HDF5 FILE *************************")
    print("***********************************************")
    
    # Initialazing the core HDF5 file
    co_hdf5 = hdf5_io.OpenHDF5(HDF5FileName, ParametersCalculation, Pressures, Temps, VMSs, WNs, Npp, Ntt, Nvms, Nwn)

    # Constructing p/T/VMS array out of values
    pTVMS, ipTVMS = initial.mergeParams(Pressures, Temps, VMSs)
    
    # print(len(pTVMS))

    print("***********************************************")
    print("*** CALCULATIONS PART *************************")
    print("***********************************************")
    
    # Calculation part
    co_hdf5 = core_calcs.ParallelPart(pTVMS,WNs,ParametersCalculation,Nwn,Npp,Ntt,Nvms,co_hdf5,METHOD)
    # Populating the HDF5 file by x-sections
    co_hdf5 = hdf5_io.UpdateHDF5(co_hdf5, pTVMS, ipTVMS, ParametersCalculation)
  
    
    
    
    
    
    
    
    
    
    
    
    
    
    # Closing the HDF5 file
    hdf5_io.CloseHDF5(co_hdf5)
    
    if (FLAG_PROFILER): 
        profiler.disable()  
        # Save profiling results to a file
        profiler.dump_stats("profiling_results.prof")
    
        # Load profiling results
        stats = pstats.Stats("profiling_results.prof")
    
        # Sort by cumulative time (default)
        print("Sorted by cumulative time:")
        stats.sort_stats("cumulative").print_stats()

    t_end = time.time()
    print('Time taken: %d seconds'%(t_end-t_begin))
    winsound.Beep(261, 400)
    winsound.Beep(329, 400)
    winsound.Beep(392 , 400)
    winsound.Beep(523, 700)
    print('\nDone.')

    # Recording the output into the file
    if (FLAG_LOG_FILE):
        sys.stdout = orig_stdout
        fLog.close()
    
