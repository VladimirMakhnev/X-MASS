import numpy as np
import time
import sys
import initial
import hdf5_io
import core_calcs

# TO RECORD OUTPUT INTO THE FILE -- TRUE
# TO KEEP IT IN THE CONSOLE      -- FALSE
FLAG_LOG_FILE = False

if __name__ == "__main__":
    
    
    
    fLog = open('output.log', 'w')
    fLog.close()
    if (FLAG_LOG_FILE):
        orig_stdout = sys.stdout
        fLog = open('output.log', 'a')
        sys.stdout = fLog

    print('**************************************************************************\nXMASSSECTION: program to calculate and store cross-sections in HDF5 files.\n**************************************************************************')
    
    #############################################################
    ### BEGIN OF MAIN PART ######################################
    #############################################################
    XMASSSEC_VERSION = '0.7'; __version__ = XMASSSEC_VERSION
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
    'COMMENTS AND STYLE (ver. 0.7)'
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
    
    print("***********************************************")
    print("*** PARAMETERS HANDLING ***********************")
    print("***********************************************")
    
    
    argvlen = len(sys.argv)
    (PARAM_FILENAME := sys.argv[1]) if (argvlen>1) else (PARAM_FILENAME := 'params.inp')
    (PRES_FILENAME  := sys.argv[2]) if (argvlen>2) else (PRES_FILENAME := 'pres_pRT.inp')
    (TEMP_FILENAME  := sys.argv[3]) if (argvlen>3) else (TEMP_FILENAME := 'temps_pRT.inp')
    (VMS_FILENAME   := sys.argv[4]) if (argvlen>4) else (VMS_FILENAME := 'vms.inp')
    (WN_FILENAME    := sys.argv[5]) if (argvlen>5) else (WN_FILENAME := 'wn.inp')
    (HDF5FileName   := sys.argv[6]) if (argvlen>6) else (HDF5FileName := '01.H2O.SDV.HITRAN2020.25wing.hdf5')
    # (METHOD         := sys.argv[7]) if (argvlen>7) else (METHOD := 'MULTITHREADING')
    # (METHOD         := sys.argv[7]) if (argvlen>7) else (METHOD := 'PLAIN')
    (METHOD         := sys.argv[7]) if (argvlen>7) else (METHOD := 'PC')

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
    
    t_end = time.time()
    print('Time taken: %d seconds'%(t_end-t_begin))
    
    print('\nDone.')

    # Recording the output into the file
    if (FLAG_LOG_FILE):
        sys.stdout = orig_stdout
        fLog.close()
    
