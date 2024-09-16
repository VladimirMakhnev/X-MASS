import numpy as np
import sys

# FLAGS
FLAG_DEBUG_PRINT = True




# opens file with necessary parameters for calculations
def openParametersFile(fname):
    class ParamsError(Exception):
        pass
    try:
        global FLAG_DEBUG_PRINT
        with open(fname,'r') as finp:
            # input as list of lines in 'lines'
            lines = finp.readlines()
            # check if input has enough parameters
            if (len(lines)<17):
                raise ParamsError
            # removing '\n'
            lines = [item.rstrip('\n') for item in lines]
            # splitting names and values
            params = [item.split(':') for item in lines]
            if (FLAG_DEBUG_PRINT):
                print('*** DEBUG: Parameters input ***')
                [print(item[0],'\t',item[1]) for item in params]
                print('*** END: Parameters input ***\n')
        print('*********\nParameters file %s is opened well\n*********'%(fname))
        return params
    except FileNotFoundError:
        print('%s file is not found!'%fname)
        sys.exit()
        # raise FileNotFoundError
    except ParamsError:
        print("Not enought parameters in %s!"%fname)
        sys.exit()
    except Exception as err:
        print('WOW, UNKNOWN ERROR: %s'%(err))
        sys.exit()



# opens pressure array file
def openPressure(fname):
    try:
        apres = np.genfromtxt(fname,dtype='float')
        apres = apres#*1.0e-3
        global FLAG_DEBUG_PRINT
        if (FLAG_DEBUG_PRINT):
            print('*** DEBUG: Pressure input ***')
            [print('%16.9f'%(item)) for item in apres]
            print('*** END: Pressure input ***\n')
        print('*********\nPressure file %s is opened well\nTotal number of lines: %d\n*********'%(fname,len(apres)))
        return apres, len(apres)
    except FileNotFoundError:
        print('%s file is not found!'%fname)
        sys.exit()
    
# opens temperature array file
def openTemp(fname,Np):
    class PxTError(Exception):
        'Corrupted relations between Np and NpxNt array'
        pass
    try:
        atemp = np.genfromtxt(fname,dtype='float', missing_values='296.15')
        (Npp, Ntt) = atemp.shape
        global FLAG_DEBUG_PRINT
        if (Npp!=Np):
            raise PxTError
        if (FLAG_DEBUG_PRINT):
            print('*** DEBUG: Temperature input ***')
            for item in atemp:
                [print(item1, end='\t') for item1 in item]
                print('')
            print('*** END: Temperature input ***\n')
        print('*********\nTemperature file %s is opened well\nTotal number of lines: %d\nNo of temperatures: %d\n*********'%(fname,Npp,Ntt))
        return atemp, Npp, Ntt
    except FileNotFoundError:
        print('%s file is not found!'%fname)
        sys.exit()
    except PxTError:
        print('Corrupted relations between Np and NpxNt array')
        sys.exit()

# opens pressure array file
def openVMS(fname):
    try:
        avms = np.genfromtxt(fname,dtype='float')
        global FLAG_DEBUG_PRINT
        if (FLAG_DEBUG_PRINT):
            print('*** DEBUG: VMS input ***')
            if (avms.shape != ()):
                [print('%12.6f'%(item)) for item in avms]
            else:
                print(avms)
            print('*** END: VMS input ***\n')
        if (len(avms)==1):
            avms = [avms]
        # avms = [avms]
        print('*********\nVMS file %s is opened well\nTotal number of lines: %d\n*********'%(fname,len(avms)))
        return avms, len(avms)
    except FileNotFoundError:
        print('%s file is not found!'%fname)
        sys.exit()

# opens and generate wn_grid
def openXgenetareWn(fname,params):
    try:
        global FLAG_DEBUG_PRINT
        with open(fname) as f:
            Nwn = int(f.readline())
            wn_begin = float(params[3][1])
            wn_end = float(params[4][1])
            WN_range = np.linspace(wn_begin,wn_end,Nwn)
#             print(WN_range[:15])
            if (FLAG_DEBUG_PRINT):
                print('*** DEBUG: WN input ***')
                [print('%12.6f'%item1) for item1 in WN_range[:10]]
                print('...')
                [print('%12.6f'%item2) for item2 in WN_range[-10:]]
                print(Nwn)
                print('*** END: WN input ***\n')
            print('*********\nWN file %s is opened well\nTotal number of wn-points: %d\n*********'%(fname,Nwn))
                
            return WN_range, Nwn
    except FileNotFoundError:
        print('%s file is not found!'%fname)
        sys.exit()

# collects the pressures/temperatures/volume mixing ration into the pTVMS array
def mergeParams(P,T,VMS):
    onevisionlist = []
    indexvisionlist = []
    for iv, tv in np.ndenumerate(VMS):
        for ip, tp in np.ndenumerate(P):
            for it, tt in np.ndenumerate(T[ip]):
                onevisionlist.append([tp,tt,tv])
                indexvisionlist.append([ip,it,iv])
    return onevisionlist, indexvisionlist
    
    





















