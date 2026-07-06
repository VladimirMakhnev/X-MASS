import numpy as np
import sys
import logging

LOG = logging.getLogger('xmass')

# FLAGS
FLAG_DEBUG_PRINT = False




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
        LOG.info('Parameters file %s is opened well'%(fname))
        return params
    except FileNotFoundError:
        LOG.error('%s file is not found!'%fname)
        sys.exit(2)
        # raise FileNotFoundError
    except ParamsError:
        LOG.error("Not enough parameters in %s!"%fname)
        sys.exit(2)
    except Exception as err:
        LOG.error('UNKNOWN ERROR reading %s: %s'%(fname,err))
        sys.exit(2)



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
        LOG.info('Pressure file %s is opened well, total number of lines: %d'%(fname,len(apres)))
        return apres, len(apres)
    except FileNotFoundError:
        LOG.error('%s file is not found!'%fname)
        sys.exit(2)
    
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
        LOG.info('Temperature file %s is opened well, %d pressure rows x %d temperatures'%(fname,Npp,Ntt))
        return atemp, Npp, Ntt
    except FileNotFoundError:
        LOG.error('%s file is not found!'%fname)
        sys.exit(2)
    except PxTError:
        LOG.error('Temperature grid shape does not match the number of pressures (expected %d rows)'%Np)
        sys.exit(2)

# opens pressure array file
def openVMS(fname):
    try:
        # atleast_1d: a single-value file yields a 0-d array, which breaks len()
        avms = np.atleast_1d(np.genfromtxt(fname,dtype='float'))
        global FLAG_DEBUG_PRINT
        if (FLAG_DEBUG_PRINT):
            print('*** DEBUG: VMS input ***')
            [print('%12.6f'%(item)) for item in avms]
            print('*** END: VMS input ***\n')
        LOG.info('VMS file %s is opened well, total number of lines: %d'%(fname,len(avms)))
        return avms, len(avms)
    except FileNotFoundError:
        LOG.error('%s file is not found!'%fname)
        sys.exit(2)

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
            LOG.info('WN file %s is opened well, total number of wn-points: %d'%(fname,Nwn))
                
            return WN_range, Nwn
    except FileNotFoundError:
        LOG.error('%s file is not found!'%fname)
        sys.exit(2)

# name-based lookup for option keys; use only with keys that appear at most
# once in the parameter file (legacy duplicated keys stay positional)
def readParamByName(params, name, default=None):
    for item in params:
        if (item[0].strip() == name) and (len(item) > 1):
            return item[1].strip()
    return default

def readSwitchByName(params, name, default=False):
    val = readParamByName(params, name)
    if (val is None):
        return default
    return val.upper() in ('ON', 'TRUE', 'YES', '1')

# override or append an option value (used by the command-line overrides)
def setParamByName(params, name, value):
    for item in params:
        if (item[0].strip() == name):
            if (len(item) > 1):
                item[1] = value
            else:
                item.append(value)
            return
    params.append([name, value])

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

# flat task list for the parallel part: one (ip, it, iv, p, T, vms) tuple of
# plain ints/floats per grid point (cheap to pickle, indices travel with the
# task and come back with the result)
def mergeParamsIndexed(P, T, VMS):
    tasks = []
    for iv, tv in np.ndenumerate(VMS):
        for ip, tp in np.ndenumerate(P):
            for it, tt in np.ndenumerate(T[ip]):
                tasks.append((ip[0], it[0], iv[0], float(tp), float(tt), float(tv)))
    return tasks
    
    





















