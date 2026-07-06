import h5py
import sys
import numpy as np
from initial import FLAG_DEBUG_PRINT, readParamByName, readSwitchByName

# spectral chunk length: 2 MB of float64 per chunk; with chunks=(1,1,1,C) each
# grid point owns whole chunks, so per-point writes need no read-modify-write
HDF5_WN_CHUNK = 262144

# creates a file HDF5 and saturate it with attributes
def OpenHDF5(fname,params,pres, temp, vms, wns, Np, Nt, Nvms, Nwn):
    global FLAG_DEBUG_PRINT
    INDEX_TEMP = '%02d'%(int(params[10][1]))
    INDEX_Qbrd_TEMP = '%02d'%(int(params[15][1]))
#    print(INDEX_TEMP)
    try:
        f = h5py.File(fname, mode='w')
        global FLAG_OPENED_HDF5
        FLAG_OPENED_HDF5 = True

        print('*********\nHDF5 file %s is opened well\n*********'%(fname))

# saturating the attributes
        [f.attrs.__setitem__(item[0],item[1]) for item in params[:6] ]

        Index_abs = INDEX_TEMP
        Index_broad = INDEX_Qbrd_TEMP
        dataset_name = 'Gas_'+Index_abs+'_Absorption'
        # ABSCO User Guide v5.0 names the broadener axis "VMR"; the historical
        # X-MASS name "VMS" is available via the Legacy_VMS_names switch
        LEGACY_VMS = readSwitchByName(params, 'Legacy_VMS_names', False)
        dataset_broadname = 'Broadener_'+Index_broad+('_VMS' if LEGACY_VMS else '_VMR')
        ds_coef = f.create_dataset(dataset_name,shape=(Np,Nt,Nvms,Nwn),dtype='float64',
                                   chunks=(1,1,1,min(Nwn,HDF5_WN_CHUNK)),
                                   compression="gzip", compression_opts=4,
                                   shuffle=True, fillvalue=np.nan)
        ds_coef.attrs.__setitem__('addl_ident', params[7][1])
        ds_coef.attrs.__setitem__('gas_name', params[8][1])
        ds_coef.attrs.__setitem__('comment', params[9][1])
        ds_coef.attrs.__setitem__('units', 'cm2/molecule')
        # failed grid points keep the fill value
        ds_coef.attrs.__setitem__('missing_value', np.nan)

        ds_index = f.create_dataset('Gas_Index',data=INDEX_TEMP)
        ds_pres = f.create_dataset('Pressure',data=pres)
        ds_pres.attrs.__setitem__('units', 'atm')
        ds_temp = f.create_dataset('Temperature',data=temp)
        ds_temp.attrs.__setitem__('units', 'K')
        ds_vms = f.create_dataset(dataset_broadname,data=vms)
        ds_vms.attrs.__setitem__('broadener_name',
                                 readParamByName(params, 'broadener_name', 'air'))
        ds_Qbrd = f.create_dataset('Broadener_Index',data=INDEX_Qbrd_TEMP)

        ds_wns = f.create_dataset('Wavenumber',data=wns)
        ds_wns.attrs.__setitem__('units', 'cm-1')
        
#        print(ds_pres[()])
        if (FLAG_DEBUG_PRINT):
            print('*** DEBUG: Attributes ***')
            print('Attributes -- root')
            [print(item, f.attrs[item]) for item in f.attrs.keys()]
            print(f.keys())
            [print(f[item]) for item in f.keys()]
            for item in f.keys():
                print('    ',f[item].attrs.keys())
                for jtem in f[item].attrs.keys():
                    print('        ',item,f[item].attrs[jtem])
                # if type(item) not ''
            #     [print(jtem, item[jtem]) for jtem in item.keys()]
            print('*** END: Attributes ***\n')
        return f
    except FileExistsError:
        print('Attempt to re-write file!')
        sys.exit()
    else:
        err = Exception
        print("Unexpected %s"%(err))
        sys.exit()

# closes the HDF5 file 
def CloseHDF5(ftype):
    try:
        global FLAG_OPENED_HDF5 
        if ((FLAG_OPENED_HDF5 != True) or (ftype.__repr__()=='<Closed HDF5 file>')):
            raise FileNotFoundError
        else:
            ftype.close()
            FLAG_OPENED_HDF5 = False
            return
    except FileNotFoundError:
        print('File to close is not found or already closed')
        sys.exit()
    else:
        err = Exception
        print("Unexpected %s"%(err))
        sys.exit()

def SaveHDF5(ftype, p,t,vms,coef_):
    try:
        global FLAG_OPENED_HDF5 
        if ((FLAG_OPENED_HDF5 != True) or (ftype.__repr__()=='<Closed HDF5 file>')):
            raise FileNotFoundError
        else:
            

            return
    except FileNotFoundError:
        print('File to work with is not found or already closed')
        sys.exit()
    else:
        err = Exception
        print("Unexpected %s"%(err))
        sys.exit()
    
def UpdateHDF5(ftype, pTVMS, ipTVMS, param):

    IndexMol = int(param[10][1])
    IndexBroad = int(param[15][1])

    Index_abs = '%02d'%(int(ftype['Gas_Index'][()]))
    dataset_name = 'Gas_'+Index_abs+'_Absorption'

    set_abs = ftype[dataset_name][()]
    for i, tptvms in enumerate(pTVMS):
        tp, tt, tv = tptvms[0], tptvms[1], tptvms[2]
        print('Opening %d file out of %d'%(i,len(pTVMS)))

        CoefFileName = './datafiles/%06.2fT_Id%02d_%06.4eatm_IdBroad%02d_%06.4fVMS_H2O_SDV_hitran2020.dat'%(tt,IndexMol,tp,IndexBroad,tv)
        coeff = ((np.loadtxt(CoefFileName)).T)[1]


        set_abs[ipTVMS[i][0][0]][ipTVMS[i][1][0]][ipTVMS[i][2][0]][:] = coeff
        coeff = []    
    ftype[dataset_name][()] = set_abs

    

    return ftype
