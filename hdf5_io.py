import h5py
import sys
import numpy as np
from initial import FLAG_DEBUG_PRINT

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
        dataset_broadname = 'Broadener_'+Index_broad+'_VMS'
        ds_coef = f.create_dataset(dataset_name,shape=(Np,Nt,Nvms,Nwn),dtype='float64', compression="gzip", compression_opts=4)
        ds_coef.attrs.__setitem__('addl_ident', '')
        ds_coef.attrs.__setitem__('gas_name', 'CO')
        ds_coef.attrs.__setitem__('comment', '')

        ds_index = f.create_dataset('Gas_Index',data=INDEX_TEMP)
        ds_pres = f.create_dataset('Pressure',data=pres)
        ds_temp = f.create_dataset('Temperature',data=temp)
        ds_vms = f.create_dataset(dataset_broadname,data=vms)
        ds_vms.attrs.__setitem__('broadener_name','h2o')
        ds_Qbrd = f.create_dataset('Broadener_Index',data=INDEX_Qbrd_TEMP)
        
        ds_wns = f.create_dataset('Wavenumber',data=wns)
        
#        print(ds_pres[()])
        if (FLAG_DEBUG_PRINT):
            print('*** DEBUG: Attributes ***')
            [print(item, f.attrs[item]) for item in f.attrs.keys()]
            print(f.keys())
            [print(f[item]) for item in f.keys()]
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

    for i, tptvms in enumerate(pTVMS):
        tp, tt, tv = tptvms[0], tptvms[1], tptvms[2]
        print('Opening %d file out of %d'%(i,len(pTVMS)))

        CoefFileName = './datafiles/%06.2fT_Id%02d_%06.4eatm_IdBroad%02d_%06.4fVMS_H2O_SDV_hitran2020.dat'%(tt,IndexMol,tp,IndexBroad,tv)
        coeff = ((np.loadtxt(CoefFileName)).T)[1]


        set_abs = ftype[dataset_name][()]
        set_abs[ipTVMS[i][0][0]][ipTVMS[i][1][0]][ipTVMS[i][2][0]][:] = coeff
        ftype[dataset_name][()] = set_abs

        coeff = []    
    

    return ftype
