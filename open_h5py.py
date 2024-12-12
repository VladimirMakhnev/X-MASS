print("\n*****\nSupplementary tool to X-MASS.\nABSCO-table-formatted file opener.\n*****\n")
print("V.Yu. Makhnev, 12 of Dec 2024\n")

# TESTING TOOL TO CHECK IF HDF5 FILE IS NOT EMPTY
import numpy as np
import sys
import matplotlib as mpl
import matplotlib.pyplot as plt




molec_id = 2

HDF5filename = 'test.LP.HITRAN2020.25wing'

# -*- coding: utf-8 -*-
import h5py
fnam = h5py.File(HDF5filename+'.hdf5',mode='r')
for i, item in enumerate(fnam.keys()):
    print(fnam[item])

# pick one cross-section 

# for ip in np.arange(10):
#     for it in np.arange(13):
#         for ivms in np.arange(3):
#             wn = fnam['Wavenumber'][()]
#             xsec = fnam['Gas_%02d_Absorption'%molec_id][()][ip,it,ivms,:]
#             print(len(xsec[xsec>0]))

# print chosen x-sec
ip, it, ivms = 1,2,1


wn = fnam['Wavenumber'][()]
xsec1 = fnam['Gas_%02d_Absorption'%molec_id][()][ip,it,ivms,:]
xsec2 = fnam['Gas_%02d_Absorption'%molec_id][()][ip+1,it,ivms,:]
xsec3 = fnam['Gas_%02d_Absorption'%molec_id][()][ip,it+1,ivms,:]
# print(wn)
# print(xsec)
p1 = fnam['Pressure'][()][ip]
T1 = fnam['Temperature'][()][ip,it]
vms = fnam['Broadener_00_VMS'][()][ivms]
p2 = fnam['Pressure'][()][ip+1]
T2 = fnam['Temperature'][()][ip+1,it]
p3 = fnam['Pressure'][()][ip]
T3 = fnam['Temperature'][()][ip,it+1]






#########################################################################################

nu_start =      fnam['Wavenumber'][()][0]
nu_end   =      fnam['Wavenumber'][()][-1]

y_start = min(xsec1)
y_end   =  max(xsec1)


figure1 = plt.figure(figsize=(16,8),dpi=500)
title_band='xsec of h2o from hdf5 file'
ax1 = figure1.add_subplot(111) #(211)
#            ax2 = figure1.add_subplot(212, sharex=ax1)

ax1.set_xlim(nu_start,nu_end)
ax1.set_ylim(y_start,y_end)
#secax = ax1.secondary_xaxis('top', functions=(cm2mum,cm2mum))
#secax.set_xlabel(r'cm$^{-1}$')
ax1.set_xlabel(r'$cm^{-1}$')
# ax1.set_xlabel(r'$cm^{-1}/(molecule*cm^{-2})$')
# ax1.set_ylabel(r'\%')

ax1.set_title(title_band, y=1.05)
# ax1.set_xscale('log')
ax1.set_yscale('log')

ax1.plot(wn,xsec1,label='Cross-section of H$_2$O, p=%4.2e bar, T=%4.2f K, VMS=%4.2f'%(p1,T1,vms),alpha=0.5)
ax1.plot(wn,xsec2,label='Cross-section of H$_2$O, p=%4.2e bar, T=%4.2f K, VMS=%4.2f'%(p2,T2,vms),alpha=0.5)
ax1.plot(wn,xsec3,label='Cross-section of H$_2$O, p=%4.2e bar, T=%4.2f K, VMS=%4.2f'%(p3,T3,vms),alpha=0.5)
# ax1.scatter(oldllm['S'],(oldllm['S']-newllm['S'])/newllm['S']*100, label='Old/New-1 *100\%',s=4.0,marker='x',color=(0.9,0.1,0.1))
# ax1.scatter(, s=5.0,marker='+', color=(0.1,0.9,0.1), label='pRT grid')
# ax1.plot(lambda_pRT,co_hapi_int,label='Interpolated',alpha=0.5)
#            ax1.scatter(lambda_pRT,co_sent/N_a*mass_amu,label='Sent',marker='+',s=150.)

# ax1.hlines(y=0.0, xmin=nu_start, xmax=nu_end)

# ax1.hlines(y=2.0, xmin=nu_start, xmax=nu_end)

# ax1.hlines(y=-2.0, xmin=nu_start, xmax=nu_end)

ax1.legend()

#            ax2.scatter(lambda_hapi, ratio_co)

name_img = "./"+HDF5filename+"_hdf5.png"
plt.savefig(name_img,bbox_inches='tight', transparent=False)
plt.close()

####################################################################################################


fnam.close()