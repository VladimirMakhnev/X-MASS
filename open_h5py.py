print("\n*****\nSupplementary tool to X-MASS.\nABSCO-table-formatted file opener.\n*****\n")
print("V.Yu. Makhnev, 16 of Sept 2024\n")

# TESTING TOOL TO CHECK IF HDF5 FILE IS NOT EMPTY
import numpy as np
import sys
import matplotlib as mpl
import matplotlib.pyplot as plt




molec_id = 6



# -*- coding: utf-8 -*-
import h5py
fnam = h5py.File('01.H2O.SDV.HITRAN2020.25wing.hdf5',mode='r')
for i, item in enumerate(fnam.keys()):
    print(fnam[item])

# pick one cross-section 

ip, it, ivms = 1,0,1
wn = fnam['Wavenumber'][()]
xsec = fnam['Gas_%02d_Absorption'%molec_id][()][ip,it,ivms,:]
    

# print chosen x-sec
print(wn)
print(xsec)

p = fnam['Pressure'][()][ip]
T = fnam['Temperature'][()][ip,it]
vms = fnam['Broadener_00_VMS'][()][ivms]






#########################################################################################

nu_start =      fnam['Wavenumber'][()][0]
nu_end   =      fnam['Wavenumber'][()][-1]

y_start = min(xsec)
y_end   =  max(xsec)


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
# ax1.set_yscale('log')

ax1.plot(wn,xsec,label='Cross-section of H$_2$O, p=%4.2f bar, T=%4.2f K, VMS=%4.2f'%(p,T,vms))
# ax1.scatter(oldllm['S'],(oldllm['S']-newllm['S'])/newllm['S']*100, label='Old/New-1 *100\%',s=4.0,marker='x',color=(0.9,0.1,0.1))
# ax1.scatter(, s=5.0,marker='+', color=(0.1,0.9,0.1), label='pRT grid')
# ax1.plot(lambda_pRT,co_hapi_int,label='Interpolated',alpha=0.5)
#            ax1.scatter(lambda_pRT,co_sent/N_a*mass_amu,label='Sent',marker='+',s=150.)

# ax1.hlines(y=0.0, xmin=nu_start, xmax=nu_end)

# ax1.hlines(y=2.0, xmin=nu_start, xmax=nu_end)

# ax1.hlines(y=-2.0, xmin=nu_start, xmax=nu_end)

ax1.legend()

#            ax2.scatter(lambda_hapi, ratio_co)

name_img = "./test_hdf5.png"
plt.savefig(name_img,bbox_inches='tight', transparent=False)
plt.close()

####################################################################################################
