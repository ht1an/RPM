#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Nov  5 09:26:12 2020
To test the reduced proper motion to select the K giant stars from dwarf stars
@author: htian
"""
#%%
from astropy.io import fits
import numpy as np
import matplotlib.pyplot as plt
#%% set the parameters
dpath = "/Users/htian/Documents/work/data/rothalo/"
dfile = "dr5GaiaDR2_v1_all_photos_EWs_dist_NODUPLICATE.fits"
bin_bprp = np.linspace(-1,7,81)
bin_rpm = np.linspace(-5,20,51)
#%% read the data
dt = fits.open(dpath + dfile)
data = dt[1].data
ind_snr = data["snrg"]>20
#%%
glon = data["glon"][ind_snr]
glat = data["glat"][ind_snr]
pmra = data["pmra"][ind_snr]
pmdec = data["pmdec"][ind_snr]
teff = data["teff"][ind_snr]
logg = data["logg"][ind_snr]
magG = data["phot_g_mean_mag"][ind_snr]
bprp = data["bp_rp"][ind_snr]
prlx = data["parallax"][ind_snr]
prlx_err = data["parallax_error"][ind_snr]
snrg = data["snrg"][ind_snr]
KGtag = data["TrueKGiants_Liu14"][ind_snr]
#%%
ind_KGiant = KGtag==True
pm_tot = np.sqrt(pmra**2 + pmdec**2)
rpm = magG + 5*np.log10(pm_tot)-10
H_RPMM, Hx, Hy = np.histogram2d(bprp, rpm, bins=(bin_bprp, bin_rpm))
H_KG, _, _ = np.histogram2d(bprp[ind_KGiant], rpm[ind_KGiant], bins=(bin_bprp, bin_rpm))
#%%

fig = plt.figure(figsize=(10,5))
ax = fig.add_axes([0.12,0.18,0.85,0.79])
plt.pcolormesh(Hx[1:]-(Hx[1]-Hx[0])*0.5, Hy[1:]-(Hy[1]-Hy[0])*0.5, np.log10(H_RPMM.T+1),cmap='jet')
plt.colorbar()
plt.contour(Hx[1:]-0.5*(Hx[1]-Hx[0]), Hy[1:]-0.5*(Hy[1]-Hy[0]),np.log10(H_KG.T+1),levels = np.linspace(0,4,10),linestyles='--')
plt.ylim(20, -5)
plt.xticks(fontsize=15)
plt.yticks(fontsize=15)
plt.ylim(Hy[-1]-(Hy[1]-Hy[0])*0.5, Hy[1]-(Hy[1]-Hy[0])*0.5)
plt.xlabel("Bp-Rp",fontsize=15)
plt.ylabel("$G+5*log_{10}(PM)-10$",fontsize=15)
plt.savefig("c_rpm.png")
#%%  this is the test on the tag from Kgiant of Liu 14
c_bprp = np.linspace(0,4,41)
c_rpm = np.linspace(0, 15, 51)
mm = np.zeros((51,41))
nn = np.zeros_like(mm)
for i in range(51):
    
    for j in range(41):
        ind_rpm = (bprp>c_bprp[j]) & (rpm<c_rpm[i])
        nn[i,j] = len(KGtag[ind_rpm][KGtag[ind_rpm]==True])
        if len(KGtag[ind_rpm])>0:
               mm[i,j] = len(KGtag[ind_rpm][KGtag[ind_rpm]==True]) / len(KGtag[ind_rpm])
#%%
cc_mat, rr_mat = np.meshgrid(c_bprp, c_rpm)
ind_max = mm[::-1,:]==np.max(mm)
cc_edge, rr_edge = cc_mat[ind_max], rr_mat[::-1][ind_max]

fig = plt.figure(figsize=(6,4))
ax = fig.add_axes([0.17,0.2,0.78,0.75])
plt.pcolormesh(c_bprp, c_rpm[::-1], mm[::-1,:], cmap="jet")
plt.colorbar()
plt.axvline(x=cc_edge,linestyle='dashed',color="black")
plt.axhline(y=rr_edge,linestyle='dashed',color="black")
plt.xticks(fontsize=15)
plt.yticks(fontsize=15)
plt.xlabel("color cut",fontsize=15)
plt.ylabel("RPM cut",fontsize=15)
plt.ylim(np.max(c_rpm), np.min(c_rpm))
plt.text(2.5,10,f"f={np.max(mm*100):.1f}%",fontsize=15,color='#ff7f0e')
plt.text(2.5,11.5,f"Bp-Rp={cc_edge[0]}",fontsize=15,color='#ff7f0e')
plt.text(2.5,13,f"RPM={rr_edge[0]}",fontsize=15,color='#ff7f0e')
plt.savefig("cut_fig.png")

fig = plt.figure(figsize=(6,4))
ax = fig.add_axes([0.17,0.2,0.78,0.75])
plt.pcolormesh(c_bprp, c_rpm[::-1], np.log10(nn[::-1,:]+1), cmap="jet")
plt.colorbar()
plt.contour(c_bprp, c_rpm[::-1], mm[::-1,:],levels=np.linspace(0,1,11))
plt.axvline(x=cc_edge,linestyle='dashed',color="black")
plt.axhline(y=rr_edge,linestyle='dashed',color="black")
plt.xticks(fontsize=15)
plt.yticks(fontsize=15)
plt.xlabel("color cut",fontsize=15)
plt.ylabel("RPM cut",fontsize=15)
plt.ylim(np.max(c_rpm), np.min(c_rpm))
plt.text(2.5,10,f"f={np.max(mm*100):.1f}%",fontsize=15,color='#ff7f0e')
plt.text(2.5,11.5,f"Bp-Rp={cc_edge[0]}",fontsize=15,color='#ff7f0e')
plt.text(2.5,13,f"RPM={rr_edge[0]}",fontsize=15,color='#ff7f0e')
plt.savefig("cut_nn_fig.png")

#%%
# aa = np.array([[1,2,3,4],[2,4,6,8],[3,6,9,12]])
# print(aa)
# plt.pcolormesh(aa[::-1,:],cmap='jet')
# plt.colorbar()

#%%  this is the test on the tag from Kgiant of Liu 14
# c_bprp = np.linspace(0,4,41)
# c_rpm = np.linspace(15, 0, 51)
mm_tg = np.zeros((51,41))
for i in range(51):
    
    for j in range(41):
        ind_rpm = (bprp>c_bprp[j]) & (rpm<c_rpm[i])
        ind_rpm_kg = (teff[ind_rpm]<5600) & (logg[ind_rpm]<4)
        if len(KGtag[ind_rpm])>0:
               mm_tg[i,j] = len(KGtag[ind_rpm][ind_rpm_kg]) / len(KGtag[ind_rpm])
#%%
# cc_mat, rr_mat = np.meshgrid(c_bprp, c_rpm)
ind_max = mm_tg[::-1]==np.max(mm_tg)
cc_edge, rr_edge = cc_mat[ind_max], rr_mat[::-1][ind_max]

fig = plt.figure(figsize=(6,4))
ax = fig.add_axes([0.17,0.2,0.78,0.75])
plt.pcolormesh(c_bprp, c_rpm[::-1], mm_tg[::-1,:], cmap="jet")
plt.colorbar()
plt.vlines(x=cc_edge,ymin=np.min(c_rpm),ymax = np.max(c_rpm),linestyle='dashed',color="black")
plt.hlines(y=rr_edge,xmin=np.min(c_bprp),xmax = np.max(c_bprp),linestyle='dashed',color="black")
plt.xticks(fontsize=15)
plt.yticks(fontsize=15)
plt.xlabel("color cut",fontsize=15)
plt.ylabel("RPM cut",fontsize=15)

plt.text(2.5,10,f"f={np.max(mm_tg*100):.1f}%",fontsize=15,color='#ff7f0e')
plt.text(2.5,11.5,f"Bp-Rp={cc_edge[0]}",fontsize=15,color='#ff7f0e')
plt.text(2.5,13,f"RPM={rr_edge[0]}",fontsize=15,color='#ff7f0e')
plt.ylim(np.max(c_rpm), np.min(c_rpm))
plt.savefig("cut_fig_onHRD.png")

