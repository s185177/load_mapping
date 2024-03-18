# -*- coding: utf-8 -*-
"""
Created on Tue Nov 14 14:11:03 2023

@author: runep
"""
import os
import re
import glob
import platform
from os.path import join
from copy import deepcopy
import numpy as np
import matplotlib.pyplot as plt
import nibabel as nib
import scipy.ndimage as ndi
import time

if platform.system() == 'Windows':
    # os.chdir('D:\Total_segmentator')
    # os.chdir(r'C:\Users\runep\OneDrive\Skrivebord\Work\total_segmentator\example_data')
    os.chdir(r'C:\Users\runep\OneDrive - Danmarks Tekniske Universitet\DTU\Thesis\Code\load_mapping')


import matplotlib.pylab as pylab
params = {'legend.fontsize': 'x-large',
           'figure.figsize': (16, 8),
         'axes.labelsize': 'xx-large',
         'axes.titlesize':'x-large',
         'xtick.labelsize':'x-large',
         'ytick.labelsize':'x-large'}
pylab.rcParams.update(params)
del params

def get_bounding_box(mask_slice, pad=0):
    locs = np.nonzero(mask_slice)
    ex_pts = [(np.min(locs[1]),np.max(locs[1])), 
              (np.min(locs[0]),np.max(locs[0]))]
    
    bounding_box = np.s_[ex_pts[1][0]-pad:ex_pts[1][1]+pad, 
                         ex_pts[0][0]-pad:ex_pts[0][1]+pad]
    
    return bounding_box

def get_xy_distribution(im_slice, tissue_slice, transpose_key, plot):
    if transpose_key:
        intensity_x = np.sum(im_slice, axis=1)
        intensity_y = np.sum(im_slice, axis=0)
        im_slice = im_slice.T
        tissue_slice = tissue_slice.T
        
    else:
        intensity_x = np.sum(im_slice, axis=0)
        intensity_y = np.sum(im_slice, axis=1)
    
    if plot:
        fig = plt.figure(layout='constrained')
        ax = fig.add_gridspec(top=0.75, right=0.75).subplots()
        ax.set(aspect=1)
        ax_x = ax.inset_axes([0, 1.05, 1, 0.25], sharex=ax)
        ax_y = ax.inset_axes([1.05, 0, 0.25, 1], sharey=ax)
        
        ax.imshow(tissue_slice)
        ax.imshow(im_slice, alpha = 0.5*im_slice, cmap='jet')
        ax_x.bar(np.arange(len(intensity_x)), intensity_x)
        ax_x.tick_params(axis='x', labelbottom=False)
        ax_x.grid()
        ax_y.barh(np.arange(len(intensity_y)), intensity_y)
        ax_y.tick_params(axis='y', labelleft=False)
        ax_y.grid()
    
        plt.show()
    
    # return intensity_x, intensity_y

def init_image(self, path):
    self.file = nib.load(path)
    self.header = self.file.header
    self.pixel_spacing = self.header['pixdim'][1]
    self.slice_thickness = self.header['pixdim'][3]    

def get_data(file):
    return np.moveaxis(np.flip(file.get_fdata(),[0,1,2]),0,1)

class PET:
    def __init__(self, path):
        init_image(self, path)
        
class CT:
    def __init__(self, path):
        init_image(self, path)

class segmentation:
    def __init__(self, path):
        pattern = re.compile('[A-Za-z](\d+)')  
        files = [f for f in os.listdir(path) if 'vertebrae' in f]
        # sorted(files,key = lambda x: (x[0], int(x[1:])))
        tissue = [pattern.search(f).group() for f in files]
        self.files = dict(zip(tissue, [nib.load(join(path,f)) for f in files]))
        self.tissue = sorted(tissue,key = lambda x: (x[0], int(x[1:])))

#%%
sub = 'data/S01'
obj_ct = CT(join(sub,'AC_CT_3_2_trunc.nii.gz'))
obj_pet = PET(join(sub,'PET_WB_NaF_resamp_trunc.nii.gz'))
obj_seg = segmentation(join(sub,'segmentations'))

obj_ct.A = get_data(obj_ct.file)
obj_pet.A = get_data(obj_pet.file)
obj_seg.A = get_data(obj_seg.files['L1'])

#%%
slices = np.unique(np.nonzero(obj_seg.A)[-1])

mask = deepcopy(obj_pet.A)
mask[obj_seg.A == 0] = 0
#%%
# %matplotlib qt
z_sum_test = np.mean(obj_seg.A, axis=2)
z_bound = get_bounding_box(z_sum_test)

fig, (ax1, ax2) = plt.subplots(1,2)

for ind in slices:
    z_0 = np.sum(obj_seg.A[...,ind], axis=0)
    z_1 = np.sum(obj_seg.A[...,ind], axis=1)
    
    ax1.plot(z_0)
    ax2.plot(z_1)
z_sum_test_x = np.sum(z_sum_test,axis=0)
z_sum_test_y = np.sum(z_sum_test,axis=1)
z_sum_test_x[z_sum_test_x == 0] = np.nan
z_sum_test_y[z_sum_test_y == 0] = np.nan
ax1.vlines(np.nanargmax(z_sum_test_x), 0, 35 ,'k')

y_locs = np.where(~np.isnan(z_sum_test_y))[0]
mid_ind = y_locs[len(y_locs)//2]

while (z_sum_test_y[mid_ind] > z_sum_test_y[mid_ind - 1]) | (z_sum_test_y[mid_ind] > z_sum_test_y[mid_ind + 1]):
    side = np.argmin((z_sum_test_y[mid_ind - 1], z_sum_test_y[mid_ind + 1]))
    
    if side == 0:
        mid_ind -= 1
    elif side == 1:
        mid_ind += 1
        

ax2.vlines(mid_ind, 0, 35 ,'k')
ax1.plot(z_sum_test_x, '-k')
ax2.plot(z_sum_test_y, '-k')
ax1.set_xlim([200, 300])
ax2.set_xlim([275, 350])


plt.show()
#%%
cor_ind = 300
sag_ind = 256
ax_ind = 200

axis = 1
bound = 1
if bound:
    cor_bound = get_bounding_box(obj_seg.A[cor_ind,...], 0)
    sag_bound = get_bounding_box(obj_seg.A[:,sag_ind,:], 0)
    ax_bound = get_bounding_box(obj_seg.A[...,ax_ind], 0)
    
    cor_bound = (cor_ind, cor_bound[0], cor_bound[1])
    sag_bound = (sag_bound[0], sag_ind, sag_bound[1])
    ax_bound = (ax_bound[0], ax_bound[1], ax_ind)
else:
    cor_bound = (cor_ind, slice(0,-1), slice(0,-1))
    sag_bound = (slice(0,-1), sag_ind, slice(0,-1))
    ax_bound = (slice(0,-1), slice(0,-1), ax_ind)

#%%
final_bound = (sag_bound[0], ax_bound[1], cor_bound[2])
z_slice_final = obj_seg.A[final_bound]
final_dim = z_slice_final.shape

z_slice_seg = obj_seg.A[final_bound]
z_slice_ct = obj_ct.A[final_bound]

#%%
np.unique(np.nonzero(z_slice_final))
z_samp = z_slice_final[...,5]
get_xy_distribution(z_samp, z_slice_ct[...,5], 0, 1)
#%%


def get_local_minima(sum_img):
    sum_img = deepcopy(sum_img)
    sum_img[sum_img == 0] = np.nan

    y_locs = np.where(~np.isnan(sum_img))[0]
    mid_ind = y_locs[len(y_locs)//2]
    
    k_nearest = len(y_locs)//4
    lowest = np.argmin(sum_img[mid_ind - k_nearest: mid_ind + k_nearest])
    if np.isnan(lowest):
        k_nearest = len(y_locs)//5
        lowest = np.argmin(sum_img[mid_ind - k_nearest: mid_ind + k_nearest])

    return mid_ind, mid_ind - k_nearest + lowest

for ind in range(5,25):
    z_samp = z_slice_final[...,ind]
    z_samp_cor_sum = np.sum(z_samp, axis=0)
    z_samp_sag_sum = np.sum(z_samp, axis=1)
    
    z_x_mid, z_x_test = get_local_minima(z_samp_cor_sum)
    z_y_mid, z_y_test = get_local_minima(z_samp_sag_sum)
    
    plt.figure()
    plt.title(ind)
    plt.plot(z_samp_cor_sum, 'b')
    plt.plot(z_samp_sag_sum, 'r')
    plt.vlines(z_x_test, 0, 35,'k')
    plt.vlines(z_y_test, 0, 35, 'k')
    
    plt.vlines(z_x_mid, 0, 35,'c')
    plt.vlines(z_y_mid, 0, 35, 'c')
    plt.show()
#%%

#%% Dette er en ny test
init_guess = ndi.center_of_mass(z_slice_final)

cor_pos = 15
sag_pos = 20
ax_pos  = 12

fig, axs = plt.subplots(1,3)
axs[0].imshow(z_slice_final[cor_pos,...].T)
axs[0].plot(init_guess[1], init_guess[2], '+r', markersize=20)
axs[0].plot(sag_pos, ax_pos, '+b', markersize=15)
axs[0].set_title('Coronal')
axs[0].set_xlabel('Sagittal')
axs[0].set_ylabel('Axial')


axs[1].imshow(z_slice_final[:,sag_pos,:].T)
axs[1].plot(init_guess[0], init_guess[2], '+r', markersize=20)
axs[1].plot(cor_pos, ax_pos, '+b', markersize=15)
axs[1].set_title('Sagittal')
axs[1].set_xlabel('Coronal')
axs[1].set_ylabel('Axial')


axs[2].imshow(z_slice_final[...,ax_pos])
axs[2].plot(init_guess[1], init_guess[0], '+r', markersize=20)
axs[2].plot(sag_pos, cor_pos, '+b', markersize=15)
axs[2].set_title('Axial')
axs[2].set_xlabel('Sagittal')
axs[2].set_ylabel('Coronal')


axs[0].set_aspect(obj_ct.slice_thickness/obj_ct.pixel_spacing)
axs[1].set_aspect(obj_ct.slice_thickness/obj_ct.pixel_spacing)


# for ax in axs.ravel():
    # ax.axis('off')
plt.show()
#%% Show Pixel distribution along all axes
z_cor_sum = np.sum(z_slice_final, axis=0) # [43,25] [sag, ax]
z_sag_sum = np.sum(z_slice_final, axis=1) # [46,25] [cor, ax]
z_ax_sum  = np.sum(z_slice_final, axis=2) # [46,43] [cor, sag]

vline_ax_pos = 20

cor_pos = 20
fig, axs = plt.subplots(3,3)
axs[0,0].imshow(z_slice_final[cor_pos,...].T)
axs[0,0].vlines(vline_ax_pos,3,20, 'k')
axs[1,0].imshow(z_slice_final[:, sag_pos,:].T)
axs[2,0].imshow(z_slice_final[...,ax_pos])

for ind in range(final_dim[0]):
    axs[0,1].plot(z_ax_sum[ind,:].T, '--', alpha=.5)
    axs[0,2].plot(z_sag_sum[ind,:].T, '--', alpha=.5)

axs[0,1].plot(np.mean(z_ax_sum, axis=0), 'k')
axs[0,1].vlines(vline_ax_pos,0,15,'k')
axs[0,1].set_xlabel('Sagittal')
axs[0,2].plot(np.mean(z_sag_sum, axis=0).T, 'k')
axs[0,2].set_xlabel('Axial')

for ind in range(final_dim[1]):
    axs[1,1].plot(z_ax_sum[:,ind], '--', alpha=.5)
    axs[1,2].plot(z_cor_sum[ind,:].T, '--', alpha=.5)

axs[1,1].plot(np.mean(z_ax_sum, axis=1), 'k')
axs[1,1].set_xlabel('Coronal')
axs[1,2].plot(np.mean(z_cor_sum, axis=0).T, 'k')
axs[1,2].set_xlabel('Axial')

for ind in range(final_dim[2]):
    axs[2,1].plot(z_cor_sum[:,ind], '--', alpha=.5)
    axs[2,2].plot(z_sag_sum[:,ind], '--', alpha=.5)

axs[2,1].plot(np.mean(z_cor_sum, axis=1), 'k')
axs[2,1].set_xlabel('Sagittal')
axs[2,2].plot(np.mean(z_sag_sum, axis=1), 'k')
axs[2,2].set_xlabel('Coronal')

axs[0,0].set_xlabel('Sagittal')
axs[0,0].set_ylabel('Axial')
axs[1,0].set_xlabel('Coronal')
axs[1,0].set_ylabel('Axial')
axs[2,0].set_xlabel('Sagittal')
axs[2,0].set_ylabel('Coronal')

axs[0,0].set_aspect(obj_ct.slice_thickness/obj_ct.pixel_spacing)
axs[1,0].set_aspect(obj_ct.slice_thickness/obj_ct.pixel_spacing)

plt.tight_layout()
plt.show()
#%% Test with sum curves along all axes
cor_ax = np.arange(0, final_dim[0])
sag_ax = np.arange(0, final_dim[1])
grid = np.meshgrid(cor_ax, sag_ax)

mean_x = np.ones((final_dim[1]))*0
mean_y = np.ones((final_dim[0]))*0
mean_sag = np.mean(z_ax_sum, axis=0)
mean_cor = np.mean(z_ax_sum, axis=1)

fig = plt.figure()
ax = plt.axes(projection='3d')
ax.plot(mean_x, np.arange(0, final_dim[1]), mean_sag, '-k')
ax.plot(np.arange(0, final_dim[0]), mean_y, mean_cor, '-k')

ax.contour3D(grid[0], grid[1], z_ax_sum.T, 500, alpha=.1)
# ax.set_ylim([final_dim[1]+5, 0])
# ax.set_xlim([final_dim[0], 0])
ax.set_xlabel('Coronal')
ax.set_ylabel('Sagittal')
ax.set_zlabel('Axial')
# plt.legend()
plt.show()

del cor_ax, sag_ax, grid, mean_x, mean_y, mean_sag, mean_cor, fig, ax
#%% Test with flood fill segmentation
import skimage.segmentation as seg

for ind in range(43):
    flood_mask_sag = seg.flood(z_slice_final[:,ind,:], 
                               (int(init_guess[0]), int(init_guess[2])))
    
    fig, axs = plt.subplots(1,2)
    axs[0].imshow(z_slice_final[:,ind,:].T)
    axs[1].imshow(flood_mask_sag.T)
    
    axs[0].set_aspect(obj_ct.slice_thickness/obj_ct.pixel_spacing)
    axs[1].set_aspect(obj_ct.slice_thickness/obj_ct.pixel_spacing)
    

    axs[0].axis('off')
    axs[1].axis('off')
    plt.show()
#%% Plot all three planes together
# Coronal
fig, (a1, a2, a3) = plt.subplots(1, 3)
a1.set_title('Coronal')
a1.imshow(obj_ct.A[cor_bound].T, cmap='gray')
# a1.imshow(obj_pet.A[cor_bound].T)
# a1.imshow(mask[cor_bound].T)
a1.imshow(obj_seg.A[cor_bound].T, alpha=.4*obj_seg.A[cor_bound].T)
a1.set_aspect(obj_ct.slice_thickness/obj_ct.pixel_spacing)
if axis:
    a1.axis('off')
elif bound:
    a1.set_yticklabels(np.arange(cor_bound[2].start, cor_bound[2].stop))
    a1.set_xticklabels(np.arange(cor_bound[1].start, cor_bound[1].stop))

# Sagital
a2.set_title('Sagittal')
a2.imshow(obj_ct.A[sag_bound].T, cmap='gray')
# a2.imshow(obj_pet.A[sag_bound].T)
# a2.imshow(mask[sag_bound].T)
a2.imshow(obj_seg.A[sag_bound].T, alpha=.4*obj_seg.A[sag_bound].T)
a2.set_aspect(obj_ct.slice_thickness/obj_ct.pixel_spacing)
if axis:
    a2.axis('off')
elif bound:
    a2.set_yticklabels(np.arange(sag_bound[2].start, sag_bound[2].stop))
    a2.set_xticklabels(np.arange(sag_bound[0].start, sag_bound[0].stop))

# Axial
a3.set_title('Axial')
a3.imshow(obj_ct.A[ax_bound], cmap='gray')
# a3.imshow(obj_pet.A[ax_bound])
# a3.imshow(mask[ax_bound])
a3.imshow(obj_seg.A[ax_bound], alpha=.4*obj_seg.A[ax_bound])
if axis:
    a3.axis('off')
elif bound:
    a3.set_yticklabels(np.arange(ax_bound[0].start, ax_bound[0].stop))
    a3.set_xticklabels(np.arange(ax_bound[1].start, ax_bound[1].stop))

plt.tight_layout()
plt.show()

# del a1, a2, a3, bound, cor_ind, sag_ind, ax_ind, fig

#%% Plot distribution
switch = 0
# zz = np.nonzero(obj_seg.A)
# np.unique(zz[switch])
slices = np.unique(np.nonzero(obj_seg.A)[switch])
# slices = np.arange(slices[0], slices[-1]+1)

for ind in slices:
    if switch == 0:
        # Coronal
        bound = get_bounding_box(obj_seg.A[ind,...], 5)
        bound = (ind, bound[0], bound[1])
    elif switch == 1:
        # Sagittal
        bound = get_bounding_box(obj_seg.A[:,ind,:], 20)
        bound = (bound[0], ind, bound[1])
    elif switch == 2:
        # Axial
        print(ind)
        bound = get_bounding_box(obj_seg.A[...,ind], 5)
        bound = (bound[0], bound[1], ind)
    
    get_xy_distribution(obj_seg.A[bound], obj_ct.A[bound], 1, 1)
