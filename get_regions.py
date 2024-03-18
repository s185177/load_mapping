# -*- coding: utf-8 -*-
"""
Created on Wed Feb  7 15:18:27 2024

@author: runep
"""
import numpy as np
import nibabel as nib

from os.path import join
from scipy.interpolate import splprep, splev
# from load_mapping import Image, Segmentation, get_bounding_box, plot_all_planes, legend
from load_mapping import plot_all_planes, load_map
import load_mapping

from time import time

#%% S02 Linux
subject = 'data/S02'
mask = load_mapping.Image(join(subject, '01', 'Segmentations.nii'))
mask.get_data()

bounding_box = load_mapping.get_bounding_box(mask.A)
mask.bound = mask.A[bounding_box]
bound_dim = mask.bound.shape
#%% S03 Windows
subject = r'E:\TotalSegmentator\Spine\S03'
scan = '02'
mask = load_mapping.Image(join(subject, scan, 'Segmentations.nii'))
mask.A = mask.get_data()

bounding_box = load_mapping.get_bounding_box(mask.A)
mask.bound = mask.A[bounding_box]
bound_dim = mask.bound.shape

mask_vert_body = load_mapping.Image(join(subject, scan, 'Segmentations_vertebrae.nii'))
mask_vert_body.A = mask_vert_body.get_data()

bounding_box_vert_body = load_mapping.get_bounding_box(mask_vert_body.A)
mask_vert_body.bound = mask_vert_body.A[bounding_box_vert_body]
bound_dim_vert_body = mask_vert_body.bound.shape
#%% Sanity check
# plane_ind = [60, 38, 150]
# plane_ind = [130, 156, 175]
# plane_ind = [300, 250, 200]
# plot_all_planes(mask.bound, 
#                 mask.slice_thickness, 
#                 mask.pixel_spacing, 
#                 plane_ind, cm='nipy_spectral', axes=1)
# plt.tight_layout()
# plt.show()
# del plane_ind
#%%
def get_slicing(array, pad=0):
    unique = np.unique(array)
    if pad == 0:
        return unique
    else:
        return unique[slice(pad,-pad)]

def split_vertebrae(img, pad_cor=0, pad_sag=0, pad_ax=[5,95]):
    slices = np.unique(np.nonzero(img)[-1])
    split_mask = np.zeros(img.shape, dtype=np.int32)
    
    center_slices = np.percentile(slices, pad_ax)
    center_ax = slice(np.where(slices == int(center_slices[0]))[0][0], 
                      np.where(slices == int(center_slices[1]))[0][0])
    # center_ax = get_slicing(slices, pad_ax)
    
    for z in slices[center_ax]:
        coords_cor, coords_sag = np.nonzero(img[...,z])
        center_cor, center_sag = get_slicing(coords_cor, pad_cor), get_slicing(coords_sag, pad_sag)
        
        for x in center_cor:
            coord_plane = np.mean(coords_sag[np.where(coords_cor == x)])
            split_mask[x, int(coord_plane), z] = 1
        for y in center_sag:
            coord_plane = np.mean(coords_cor[np.where(coords_sag == y)])
            split_mask[int(coord_plane), y, z] = 2
    return split_mask

def fit_spine(sum_array, smoothing_scalar=1e3, n_points=1000):
    y_coords, x_coords = np.nonzero(sum_array)
    x_list = np.unique(x_coords)
    y_list = np.zeros(len(x_list))
    
    for n, x in enumerate(x_list):
        max_val = 0
        ind = 0
        
        for y in list(np.where(x_coords == x)[0]):
            value = sum_array[y_coords[y], x]
            
            if value > max_val:
                max_val = value
                ind = y            
        y_list[n] = y_coords[ind]
    
    # Fit the found coordinates
    x_new = np.linspace(0,1.01, n_points)
    (tck,u), fp, _,_ = splprep([x_list, y_list], s=smoothing_scalar, full_output=1)
    new_points = splev(x_new, tck)
    return new_points

def save_data(img, name, orig_size=False):
    if orig_size:
        img = load_mapping.remove_bounding_box(img, bounding_box, mask.A.shape).astype(np.int32)
        
    nib.save(nib.Nifti1Image(np.moveaxis(np.flip(img,[0,1,2]),0,1), mask.file.affine), name + '.nii.gz')

#% Processing
mask.split_coords = np.zeros(mask.bound.shape, dtype=np.int32)
mask.spine = np.zeros(mask.bound.shape, dtype=np.int32)

t1 = time()

for key, label in load_mapping.legend.items():
    # if label == 32:
    print(key, label)
    vertebrae = np.where(mask.bound.astype(np.int32) == label, 
                             load_map._PROCESS_KEY, 
                             0)
    vertebrae = np.where((mask_vert_body.A[bounding_box].astype(np.int32) != 0) & 
                             (vertebrae == load_map._PROCESS_KEY), 
                             load_map._VERTEBRAE_KEY, 
                             vertebrae)
    mask.spine += vertebrae
    
    # Change axial slice percentile to reduce number of slices used for fitting
    mask.split_coords += split_vertebrae(np.where(vertebrae == load_map._VERTEBRAE_KEY, 1,0),
                                         0,0,[30,70])
print(time() - t1)
#%% Fit spine
t1 = time()
print("Fitting spine")
x_list_cor, y_list_cor = fit_spine(np.sum(np.where(mask.split_coords == 2, 1,0), 1), 
                                   smoothing_scalar=300)
x_list_sag, y_list_sag = fit_spine(np.sum(np.where(mask.split_coords == 1, 1,0), 0), 
                                   smoothing_scalar=10)

coords_cor = np.array([x_list_cor, y_list_cor], dtype=np.int32).T
coords_sag = np.array([x_list_sag, y_list_sag], dtype=np.int32).T

split_cor = np.zeros(bound_dim, dtype=np.int32)
split_sag = np.zeros(bound_dim, dtype=np.int32)
print("Assigning mask values")
for n in range(len(coords_cor)):
    split_sag[:coords_cor[n,1],:,coords_cor[n,0]] = load_map._ANTERIOR_KEY # 4
    split_sag[coords_cor[n,1]:,:,coords_cor[n,0]] = load_map._POSTERIOR_KEY # 8
    
    split_cor[:,:coords_sag[n,1],coords_sag[n,0]] = load_map._LEFT_KEY # 16
    split_cor[:,coords_sag[n,1]:,coords_sag[n,0]] = load_map._RIGHT_KEY # 32

mask.spine += (split_sag + split_cor)
print("Saving")
# save_data(mask.spine, 'data/split_map_S03_02')

print("Finished time taken {:.2f} seconds".format(time() - t1))
# plot_all_planes(mask.spine, mask.slice_thickness, mask.pixel_spacing, [120, 150, 180], axes=1)



# img = load_mapping.remove_bounding_box(mask.spine, bounding_box, mask.A.shape).astype(np.int32)
# nib.save(nib.Nifti1Image(np.moveaxis(np.flip(img,[0,1,2]),0,1), mask.file.affine), 'data/split_map_S03_02' + '.nii.gz')
#%%
import matplotlib.pyplot as plt
from scipy.ndimage import center_of_mass
#%%

plot_all_planes(masked_bound, mask.slice_thickness, mask.pixel_spacing, [7, 15, 10], axes=1)
plt.show()

#%%
# masked array needs to be only vertebral body
def get_coronal_split(masked_array):
    bounding_box = load_mapping.get_bounding_box(masked_array)
    bound_array = masked_array[bounding_box]
    com = center_of_mass(bound_array)
    x, y = int(com[0]), int(com[1])
    
    bound_coords = np.nonzero(bound_array[:,y,:])
    bound_coords_ant = bound_coords[1][np.where(bound_coords[0] >= x)]
    bound_coords_pos = bound_coords[1][np.where(bound_coords[0] < x)]
    mid_ant = (max(bound_coords_ant) + min(bound_coords_ant))/2
    mid_pos = (max(bound_coords_pos) + min(bound_coords_pos))/2
    
    m, c = np.linalg.lstsq(np.vstack([np.array([0, max(bound_coords[0])]), np.ones(2)]).T, 
                           np.array([mid_ant, mid_pos]), rcond=None)[0]
    # x_range = np.arange(0, np.max(bound_coords[0]) + 1)
    x_range = np.linspace(0, np.max(bound_coords[0]) + 1, 200)
    y_points = m*x_range + c
    
    # out_array = np.zeros(bound_array.shape, dtype=np.int32)
    # out_array[x_range, :, np.ceil(y_points).astype(np.int32)] = 1
    # out_array = load_mapping.remove_bounding_box(out_array, bounding_box, masked_array.shape)
    # return out_array, y
    return [np.ceil(x_range).astype(np.int32), np.ceil(y_points).astype(np.int32)], y
    

def get_coronal_split_1(masked_array):
    bounding_box = load_mapping.get_bounding_box(masked_array)
    bound_array = masked_array[bounding_box]
    com = center_of_mass(bound_array)
    cor, sag = int(com[0]), int(com[1])
    
    bound_coords = np.nonzero(bound_array[cor,:,:]) # [sag, ax]
    bound_coords_left = bound_coords[1][np.where(bound_coords[0] < sag)]
    bound_coords_right = bound_coords[1][np.where(bound_coords[0] >= sag)]
    mid_left = (max(bound_coords_left) + min(bound_coords_left))/2
    mid_right = (max(bound_coords_right) + min(bound_coords_right))/2
    
    m, c = np.linalg.lstsq(np.vstack([np.array([0, max(bound_coords[0])]), np.ones(2)]).T, 
                           np.array([mid_left, mid_right]), rcond=None)[0]
    sag_range = np.arange(0, np.max(bound_coords[0]))
    ax_points = m*sag_range + c
    
    out_array = np.zeros(bound_array.shape, dtype=np.int32)
    out_array[:, sag_range, np.ceil(ax_points).astype(np.int32)] = 1
    out_array = load_mapping.remove_bounding_box(out_array, bounding_box, masked_array.shape)
    return out_array, cor


# roi = 32
for roi in range(32, 33):
    z_array = np.logical_and(mask.bound == roi, mask.spine == 1)

    # plot_all_planes(z_array, mask.slice_thickness, mask.pixel_spacing, [*com], axes=1)
    
    z_dim = load_mapping.get_bounding_box(z_array)
    z_arr_bound = z_array[z_dim]

    z_test_new_func, zy = get_coronal_split(z_array)
    z_test_new_func_1, zy_1 = get_coronal_split_1(z_array)
    # zzzz = z_test_new_func[z_dim]
    # zzzz1 = z_test_new_func_1[z_dim]
    plt.imshow(z_arr_bound[:,zy,:].T)
    plt.plot(z_test_new_func[0], z_test_new_func[1])
    
    zz_mask = np.zeros(z_arr_bound.shape, dtype=bool)
    zz_mask[z_test_new_func[0], :, z_test_new_func[1]] = 1
    
    # for zx,zy in zip(z_test_new_func[0], z_test_new_func[1]):
    #     zz_mask[:zx+1, :, :zy+1] = 1
    #     zz_mask[:zx+1, :, :zy+1] = 1
    
    plt.imshow(zz_mask[:,zy, :].T)
    break
    plt.figure()
    plt.subplot(1,2,1)
    plt.title(roi)
    plt.imshow(z_arr_bound[:,zy,:].T)
    plt.imshow(zzzz[:,zy,:].T, alpha=zzzz[:,zy,:].T, cmap='hot')

    plt.subplot(1,2,2)
    plt.title(roi)
    plt.imshow(z_arr_bound[zy_1,:,:].T)
    plt.imshow(zzzz1[zy_1,:,:].T, alpha=zzzz1[zy_1,:,:].T, cmap='hot')
    

    plt.show()

#%%



#%%

