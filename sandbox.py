# -*- coding: utf-8 -*-
"""
Created on Wed Jan 31 12:04:12 2024

@author: runep
"""

from skimage import measure
from skimage.feature import canny
from copy import deepcopy

plt.imshow(measure.label(obj_seg.A[sag_bound], connectivity=2).T)
plt.imshow(measure.label(obj_seg.A[cor_bound], connectivity=2).T)
#%%
zz_test = np.zeros((46,25))
for i in range(final_bound[1].start, final_bound[1].stop):
    # zzzz = tuple((slice(287,347), i, slice(181,215)))
    zzzz = tuple((slice(final_bound[0].start, final_bound[0].stop), 
                  i,
                  slice(final_bound[2].start, final_bound[2].stop)
                  ))
    zz_labels = measure.label(obj_seg.A[zzzz], connectivity=2)
    zz_test += zz_labels
    plt.figure()
    plt.title(i)
    plt.imshow(obj_ct.A[zzzz].T, cmap='gray')
    plt.imshow(zz_labels.T, alpha=.3*zz_labels.T)
    plt.axis('off')
    plt.show()

#%%
for i in range(231, 274):
    plt.figure()
    plt.axvline(init_guess[0], ymin=0, ymax=20, color='k')
    plt.axhline(0,0,40, color='k')
    
    plt.plot(np.sum(obj_seg.A[285:331, i, 181:206], axis=1))
    plt.plot(np.gradient(np.sum(obj_seg.A[285:331, i, 181:206], axis=1)))
    plt.plot(np.mean(np.sum(obj_seg.A[285:331, :, 181:206], axis=1), axis=1))
    plt.plot(np.gradient(np.mean(np.sum(obj_seg.A[285:331, :, 181:206], axis=1), axis=1)))
    plt.show()
#%%

regions_mask = np.zeros(final_dim)
# regions_mask = np.zeros((512,512,400))

cor_mean = np.mean(np.sum(obj_seg.A[final_bound], axis=0), axis=1) # [43,25]
sag_mean = np.mean(np.sum(obj_seg.A[final_bound], axis=1), axis=0) # [46,25]
ax_sum = np.sum(obj_seg.A[final_bound], axis=2) # [46,43]
ax_gradients = [np.gradient(ax_sum[:,i]) for i in range(ax_sum.shape[1])]

# Different method could just take half of the shape / bounding box dimension
ax_split = np.where(sag_mean > np.percentile(sag_mean, 50))[0]
ax_split = np.diff(ax_split[[0,-1]])[0]//2 + ax_split[0]

sag_split = np.where(cor_mean > np.percentile(cor_mean, 50))[0]
sag_split = np.diff(sag_split[[0,-1]])[0]//2 + sag_split[0]

def get_coronal_indices(gradient_list, dimensions, com):
    process = np.zeros(dimensions)
    com = int(np.ceil(com[0]))
    
    for i in range(len(gradient_list)):
        sub_gradient = gradient_list[i][com:]
        
        cor_split = 0
        for j in range(len(sub_gradient)):
            if (j < 4) & (sub_gradient[j] >= 0):
                continue
            elif sub_gradient[j] >= 0:
                cor_split = com + j
                break
        process[cor_split:, i, :] = 1
    return process

process = get_coronal_indices(ax_gradients, final_dim, init_guess)
# process_big = np.zeros((512,512,400))

# process_big[final_bound[0].start:final_bound[0].stop,
#             final_bound[1].start:final_bound[1].stop,
#             final_bound[2].start:final_bound[2].stop] = process

# plt.imshow(zz_test[:,256,:].T)
# plt.imshow(z_slice_final[:,15,:].T)

zz_test = deepcopy(z_slice_final)
zz_test[process == 1] = 0

ax_mean = np.mean(np.sum(zz_test, axis=2), axis=1) # [46,43]
cor_split = np.where(ax_mean > np.percentile(ax_mean, 50))[0]
cor_split = np.diff(cor_split[[0,-1]])[0]//2 + cor_split[0]




regions_mask[z_slice_final == 1] = 1
regions_mask[:cor_split,sag_split:,:ax_split] *= 1
regions_mask[:cor_split,:sag_split,:ax_split] *= 2
regions_mask[cor_split:,:sag_split,:ax_split] *= 3
regions_mask[cor_split:,sag_split:,:ax_split] *= 4
regions_mask[:cor_split,sag_split:,ax_split:] *= 5
regions_mask[:cor_split,:sag_split,ax_split:] *= 6
regions_mask[cor_split:,:sag_split,ax_split:] *= 7
regions_mask[cor_split:,sag_split:,ax_split:] *= 8
regions_mask[np.logical_and(process, z_slice_final)] = 9

regions_big = np.zeros((512,512,400))
regions_big[final_bound[0].start:final_bound[0].stop,
            final_bound[1].start:final_bound[1].stop,
            final_bound[2].start:final_bound[2].stop] = regions_mask

#%%
cor_pos = 20
sag_pos = 30
ax_pos  = 10
# cor_pos = cor_ind
# sag_pos = sag_ind
# ax_pos = ax_ind

fig, axs = plt.subplots(1,3)
axs[0].imshow(regions_mask[cor_pos,...].T, alpha = .1*regions_mask[cor_pos,...].T)
axs[0].plot(sag_pos, ax_pos, '+r', markersize=15)
axs[0].set_title('Coronal')
axs[0].set_xlabel('Sagittal')
axs[0].set_ylabel('Axial')


axs[1].imshow(regions_mask[:,sag_pos,:].T, alpha = .1*regions_mask[cor_pos,...].T)
axs[1].plot(cor_pos, ax_pos, '+r', markersize=15)
axs[1].set_title('Sagittal')
axs[1].set_xlabel('Coronal')
axs[1].set_ylabel('Axial')

# plt.imshow(regions_big[:,246,:].T)
axs[2].imshow(obj_ct.A[...,ax_ind], cmap='gray')
axs[2].imshow(regions_big[...,ax_ind], alpha = .1*regions_big[...,ax_ind])
# axs[2].plot(sag_pos, cor_pos, '+r', markersize=15)
axs[2].set_title('Axial')
axs[2].set_xlabel('Sagittal')
axs[2].set_ylabel('Coronal')


axs[0].set_aspect(obj_ct.slice_thickness/obj_ct.pixel_spacing)
axs[1].set_aspect(obj_ct.slice_thickness/obj_ct.pixel_spacing)


# for ax in axs.ravel():
    # ax.axis('off')
plt.show()


img = nib.Nifti1Image(np.moveaxis(np.flip(regions_big,[0,1,2]),0,1), obj_seg.files['L1'].affine)
nib.save(img, join(sub, 'test_regions_mask.nii.gz'))
#%%

zz_locs = np.nonzero(obj_seg.A)
zz_ex_pts = [(np.min(zz_locs[0]),np.max(zz_locs[0])), 
             (np.min(zz_locs[1]),np.max(zz_locs[1])),
             (np.min(zz_locs[2]),np.max(zz_locs[2]))]

zz_bounding_box = np.s_[zz_ex_pts[0][0]:zz_ex_pts[0][1],
                        zz_ex_pts[1][0]:zz_ex_pts[1][1],
                        zz_ex_pts[2][0]:zz_ex_pts[2][1]]
zz_cor = (np.min(zz_locs[0]), np.max(zz_locs[0]))
zz_sag = (np.min(zz_locs[1]), np.max(zz_locs[1]))
zz_ax = (np.min(zz_locs[2]), np.max(zz_locs[2]))
#%% Coronal sums over axial axis
zz_vertebra = np.zeros((46,43,25))

# zz_sag_sums = np.sum(obj_seg.A[285:331, 231:274, 181:206], axis=2)
# zz_sag_sum_gradients = [np.gradient(zz_sag_sums[:,i]) for i in range(zz_sag_sums.shape[1])]
# for j in range(43):
    # zz_sag_sub_gradients = zz_sag_sum_gradients[j][int(np.ceil(init_guess[0])):]
    
    # zz_loc = 0
    
    # for i in range(len(zz_sag_sub_gradients)):
    #     if (i < 6) & (zz_sag_sub_gradients[i] >= 0):
    #         # print(i)
    #         continue
    #     elif zz_sag_sub_gradients[i] >= 0:
    #         # print(i)
    #         zz_loc = int(np.ceil(init_guess[0])) + i
    #         break
    # plt.figure()
    # plt.axvline(init_guess[0], ymin=0, ymax=20, color='k')
    # plt.axvline(zz_loc, ymin=0, ymax=20, color='k')
    # plt.axhline(0,0,40, color='k')
    # plt.plot(np.mean(np.sum(obj_seg.A[285:331, :, 181:206], axis=1), axis=1))
    # plt.plot(zz_sag_sums[:,j])
    # plt.plot(zz_sag_sum_gradients[j])
    # plt.show()
    # zz_vertebra[zz_loc:, j, :] = 1
    