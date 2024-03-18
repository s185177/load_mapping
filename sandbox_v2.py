# -*- coding: utf-8 -*-
"""
Created on Wed Mar  6 11:16:18 2024

@author: runep
"""


#%%
from scipy.optimize import curve_fit
from skimage.feature import corner_harris, corner_peaks
z = np.sum(vert_samp[vert_bound], axis=1)    

z_coords = corner_peaks(corner_harris(z), min_distance=5)

n=0
for i in range(2,33):
    z = vert_samp[vert_bound][i,:,:]
    z_coords = corner_peaks(corner_harris(z), min_distance=4, threshold_rel=0.08, num_peaks=4)    
    if len(z_coords) < 4:
        n += 1
        print("Skipping %", i, n)
        continue
    elif i == 14:
        break
    print(i, z_coords)
    zr = z_coords[1,1] - (z_coords[1,1] - z_coords[0,1])//2
    zl = z_coords[3,1] - (z_coords[3,1] - z_coords[2,1])//2
    
    plt.figure()
    plt.title(i)
    plt.imshow(z.T)
    plt.plot(z_coords[:,0], z_coords[:,1], 'b.')
    plt.plot([min(z_coords[:,0]),max(z_coords[:,0])],[zl,zr], 'r')
    
    plt.show()

#%%

def z_func(xy, a, b, c, d, e, f): 
    x, y = xy 
    return a + b*x + c*y + d*x**2 + e*y**2 + f*x*y


zx,zy = np.meshgrid(z_coords[:,0], z_coords[:,1])

popt,_ = curve_fit(z_func, (zx,zy), z)

fig = plt.figure(figsize=(10, 6))
ax = fig.add_subplot(111, projection='3d')
X1, X2 = np.mgrid[:27, :25]
ax.plot_surface(X1, X2, z_samp[...,10].T, alpha=.2, label='Data Points', cmap='jet')

ax.plot_surface(xx, yy, zz, alpha=0.5, label='Fitted Plane', color='b')
#%% Plot different spine fitting methods
# So far splprep show most promising results!

# x_new = np.linspace(0,split_sub.shape[-1]-1, split_sub.shape[-1]*2)
# cm = 'nipy_spectral'
# al = .9

# fig, axs = plt.subplots(1,2)
# fit_sag = np.polyfit(x_list_sag, y_list_sag, deg=4)
# p_sag = np.poly1d(fit_sag)
# axs[0].imshow(cor_sum.T, cmap=cm)
# axs[0].plot(p_sag(x_new),x_new, alpha=al)
# axs[0].set_xlim([0,314])

# fit_cor = np.polyfit(x_list_cor, y_list_cor, deg=4)
# p_cor = np.poly1d(fit_cor)
# axs[1].imshow(sag_sum.T, cmap=cm)
# axs[1].plot(p_cor(x_new),x_new, 'r', alpha=al)
# axs[1].set_xlim([0,175])
# plt.tight_layout()
# plt.show()

# fig, axs = plt.subplots(1,2)
# tck,u = splprep([x_list_sag, y_list_sag], s=1e3)
# new_points = splev(u,tck)
# axs[0].imshow(cor_sum.T, cmap=cm)
# axs[0].plot(new_points[1], new_points[0], alpha=al)

# tck,u = splprep([x_list_cor, y_list_cor], s=1e3)
# new_points = splev(u,tck)
# axs[1].imshow(sag_sum.T, cmap=cm)
# axs[1].plot(new_points[1], new_points[0], 'r', alpha=al)
# plt.tight_layout()
# plt.show()

# fig, axs = plt.subplots(1,2)
# spl = splrep(x_list_sag, y_list_sag, s=1e3)
# y_new = splev(x_new, spl)
# axs[0].imshow(cor_sum.T, cmap=cm)
# axs[0].plot(y_new, x_new, alpha=al)

# spl = splrep(x_list_cor, y_list_cor, s=1e3)
# y_new = splev(x_new, spl)
# spl = BSpline(*spl)
# axs[1].imshow(sag_sum.T, cmap=cm)
# axs[1].plot(spl(x_new), x_new, 'r', alpha=al)
# axs[1].set_xlim([0,175])
# plt.tight_layout()
# plt.show()

#%% Sanity check of fit spine

cor_sum = np.sum(np.where(mask.split_coords == 1, 1,0), axis=0)
sag_sum = np.sum(np.where(mask.split_coords == 2, 1,0), axis=1)

cm, al = 'nipy_spectral', .9
fig, axs = plt.subplots(1,2)
axs[0].set_title('Cor_sum')
axs[0].imshow(cor_sum.T, cmap=cm)
axs[0].plot(y_list_sag, x_list_sag, alpha=al)
axs[0].set_xlabel('Sagittal')
axs[0].set_ylabel('Axial')

axs[1].set_title('Sag_sum')
axs[1].imshow(sag_sum.T, cmap=cm)
axs[1].plot(y_list_cor, x_list_cor, 'r', alpha=al)
axs[1].set_xlabel('Coronal')
axs[1].set_ylabel('Axial')
plt.tight_layout()
plt.show()


del split_sag, split_cor

#%% How to access spine_map
# roi = (load_map._ANTERIOR_KEY)
# roi = (load_map._ANTERIOR_KEY|load_map._RIGHT_KEY)
roi = (load_map._VERTEBRAE_KEY|load_map._ANTERIOR_KEY|load_map._RIGHT_KEY)
z_test = np.array(np.bitwise_and(mask.spine, roi) == roi)
zz_test = np.where(mask.bound == 34, mask.bound, 0)
zzz_test = np.array(np.logical_and(zz_test, z_test))


plot_all_planes(z_test, mask.slice_thickness, mask.pixel_spacing, [130, 156, 175], axes=1)
plt.show()
# plot_all_planes(zz_test, mask.slice_thickness, mask.pixel_spacing, [130, 156, 175], axes=1)
# plt.show()
# plot_all_planes(zzz_test, mask.slice_thickness, mask.pixel_spacing, [130, 156, 175], axes=1)
# plt.show()
