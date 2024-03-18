# -*- coding: utf-8 -*-
"""
Created on Mon Oct 23 10:13:51 2023

@author: runep
"""

#%%
import SimpleITK as sitk

ref_img = sitk.GetImageFromArray(pet_img)
ref_dims = ref_img.GetSize()
samp_img = sitk.GetImageFromArray(ct_img)

#%%
sitk.AffineTransform(3)
euler3d = sitk.Euler3DTransform()
euler3d.SetCenter(ref_dims)
euler3d.SetTranslation(ref_dims)
# inv_euler3d = euler3d.GetInverse()

# extreme_points = [(0.0,0.0,0.0),
#                   (0.0,0.0,0.0),
#                   (0.0,0.0,0.0),
#                   (0.0,0.0,0.0),
#                   (0.0,0.0,0.0),
#                   (0.0,0.0,0.0),
#                   (0.0,0.0,0.0),
#                   (0.0,0.0,0.0)]
# extreme_points_transform = [inv_euler3d.TransformPoint(pnt) for pnt in extreme_points]



resampled = sitk.Resample(ref_img, euler3d)