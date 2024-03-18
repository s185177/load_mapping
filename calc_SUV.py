# -*- coding: utf-8 -*-
"""
Created on Mon Feb  5 17:57:01 2024

@author: runep
"""

import os
import pydicom
import numpy as np
import nibabel as nib
import SimpleITK as sitk
import matplotlib.pyplot as plt

from os.path import join
# from load_mapping import Image, get_bounding_box, plot_all_planes
import load_mapping

# dir_path = '/media/klinfys/KINGSTON/TotalSegmentator/Spine/S03/01'
dir_path = r'E:\TotalSegmentator\Spine\S03\02'


files_PET = load_mapping.dir_walker(dir_path, 'PET/WB')
PET = load_mapping.resample(join(dir_path, 'Segmentations.nii'), files_PET)

PET = np.moveaxis(np.flip(PET,[0,1,2]), 0, 2)

Spine = load_mapping.Image(join(dir_path, 'Segmentations.nii'))
Spine.A = Spine.get_data()
Spine.A.shape

bounding_box = load_mapping.get_bounding_box(Spine.A)
PET_bound = PET[bounding_box]


roi = 27
masked = PET * np.where(Spine.A == roi, 1,0)

vert_bound = load_mapping.get_bounding_box(masked)
bound_dim = masked[vert_bound].shape
#%%
load_mapping.plot_all_planes(PET, Spine.slice_thickness, Spine.pixel_spacing, [275,260,120], axes=1)
plt.show()
load_mapping.plot_all_planes(Spine.A, Spine.slice_thickness, Spine.pixel_spacing, [275,260,120], axes=1)
plt.show()
load_mapping.plot_all_planes(masked, Spine.slice_thickness, Spine.pixel_spacing, [275,260,120], axes=1)
plt.show()
load_mapping.plot_all_planes(masked[vert_bound], Spine.slice_thickness, Spine.pixel_spacing, [10,19,10], axes=1)
plt.show()
#%%

PET_header = pydicom.dcmread(files_PET[0])

radio_pharmacy = PET_header[0x00540016][0]
# time_injection = float(radio_pharmacy[0x00181072].value)
# time_acquisition = float(PET_header.AcquisitionTime)
# half_life = radio_pharmacy[0x00181075].value / 60 # min
dose_total = radio_pharmacy[0x00181074].value #MBq
patient_weight = PET_header.PatientWeight.real # kg

# PET_header.RadiopharmaceuticalInformationSequence[0].RadionuclideTotalDose
# time_delta = (time_acquisition - time_injection) / 100 # min
# dose_corrected = dose_total * np.exp( -time_delta * np.log(2) / half_life) # Bq
# SUV_factor = (PET_header.RescaleSlope * PET_header.PatientWeight.real * 1e3) / dose_corrected # [g/Bq] = [] * [Kg]* 1000[g/kg] / [Bq]
SUV_factor = dose_total /  patient_weight # MBq/kg
#%%

PET_SUV = masked[vert_bound] / SUV_factor # MBq/mL / 
load_mapping.plot_all_planes(PET_SUV, Spine.slice_thickness, Spine.pixel_spacing, [10,19,10], axes=1)
plt.show()
# np.mean(PET_SUV)
#%%
img = load_mapping.remove_bounding_box(PET_SUV, vert_bound, Spine.A.shape).astype(np.float32)
nib.save(nib.Nifti1Image(np.moveaxis(np.flip(img,[0,1,2]),0,1), Spine.file.affine), 'data/PET_SUV_S03_02' + '.nii.gz')
#%%

# load_mapping.plot_all_planes(PET_SUV * np.where(PET_SUV > 20,1,0), Spine.slice_thickness, Spine.pixel_spacing, [10,19,10], axes=1)
# plt.show()

# plt.imshow(np.mean(PET_SUV, axis=2), cmap='nipy_spectral')
