# -*- coding: utf-8 -*-
"""
Created on Wed Mar  6 17:33:21 2024

@author: runep
"""

import pydicom
import os
import numpy as np
import nibabel as nib
import SimpleITK as sitk
import matplotlib.pyplot as plt
from scipy import ndimage
from os.path import join

def directory_walker(modality):
    files = []
    for x,_,z in os.walk(join(subject, modality)):
        if (len(z) > 0) & ('DICOMDIR' not in z):
            sub_files = [join(x,f) for f in z]
            files += sub_files
    return files


subject = r'D:\TotalSegmentator\Spine\S02\01'

files_CT = directory_walker('CT')
files_PET = directory_walker('PET')


series_reader = sitk.ImageSeriesReader()
series_reader.SetFileNames(files_PET)
series_reader.MetaDataDictionaryArrayUpdateOn()
series_reader.LoadPrivateTagsOn()
image3D = series_reader.Execute()

z_test = sitk.ReadImage(join(subject, 'Segmentations.nii'), imageIO='NiftiImageIO')


# nib_img = nib.load(join(subject, 'Segmentations.nii')).get_fdata()
# image_out_new = nib.orientations.apply_orientation(nib_img, nib.orientations.axcodes2ornt(('LPS')))

z_test1 = sitk.DICOMOrient(z_test, 'LPS')
z_test1 = sitk.Resample(image3D, z_test1)

sitk.GetArrayFromImage(z_test1).shape
plt.imshow(nib_img[...,100])
plt.imshow(image_out_new[...,100])
# PET = sitk.ReadImage(sitk.ImageSeriesReader_GetGDCMSeriesFileNames(join(subject,'PET/DICOM')))
#%%
CT = pydicom.dcmread(files_CT[0])
CT_imgs = np.empty((CT.Rows, CT.Columns,len(files_CT)), dtype=np.float64)
PET = pydicom.dcmread(files_PET[0])
PET_imgs = np.empty((PET.Rows, PET.Columns, len(files_PET)), dtype=np.float64)

# for i in range(len(files_CT)):
    # CT_imgs[...,i] = pydicom.dcmread(files_CT[i]).pixel_array

for i in range(len(files_PET)):
    PET_imgs[...,i] = pydicom.dcmread(files_PET[i]).pixel_array


