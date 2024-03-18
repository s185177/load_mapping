# -*- coding: utf-8 -*-
"""
Created on Tue Feb 13 18:06:21 2024

@author: runep
"""
import os
import nibabel as nib
import SimpleITK as sitk
import matplotlib.pyplot as plt
from tkinter.filedialog import askdirectory, askopenfilename


dir_path = r'C:\Users\runep\OneDrive\Skrivebord\Work\total_segmentator\Feet\S01'
#%%

ct = sitk.ReadImage(sitk.ImageSeriesReader_GetGDCMSeriesFileNames(r'C:\Users\runep\OneDrive\Skrivebord\Work\total_segmentator\Feet\S01\CT'))
pet = sitk.ReadImage(r'C:\Users\runep\OneDrive\Skrivebord\Work\total_segmentator\Feet\S01\PET\BoneWBTomo_TX_IRACRR_FAME_PASTED_Tc99m001_DS.dcm')

pet_resamp = sitk.Resample(pet, ct)


nib_img = nib.load(r'C:\Users\runep\OneDrive\Skrivebord\Work\total_segmentator\Feet\S01\segmentations\tibia.nii.gz').get_fdata()
image_out_new = nib.orientations.apply_orientation(nib_img, nib.orientations.axcodes2ornt(('LPS')))

#%%
os.chdir(askdirectory(title='Choose data directory'))
file_ref = askopenfilename(initialdir=os.getcwd(), title='Select reference NIFTI', filetypes=[('NIfTI','.gz')])
file_in = askopenfilename(initialdir=os.getcwd(), title='Select input DICOM', filetypes=[('DICOM','.dcm')])

class sitk_io:
    def __init__(self, image_io : str, path : str):
        self.image  = sitk.ReadImage(path, imageIO=image_io)
        # reader.SetImageIO(image_io)
        # reader.SetFileName(path)
        # self.image = reader.Execute()
    
    def write(self, file_name : str):
        sitk.WriteImage(self.resamp, file_name)
        
z_ref = sitk_io("NiftiImageIO", file_ref)
z_in = sitk_io("", file_in)

z_ref.image_ornt = sitk.DICOMOrient(z_ref.image, 'LPS')
z_in.resamp = sitk.Resample(z_in.image, z_ref.image_ornt)
z_in.write('test.nii.gz')
#%%


dir_ref = askdirectory()
series_IDs = sitk.ImageSeriesReader.GetGDCMSeriesIDs(dir_ref)
series_file_names = sitk.ImageSeriesReader.GetGDCMSeriesFileNames(dir_ref, series_IDs[0])

series_reader = sitk.ImageSeriesReader()
series_reader.SetFileNames(series_file_names)
series_reader.MetaDataDictionaryArrayUpdateOn()
series_reader.LoadPrivateTagsOn()
image3D = series_reader.Execute()

for k in series_reader.GetMetaDataKeys(0):
       v = series_reader.GetMetaData(k)
       print(f'({k}) = = "{v}"')
