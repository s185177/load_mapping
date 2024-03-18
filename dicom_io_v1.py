#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Feb  2 12:12:25 2024

@author: klinfys
"""

from os.path import join
import numpy as np
import os
from os.path import join
# import nibabel as nib
# import matplotlib.pyplot as plt
import pydicom
import load_mapping as l_map

global scanner
global CT_series_name

class DICOM:
    def __init__(self, files):
        self.keys = list(files.keys())
        self.files = files
        self.study = files[self.keys[0]].StudyDescription
        if '_8' in self.study:
            DICOM.scanner = 8
        elif '_10' in self.study:
            DICOM.scanner = 10
        else:
            DICOM.scanner = None
    
    def get_data(sorted_files):
        print(len(sorted_files))
        data_array = np.array([ds for ds in sorted_files], dtype=np.float64)
        return data_array
        

class CT(DICOM):
    def __init__(self, files):
        super().__init__(files)
        if self.scanner == 8:
            self.series_description = ''
        elif self.scanner == 10:
            self.series_description = 'AC_ULDCT 2mm'
        
        self.files = {key:value for key,value in self.files.items() if self.series_description in value.SeriesDescription}
        self.series_numbers = sorted(list(set([f.SeriesNumber for f in self.files.values()])))
        self.scan_baseline = sorted([f for f in self.files.values() if f.SeriesNumber == self.series_numbers[0]],
                                    key=lambda x: x.InstanceNumber)
        
        self.scan_intervention = sorted([f for f in self.files.values() if f.SeriesNumber == self.series_numbers[1]], 
                                        key= lambda x: x.InstanceNumber)
        
class PET(DICOM):
    def __init__(self, files):
        super().__init__(files)
    # pass
#%%
subject = r'E:\TotalSegmentator\Spine\S02'
files = l_map.dir_walker(subject, 'all_dicom')

CT_files = {}
CT_series = []
CT_series_num = []
CT_acq_num = []

PET_files = {}
PET_series = []
PET_series_num = []
PET_acq_num = []

tot_files = len(files)
for n, f in enumerate(files):
    if n % 1000 == 0:
        print("{:d} loaded out of {:d}".format(n,tot_files))
    ds = pydicom.dcmread(f)
    
    if ds[0x00080060].value == 'CT':
        CT_files[n] = ds
        CT_series.append(ds.SeriesDescription)
        CT_series_num.append(ds.SeriesNumber)
        # CT_acq_num.append(ds.AcquisitionNumber)
        
    elif ds[0x00080060].value == 'PT':
        PET_files[n] = ds
        PET_series.append(ds.SeriesDescription)
        PET_series_num.append(ds.SeriesNumber)
        # PET_acq_num.append(ds.AcquisitionNumber)
#%%
CT = CT(CT_files)
CT.series_numbers

len(CT.scan_baseline)
len(CT.scan_intervention)

DICOM.get_data(CT.scan_baseline)
DICOM.get_data()

CT.scan_baseline[0].pixel_array
CT.A_baseline = CT.get_data(CT.scan_baseline)

PET = PET(PET_files)


#%%
n = 0
for key, value in CT.files.items():
    n+=1
    print(value.SeriesDescription)
    if n == 100:
        break
    

#%%
def get_correct_series_num(files_series_name, files_series_num, series_name):
    num_list = []
    indeces = []
    for i in range(len(files_series_name)):
        if series_name in files_series_name[i]:
            num_list.append(files_series_num[i])
            indeces.append(i)
    return sorted(list(set(num_list))), indeces

unique_series, series_indeces = get_correct_series_num(CT_series, CT_series_num, 'AC_ULDCT')
#%%

        
#%%
CT_series_uniq = list(set(CT_series))
PET_series_uniq = list(set(PET_series))

CT_keys = list(CT_files.keys())
CT_sorted = []

for i in range(len(CT_files)):
    if 'AC_ULDCT' in CT_series[i]:
        
        


#%%
indeces = CT_files.keys()
# n = 0
asd = []
# CT_sorted = []
for f in indeces:
    if CT_files[f].SeriesDescription == CT_series[3]:
        print(CT_files[f].SeriesNumber)
        asd.append(CT_files[f].SeriesNumber)
        # CT_sorted.append(CT_files[f])
#%%
files[4]
files[5306]

#%%
# os.system("TotalSegmentator -i load_mapping/data/S01/AC_CT_3_2.nii.gz -o test_1 --roi_subset vertebrae_L1 --device gpu --fast")
# os.system("flirt \
#           -in load_mapping/data/S01/PET_WB_NaF_TrueX_+_TOF_2,5_mm_4_3001.nii.gz  \
#           -ref load_mapping/data/S01/AC_CT_3_2.nii.gz \
#           -out flirt.nii.gz \
#           -applyxfm \
#           -usesqform")

# file = nib.load('data/S01/test_1/vertebrae_body.nii.gz')

# data = file.get_fdata()
# np.unique(data)


# slices = np.unique(np.nonzero(data)[-1])
# #%%
# for i in slices:
#     plt.figure()
#     plt.imshow(data[...,i], cmap='gray')
#     plt.axis('off')
#     plt.show()

# from sys import argv
# print(argv[1:])
#%%
# import os
# import subprocess
# os.chdir(r'C:\Users\runep\OneDrive - Danmarks Tekniske Universitet\DTU\Thesis\Code\load_mapping')

# mnt_path = '/mnt/c/Users/runep/OneDrive - Danmarks Tekniske Universitet/DTU/Thesis/Code/load_mapping'


# p = subprocess.run('wsl /home/rped0143/fsl/share/fsl/bin/flirt -in', capture_output=True, shell=True)
# print(p)
# print(p.args)
# print(p.returncode)
# print(p.stdout)
#%%
# import time

ct = []
pet = []
rest = []
path = r'D:\TotalSegmentator\Spine\S02\01\CT\DICOM\24020613'
# path = r'D:\TotalSegmentator\Feet'
for root, sub, subsub in os.walk(path):
    if len(subsub) > 1:
        for n, i in enumerate(subsub):
            file = pydicom.dcmread(join(root,i), force=True)
            if all([tag in file._dict.keys() for tag in CTImageStorage_required_tags]):
                ct.append(file)
                print("CT stored\n")
            elif all([tag in file._dict.keys() for tag in NuclearMedicineImageStorage_required_tags]):
                pet.append(file)
                print("PET stored\n")
            else:
                rest.append(join(root,i))

del root, sub, subsub, i
#%%

# files_dcm = [pydicom.dcmread(f) for f in files]

for i in range(20):
    dict_keys = pydicom.dcmread(files[i])._dict.keys()
    passed = all([tag in dict_keys for tag in CTImageStorage_required_tags])
    passed_1 = all([tag in dict_keys for tag in NuclearMedicineImageStorage_required_tags])
    print(passed)
    print(passed_1)
    print()
    
# z = files_dcm[0]._dict.keys()
# files_dcm[0].SeriesDescription
# series = [f.SliceLocation for f in files_dcm]

#%%

Patient_required_tags: list[int] = [
  0x00100010, # PatientName
  0x00100020, # PatientID
  0x00100030, # PatientsBirthDate
  0x00100040, # PatientSex
]
General_Study_required_tags: list[int] = [
  0x00080020, # StudyDate
  0x00080030, # StudyTime
  0x00080050, # AccessionNumber
  0x00080090, # ReferringPhysicianName
  0x0020000D, # StudyInstanceUID
  0x00200010, # StudyID
]
General_Series_required_tags: list[int] = [
  0x00080060, # Modality
  0x0020000E, # SeriesInstanceUID
  0x00200011, # SeriesNumber
]
Frame_of_Reference_required_tags: list[int] = [
   0x00200052, # FrameOfReferenceUID
   0x00201040  # PositionReferenceIndicator
]
General_Equipment_required_tags: list[int] = [
  0x00080070, # Manufacturer
]
General_Image_required_tags: list[int] = [
  0x00200013 # InstanceNumber
]
Image_Plane_required_tags: list[int] = [
  0x00180050, # Slice thickness
  0x00200032, # ImagePosition
  0x00200037,  # ImageOrientation
  0x00280030  # PixelSpacing
]
Image_Pixel_required_tags: list[int] = [
  0x00280002, # Samples per Pixel
  0x00280004, # Photometric Interpretation
  0x00280010, # Rows
  0x00280011, # Columns
  0x00280100, # BitsAllocated
  0x00280101, # BitsStored,
  0x00280102, # HighBit
  0x00280103, # PixelRepresentation
  0x7FE00010, # PixelData
]
CT_Image_required_tags: list[int] = [
  0x00080008, # ImageType
  0x00180060, # KVP
  0x00280100, # BitsAllocated
  0x00280101, # BitsStored,
  0x00280102, # HighBit
  0x00281052, # RescaleIntercept
  0x00281053, # RescaleType
]
SOP_Common_required_tags: list[int] = [
  0x00080016, # SOPInstanceUID
  0x00080018  # SOPClassUID
]

CTImageStorage_required_tags: list[int] = Patient_required_tags \
  + General_Study_required_tags \
  + General_Series_required_tags \
  + Frame_of_Reference_required_tags \
  + General_Equipment_required_tags \
  + General_Image_required_tags \
  + Image_Plane_required_tags \
  + Image_Pixel_required_tags \
  + CT_Image_required_tags \
  + SOP_Common_required_tags
  
NM_PET_Patient_Orientation_required_tags: list[int] = [
  0x00540410, # PatientOrientationCodeSequence
  0x00540414  # PatientGantryRelationshipCodeSequence
]  
NM_Image_Pixel_required_tags: list[int] = [
  0x00280002, # Samples per Pixel
  0x00280004, # Photometric Interpretation
  0x00280030, # PixelSpacing
  0x00280100, # BitsAllocated
  0x00280101, # BitsStored,
  0x00280102, # HighBit
]
MultiFrame_required_tags = [
  0x00280008, # NumberOfFrames
  0x00280009, # FrameIncrementPointer
]
NM_Image_required_tags: list[int] = [
  0x00080008 # ImageType
]
NM_Isotope_required_tags: list[int] = [
  0x00540012, # EnergyWindowInformationSequence
  0x00540016, # RadiopharmaceuticalInformationSequence
]

NM_Detector_required_tags: list[int] = [
  0x00540022,  # DetectorInformationSequence
]
NuclearMedicineImageStorage_required_tags: list[int] = Patient_required_tags \
  + General_Study_required_tags \
  + General_Series_required_tags \
  + NM_PET_Patient_Orientation_required_tags \
  + Frame_of_Reference_required_tags \
  + General_Equipment_required_tags \
  + General_Image_required_tags \
  + Image_Pixel_required_tags \
  + NM_Image_Pixel_required_tags \
  + MultiFrame_required_tags \
  + NM_Image_required_tags \
  + NM_Isotope_required_tags \
  + NM_Detector_required_tags \
  + SOP_Common_required_tags
  
del Patient_required_tags, General_Study_required_tags, General_Series_required_tags, \
    Frame_of_Reference_required_tags, General_Equipment_required_tags, General_Image_required_tags, \
    Image_Plane_required_tags, Image_Pixel_required_tags, CT_Image_required_tags, SOP_Common_required_tags, \
    NM_PET_Patient_Orientation_required_tags, MultiFrame_required_tags, NM_Image_Pixel_required_tags, NM_Detector_required_tags, \
        NM_Image_required_tags, NM_Isotope_required_tags

