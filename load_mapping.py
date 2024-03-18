# -*- coding: utf-8 -*-
"""
Created on Mon Feb  5 18:16:22 2024

@author: runep
"""
import re
import os
import numpy as np
import nibabel as nib
import SimpleITK as sitk
import matplotlib.pyplot as plt
from os.path import join
import matplotlib.pylab as pylab
params = {'legend.fontsize': 'x-large',
          'figure.figsize': (16, 8),
          'axes.labelsize': 'xx-large',
          'axes.titlesize':'x-large',
          'xtick.labelsize':'x-large',
          'ytick.labelsize':'x-large'}
pylab.rcParams.update(params)
del params

class Image:
    def init_image(self, path):
        self.file = nib.load(path)
        self.header = self.file.header
        self.pixel_spacing = self.header['pixdim'][1]
        self.slice_thickness = self.header['pixdim'][3] 

    def get_data(self):
        # self.A = np.moveaxis(nib.orientations.apply_orientation(
        #     self.file.get_fdata(), 
        #     nib.orientations.axcodes2ornt('LPI')), # LPI = [R->L, A->P, S->I]
        #     0,1)
        return np.moveaxis(np.flip(self.file.get_fdata().astype(np.int32),[0,1,2]),0,1)

    def __init__(self, path):
        self.init_image(path)
    
class Segmentation:
    def get_data(self, tissue):
        self.A = np.moveaxis(
            nib.orientations.apply_orientation(self.files[tissue].get_fdata(), 
                                               nib.orientations.axcodes2ornt('LPI')), 0,1)
        self.affine = self.files[tissue].affine
        # self.A = np.moveaxis(np.flip(self.files[tissue].get_fdata(),[0,1,2]),0,1)

    def __init__(self, path):
        pattern = re.compile('[A-Za-z](\d+)')  
        files = [f for f in os.listdir(path) if 'vertebrae' in f]
        tissue = [f.split('.')[0] for f in files]
        # sorted(files,key = lambda x: (x[0], int(x[1:])))
        # tissue = [pattern.search(f).group() for f in files]
        self.files = dict(zip(tissue, [nib.load(join(path,f)) for f in files]))
        # self.tissue = sorted(tissue,key = lambda x: (x[0], int(x[1:])))
        self.tissue = tissue
        self.pixel_spacing = self.files[tissue[0]].header['pixdim'][1]
        self.slice_thickness = self.files[tissue[0]].header['pixdim'][3]

class load_map:
    _VERTEBRAE_KEY = 2 ** 0
    _PROCESS_KEY = 2 ** 1
    # Coronal Keys
    _ANTERIOR_KEY = 2 ** 2
    _POSTERIOR_KEY = 2 ** 3
    _CORONAL_KEYS = [_ANTERIOR_KEY, _POSTERIOR_KEY, _PROCESS_KEY]
    
    # Sagittal Keys
    _LEFT_KEY = 2 ** 4
    _RIGHT_KEY = 2 ** 5
    _SAGITTAL_KEYS = [_RIGHT_KEY, _LEFT_KEY]
    
    # Axial Keys
    _SUPERIOR_KEY = 2 ** 6
    _INFERIOR_KEY = 2 ** 7
    _AXIAL_KEYS = [_SUPERIOR_KEY, _INFERIOR_KEY]

def resample(ref_img_path:str , input_img_path:[str,list]):
    ref_img = sitk.ReadImage(ref_img_path, imageIO='NiftiImageIO')
    
    if isinstance(input_img_path, list):
        series_reader = sitk.ImageSeriesReader()
        series_reader.SetFileNames(input_img_path)
        series_reader.LoadPrivateTagsOn()
        input_img = series_reader.Execute()
    else:
        input_img = sitk.ReadImage(input_img_path)
    
    return sitk.GetArrayFromImage(sitk.Resample(input_img, ref_img))

def dir_walker(dir_path, subject):
    files = []
    for x,_,z in os.walk(join(dir_path, subject)):
        if (len(z) > 0) & ('DICOMDIR' not in z):
            sub_files = [join(x,f) for f in z]
            files += sub_files
    return files

def get_bounding_box(mask, pad=0):
    locs = np.nonzero(mask)
    
    if len(locs[0]) == 0:
        return None
    else:
        ex_pts = [(np.min(locs[0]),np.max(locs[0])), 
                  (np.min(locs[1]),np.max(locs[1])),
                  (np.min(locs[2]),np.max(locs[2])),]
        
        bounding_box = np.s_[ex_pts[0][0]-pad:ex_pts[0][1]+pad, 
                             ex_pts[1][0]-pad:ex_pts[1][1]+pad,
                             ex_pts[2][0]-pad:ex_pts[2][1]+pad]
    return bounding_box

def remove_bounding_box(bounded_image, bounding_box, orig_dims):
    out_mask = np.zeros(orig_dims)
    out_mask[bounding_box[0].start:bounding_box[0].stop,
             bounding_box[1].start:bounding_box[1].stop,
             bounding_box[2].start:bounding_box[2].stop] = bounded_image
    return out_mask

legend = {'vertebrae_S1': 26,
          'vertebrae_L5': 27,
          'vertebrae_L4': 28,
          'vertebrae_L3': 29,
          'vertebrae_L2': 30,
          'vertebrae_L1': 31,
          'vertebrae_T12':32,
          'vertebrae_T11':33,
          'vertebrae_T10':34,
          'vertebrae_T9': 35,
          'vertebrae_T8': 36,
          'vertebrae_T7': 37,
          'vertebrae_T6': 38,
          'vertebrae_T5': 39,
          'vertebrae_T4': 40,
          'vertebrae_T3': 41,
          'vertebrae_T2': 42,
          'vertebrae_T1': 43,
          'vertebrae_C7': 44,
          'vertebrae_C6': 45,
          'vertebrae_C5': 46,
          'vertebrae_C4': 47,
          'vertebrae_C3': 48,
          'vertebrae_C2': 49,
          'vertebrae_C1': 50}


def plot_all_planes(img, st, ps, plane_indeces, cm='nipy_spectral', axes=0):
    fig, axs = plt.subplots(1,3)
    axs[0].imshow(img[plane_indeces[0],...].T, cmap=cm)
    axs[0].plot(plane_indeces[1], plane_indeces[2], '+b', markersize=15)
    axs[0].set_title('Coronal')
    axs[0].set_xlabel('Sagittal')
    axs[0].set_ylabel('Axial')
    
    
    axs[1].imshow(img[:,plane_indeces[1],:].T, cmap=cm)
    axs[1].plot(plane_indeces[0], plane_indeces[2], '+b', markersize=15)
    axs[1].set_title('Sagittal')
    axs[1].set_xlabel('Coronal')
    axs[1].set_ylabel('Axial')
    
    
    axs[2].imshow(img[...,plane_indeces[2]], cmap=cm)
    axs[2].plot(plane_indeces[1], plane_indeces[0], '+b', markersize=15)
    axs[2].set_title('Axial')
    axs[2].set_xlabel('Sagittal')
    axs[2].set_ylabel('Coronal')
    
    
    axs[0].set_aspect(st/ps)
    axs[1].set_aspect(st/ps)
    
    if axes == 0:
        for ax in axs.ravel():
            ax.axis('off')
    # plt.show()
    return axs

# Depcrecated
def collect_individual_segmentation_files(mask, path, file_name='collected_spine.nii.gz'):
    collected_mask = np.zeros(mask.files[mask.tissue[0]].shape)
    affine_test = mask.files[mask.tissue[0]].affine
    
    for tissue in mask.tissue:
        if not np.allclose(affine_test, mask.files[tissue].affine):
            break
        else:
            print(tissue, mask.legend[tissue])
            collected_mask[mask.files[tissue].get_fdata() == 1] = mask.legend[tissue]
    
    img = nib.Nifti1Image(collected_mask, affine_test)
    nib.save(img, join(path, file_name))