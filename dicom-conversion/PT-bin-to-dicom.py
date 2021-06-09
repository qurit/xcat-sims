# -*- coding: utf-8 -*-
"""
Summary:

- File converts binary files to dicom and applies header information
        Extracts binary from folder called "runs"
        Creates dicom file in "dicom" folder
        
- File is intended to be used after ground truth image has been generated from
        4D-extended cardiac torso (XCAT) phantom and PET image is subsequently simulated
        
- This code accepts a binary file with 128x128x301 matrix and 3x 3x 1.5mm voxel sizes
        File must be modified for different file structures

- File utilizes sample code from PyDICOM:
        https://pydicom.github.io/pydicom/stable/auto_examples/index.html
        

Reference:
Please use the following reference if you publish results obtained using this software:

    R. Fedrigo, et al., 
    “Development of the Lymphatic System in the 4D XCAT Phantom for Improved Multimodality Imaging Research”, 
    J. Nuc. Med., vol. 62, publication 113, 2021.

Contact
For support, contributions, and questions, please contact rfedrigo (at) student (dot) ubc (dot) ca
"""

from __future__ import print_function
import SimpleITK as sitk
import time, os
import numpy as np
import shutil

################################################
# INPUT FILE/DIRECTORY NAMES
################################################
myfilename = "snm1_act_1_R1.bin"

mydirectory=os.getcwd()
mydirectory=mydirectory+"/"
mynewfolder = "snm1_R1"

desc = "SNM1_R1" # series description
name = "SNM1" # patient name/ID

################################################
# OPEN BINARY FILE AND CONVERT TO NUMPY ARRAY
################################################
full_file = "input/"+myfilename
out_dir = mydirectory + "../output/" + mynewfolder

# Create a new series from a numpy array
new_arr = np.fromfile(full_file, dtype='single').astype(np.int16)
new_arr = np.reshape(new_arr,(301,128,128)).astype(np.int16)
    
new_img = sitk.GetImageFromArray(new_arr)
new_img.SetSpacing([3.0,3.0, 1.5])

################################################
# CREATE DICOM FILE AND APPLY HEADER INFO
################################################
z = mydirectory + "../output/" + mynewfolder
shutil.rmtree(out_dir)
os.mkdir(z)
    
# Write the 3D image as a series
writer = sitk.ImageFileWriter()
writer.KeepOriginalImageUIDOn()
modification_time = time.strftime("%H%M%S")
modification_date = time.strftime("%Y%m%d")
    
direction = new_img.GetDirection()
series_tag_values = [("0008|0031",modification_time), # Series Time
                      ("0008|0021",modification_date), # Series Date
                      ("0008|0008","DERIVED\\SECONDARY"), # Image Type
                      ("0020|000e", "1.2.826.0.1.3680043.2.1125."+modification_date+".1"+modification_time), # Series Instance UID
                      ("0020|0037", '\\'.join(map(str, (direction[0], direction[3], direction[6],# Image Orientation (Patient)
                                                        direction[1],direction[4],direction[7])))),
                      ("0008|103e", desc), # Series Description
                      ("0054|1001", "BQML"), # radioactivity units
                      
                      #patient information
                      ("0010|0010", name), # Patient name
                      ("0010|0020", name), # Patient ID
                      
                      # image space information
                      ("0028|0010", "128"), # rows
                      ("0028|0011", "128"), # columns
                      ("0018|0050", "1.5"), # slice thickness
                      ("0054|0081", "301")  # number of slices
                      ]
    
for i in range(new_img.GetDepth()):
    image_slice = new_img[:,:,i]
    # Tags shared by the series.
    for tag, value in series_tag_values:
        image_slice.SetMetaData(tag, value)
            
    # Slice specific tags.
    image_slice.SetMetaData("0008|0012", time.strftime("%Y%m%d")) # Instance Creation Date
    image_slice.SetMetaData("0008|0013", time.strftime("%H%M%S")) # Instance Creation Time
        
    # Setting the type to CT preserves the slice location.
    image_slice.SetMetaData("0008|0060", "PT") # set the type as a PET image
        
    # (0020, 0032) image position patient determines the 3D spacing between slices.
    image_slice.SetMetaData("0020|0032", '\\'.join(map(str,new_img.TransformIndexToPhysicalPoint((0,0,i))))) # Image Position (Patient)
    image_slice.SetMetaData("0020,0013", str(i)) # Instance Number
    
    # Write to the output directory and add the extension dcm, to force writing in DICOM format.
    writer.SetFileName(os.path.join(out_dir,str(i)+'.dcm'))
    writer.Execute(image_slice)