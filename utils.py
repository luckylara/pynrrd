import nibabel as nib
import numpy as np
from nrrd import *
import os

def nifti_to_nrrd(nifti_filename, nrrd_filename=None):
    if not nrrd_filename:
        basename = os.path.basename(nifti_filename).split('.')[0]
        nrrd_filename = basename+'.nrrd'

    nifti_img = nib.load(nifti_filename)
    data = nifti_img.get_data()

    nrrd_header = NrrdHeader.fromNiftiHeader(nifti_img.get_header())
    


    return nrrd_filename
