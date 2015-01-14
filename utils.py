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

def nrrd_to_nifti(nrrd_filename, nifti_filename=None):
    if not nifti_filename:
        basename = os.path.basename(nrrd_filename).split('.')[0]
        nrrd_filename = basename+'.nii.gz'
    else:
        basename = nifti_filename.split('.')[0]

    hdr, data = NrrdReader().load(nrrd_filename)
    hdr.correctSpaceRas()

    newimg = nib.Nifti1Image(data, hdr.getAffine())
    nib.save(newimg, nifti_filename)

    if hdr.isDTMR():
        bvecs = np.array(hdr.getDwiGradients())
        np.savetxt(basename+'.bvec', bvecs, fmt='%.6f')

        bvals = np.array(hdr.getBvals())
        np.savetxt(basename+'.bval', bvals, fmt='%.2f')
        
    return nrrd_filename    