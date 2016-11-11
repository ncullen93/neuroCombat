"""
Functions for processing neuroimaging files in order
to be fed into neuroCombat function


from nilearn.input_data import NiftiMasker
import os

mask_img = '/users/nick/desktop/pbac_images/mask/mask.nii.gz'

masker = NiftiMasker()
masker.fit(mask_img)

img_dir = '/users/nick/desktop/pbac_images/images/'

img_paths = os.listdir(img_dir)
try:
	img_paths.remove('.DS_Store')
except:
	pass
img_paths = [img_dir+p for p in img_paths]

masked_data = masker.transform(img_paths)

### run combat on masked_data ###
#combat_data = neuroCombat(masked_data, ...)

### inverse transform to get 3D images back from numpy array
#corrected_imgs = masker.inverse_transform(combat_data)

## save the images
"""

import os

import nibabel as nib
import nilearn
from nilearn.input_data import NiftiMasker

def load_imgs_from_dir(load_dir, mask):
	if load_dir[-1] != '/': load_dir += '/'
	img_names 	= os.listdir(load_dir)
	try:
		img_names.remove('.DS_Store')
	except:
		pass
	
	img_paths = [load_dir+name for name in img_names]

	assert (mask is not None), 'Must include mask'

	masker = NiftiMasker()
	masker.fit(mask)
	return_data = masker.transform(img_paths)
	print 'Found %i images w/ %i unmasked voxels' % \
				(return_data.shape[0], return_data.shape[1])
	
	return return_data

def save_imgs_to_dir(data, mask, save_dir, load_dir):
	if not os.path.isdir(save_dir):
		os.mkdir(save_dir)
	
	if load_dir[-1] != '/': load_dir += '/'
	if save_dir[-1] != '/': save_dir += '/'

	orig_names 	= os.listdir(load_dir)
	try:
		orig_names.remove('.DS_Store')
	except:
		pass

	masker = NiftiMasker()
	masker.fit(mask)
	save_imgs = []
	for subject in range(data.shape[0]):
		save_imgs.append(masker.inverse_transform(data[subject,:]))

	for i, orig_name in enumerate(orig_names):
		save_path = save_dir+'corrected_'+orig_name
		img = save_imgs[i]
		nib.save(img, save_path)











