"""
Functions for processing neuroimaging files in order
to be fed into neuroCombat function
"""

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

masked_imgs = masker.transform(img_paths)
