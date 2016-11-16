"""
Functions for processing neuroimaging files in order
to be fed into neuroCombat function
"""

import os

import nibabel as nib
import nilearn
from nilearn.input_data import NiftiMasker
import nilearn.masking as nim


def load_imgs_from_dir(load_dir, mask):
	load_files = [l for l in os.listdir(load_dir) if l[0] !='.']

	# check if images are all in one directory
	imgs_in_top_dir = np.sum(['.nii.gz' in l for l in load_files]) > (0.8*len(load_files))
	if imgs_in_top_dir:
		data_names 	= load_files
		mask_names  = mask
	else:
		# if imgs in their own subject dir
		mask_names = []
		data_names = []
		for subject_dir in load_files:
			subject_imgs = os.listdir(os.path.join(load_dir,subject_dir))
			mask_img = [img for img in subject_imgs if 'mask' in img][0]
			data_img = [img for img in subject_imgs if 'nii.gz' in img and 'mask' not in img][0]
			mask_names.append(os.path.join(subject_dir,mask_img))
			data_names.append(os.path.join(subject_dir,data_img))

	data_paths = [os.path.join(load_dir,name) for name in data_names]
	if isinstance(mask_names, list):
		mask_paths = [os.path.join(load_dir,name) for name in mask_names]
	else:
		mask_paths = os.path.join(load_dir, mask_names)

	return_data = []
	if not isinstance(mask_paths, list):
		masker = NiftiMasker()
		masker.fit(mask_paths)
		return_data = masker.transform(data_paths)
	else:
		count = 0

		assert len(data_paths) == len(mask_paths), 'Found unequal number of masks and images'
		for data_path,mask_path in zip(data_paths,mask_paths):
			print 'loading and masking img %i/%i' % (count,len(data_paths))
			masker = NiftiMasker()
			masker.fit(mask_path)
			data = masker.transform(data_path)
			return_data.append(data.flatten())
			count += 1

		# pad zeros to make up for shorter masks
		max_len = np.max([len(i) for i in return_data])
		temp_data = np.zeros((len(return_data),max_len))
		for i,img in enumerate(return_data):
			temp_data[i,:len(img)] = img
		return_data = temp_data

	print 'Found %i images' % return_data.shape[0]
	
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











