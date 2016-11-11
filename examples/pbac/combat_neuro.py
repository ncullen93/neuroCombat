
from neuroCombat.neuroCombat import neuroCombat

X_dir 	= 'images/'
mask 	= 'mask/mask.nii.gz'
# make random batch effect
Y 		= np.zeros((80,1))
Y[:20] 	= 1
Y[20:40]= 2
# y labels .. no confounds/targets
y_labels = np.array(['batch'])
batch_var = 'batch'

neuroCombat(X=X_dir, mask=mask,
	Y=Y, y_feature_labels=y_labels,
	batch_var='batch', save_dir='combat_images/')











