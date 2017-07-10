# neuroCombat

## Overview
This repository contains the ComBat algorithm for correcting batch effects in neuroimaging (or microarray) data. This code runs ~10-50 times faster than the R version, and is incredibly simplified in that you do NOT have to create any design matrices, etc. All you have to do is pass in TWO numpy arrays or pandas DataFrames (the dataset to correct, and the dataset containing the batch/confound/target variables).

Combining the ease-of-use, the ability to read neuroimages directly, the significant computational speedup, and the Python language, neuroCombat should fulfill an important niche in the neuroscience community.

To note, if you need to read a certain type of neuroimage (e.g. fMRI) that is not supported yet, you can simply write a description of the data and submit an <b>Issue</b> and I'll implement it.


## Correspondance to SVA (R) version of ComBat

<b> The most significant difference between the SVA and neuroCombat algorithms is that the neuroCombat version accepts input data X of shape (n_samples, n_features) while SVA accepts (n_features, n_samples).</b>

In SVA's version of ComBat, you might do the following to correct for dataset `X` while adjusting for covariates `Y$c1` and `Y$c2`:

```R
batch <- Y$batch
model <- model.matrix(~ c1 + c2, data=Y)
combat_data <- ComBat(dat=X,batch=batch, mod=model)
```

In this Python version, you would do the following (assuming Y is a Pandas Dataframe with appropriate column labels):

```Python
categorical_targets = ['c1','c2']
batch_var = 'batch'
combat_data = neuroCombat(X=X, Y=Y, batch_var=batch_var, categorical_targets=categorical_targets)
```

As you see, there is no need for a model matrix. As we will see below, there is also the possibility of reading in neuroimages (e.g. in raw Nifti format) from a file directory. 

## Installation
1. Download zipped repository
2. Unpack
3. cd neuroCombat-master
4. run `python setup.py install`
5. To use neuroCombat:
    - `from neuroCombat.neuroCombat import neuroCombat`
    - Now you can use the function

## Examples

### Correcting from Pandas Dataframes
```python
import pandas as pd
from neuroCombat.neuroCombat import neuroCombat

pheno = pd.read_table('examples/bladder/bladder-pheno.txt', index_col=0) # Y (cognitive) data
dat = pd.read_table('examples/bladder/bladder-expr.txt', index_col=0) # X (imaging) data)
dat = dat.T

categorical_targets = ["cancer"]
numerical_targets   = ["age"]
batch_var           = "batch"

result = neuroCombat(X=dat, Y=pheno, batch_var=batch_var,
    categorical_targets=categorical_targets, numerical_targets=numerical_targets)
```



### Correcting from Numpy Arrays
NOTE: If you read in the Y dataset as a numpy array, you MUST include `y_feature_labels`, which is a numpy array (or python list) of strings containing the Y column labels.

```python
import numpy as np
from neuroCombat.neuroCombat import neuroCombat

X = np.load('examples/bladder/bladder-expr.npy')
X = X.T # neuroimaging people like our data to be shape = (samples, features)
Y = np.load('examples/bladder/bladder-pheno.npy') # shape = (samples, features)
y_feature_labels = np.load('examples/sva/feature_labels.npy')

categorical_targets = ["cancer"]
numerical_targets   = ["age"]
batch_var           = "batch"

result1 = neuroCombat(X=X, Y=Y, batch_var=batch_var,
    categorical_targets=categorical_targets, numerical_targets=numerical_targets,
    y_feature_labels=y_feature_labels)
```

### Correcting from directory of Structural Nifti images with One Mask
The beauty of neuroCombat is that you can correct neuroimages without even loading them yourself. All that work is done for you, and all you have to do is give the directory to the set of images.

Take this example, with a set of T1 images found in the 'pbac/images' directory. We give the path to those images, the path to a Mask file, and the batch data. That's it. The corrected images will be save to the directory of your choice and will be of the same name as the originals but prefixed with 'corrected_'.

```python
from neuroCombat.neuroCombat import neuroCombat

X_dir   = 'examples/T1/images/'
mask    = 'examples/T1/mask/mask.nii.gz'
# make random batch effect
Y       = np.zeros((80,1))
Y[:20]  = 1
Y[20:40]= 2
# y labels .. no confounds/targets
y_labels = np.array(['batch'])
batch_var = 'batch'

neuroCombat(X=X_dir, mask=mask,
    Y=Y, y_feature_labels=y_labels,
    batch_var='batch', save_dir='combat_images/')
```

### Correcting from directory of Structural Nifti images with Individual Masks
If every subject has their own masks, that's completely fine. Instead of passing in the path to a single Nifti mask file, simply pass in the mask directory.

In this case, the images for each subject have their own directory which includes the masks. Note that when you pass in a value for the `save_dir` argument, then the ComBat-corrected images get re-masked and saved to that save directory in Nifti format. Otherwise, we return the corrected images (NOT re-masked) as an array.

e.g:
    images/
        subject1/
            t1.nii.gz
            mask.nii.gz
        subject2/
            t1.nii.gz
            mask.nii.gz
        ..
```python
img_dir     = 'T1/'
mask_dir    = 'T1/'

# read in demographic data as usual
Y = pd.read_csv('subject_info/demographics.csv')
cat_covars  = ['gender','group']
num_covars  = ['age']
batch   = 'scanner'
save_dir= T1_combat'

neuroCombat(X=img_dir, mask=mask_dir, Y=Y, batch=batch, cat_covars=cat_covars,
    num_covars=num_covars, save_dir=save_dir)
```

### Loading Images from directory as Quality Assurance
If you want to make sure that your images are being loaded correctly, you can import the functionality which neuroCombat uses and load the images from directory yourself. Then, you can just pass in that data for the 'X' argument.

```python
from neuroCombat.neuroCombat import neuroCombat
from neuroCombat.neuroimage_process import load_imgs_from_dir

img_dir = 'T1/'
mask_dir = 'T1/'
# load the neuroimaging data from directory yourself
image_data = load_imgs_from_dir(img_dir, mask_dir)

# read in demographic data as usual
Y = pd.read_csv('subject_info/demographics.csv')
cat_covars  = ['gender','group']
num_covars  = ['age']
batch   = 'scanner'
save_dir= T1_combat'

combat_data = neuroCombat(X=image_data, Y=Y, batch=batch, cat_covars=cat_covars,
    num_covars=num_covars)
```

## Performance
On the bladder example above, we report the following:

- Combat.py repository (not mine): ~0.6s
- R's SVA package : ~10s
- neuroCombat repository (mine): ~0.1 s