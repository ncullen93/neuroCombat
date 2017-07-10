# neuroCombat

## Overview
This repository contains the ComBat algorithm for correcting batch effects in neuroimaging (or microarray) data. This code runs ~10-50 times faster than the R version, and is incredibly simplified in that you do NOT have to create any design matrices, etc. All you have to do is pass in TWO numpy arrays or pandas DataFrames (the dataset to correct, and the dataset containing the batch/confound/target variables).

Combining the ease-of-use, the ability to read neuroimages directly, the significant computational speedup, and the Python language, neuroCombat should fulfill an important niche in the neuroscience community.

To note, if you need to read a certain type of neuroimage (e.g. fMRI) that is not supported yet, you can simply write a description of the data and submit an <b>Issue</b> and I'll implement it.


## Installation
1. Download zipped repository
2. Unpack
3. cd neuroCombat-master
4. run `python setup.py install` or `pip install .`
5. To use the `neuroCombat` function once installed:
    - `from neuroCombat import neuroCombat`

## neuroCombat Function

```
neuroCombat(data, covars, batch_col, discrete_cols=None, continuous_cols=None)

Docstring:
Run ComBat to correct for batch effects in neuroimaging data

Arguments
---------
data : a pandas data frame or numpy array
    neuroimaging data to correct with shape = (samples, features)
    e.g. cortical thickness measurements, image voxels, etc

covars : a pandas data frame w/ shape = (samples, features)
    demographic/phenotypic/behavioral/batch data 
    
batch_col : string
    - batch effect variable
    - e.g. scan site

discrete_cols : string or list of strings
    - variables which are categorical that you want to predict
    - e.g. binary depression or no depression

continuous_cols : string or list of strings
    - variables which are continous that you want to predict
    - e.g. depression sub-scores

Returns
-------
- A numpy array with the same shape as `data` which has now been ComBat-corrected
```

## Correspondance to SVA (R) version of ComBat

<b> The most significant difference between the SVA and neuroCombat algorithms is that the neuroCombat version accepts input data X of shape (n_samples, n_features) while SVA accepts (n_features, n_samples).</b>

In SVA's version of ComBat, you might do the following to correct for dataset `X` while adjusting for covariates `Y$c1` and `Y$c2`:

```R
batch <- Y$batch
model <- model.matrix(~ c1 + c2, data=Y)
combat_data <- ComBat(dat=X,batch=batch, mod=model)
```

In this Python version, you would do the following (assuming `covars` is a Pandas Dataframe with appropriate column labels):

```Python
combat_data = neuroCombat(data=X, 
                          covars=Y, 
                          batch_col='batch', 
                          discrete_cols=['c1','c2'])
```

As you see, there is no need for a model matrix. 


## Examples

No matter what, the `covars` argument must be a pandas.DataFrame. However, for the `data` argument, you may
either give a numpy array or a pandas dataframe. Here is the difference between the two (note the result
will be the same):

### Correcting from Numpy Array as Data
```python
from neuroCombat import neuroCombat

import pandas as pd
import numpy as np

data = np.load('data/bladder-expr.npy')
covars = pd.read_csv('data/bladder-pheno.txt', delimiter='\t')

discrete_cols = ['cancer']
continuous_cols = ['age']
batch_col = 'batch'

data_combat = neuroCombat(data=data,
                          covars=covars,
                          batch_col=batch_col,
                          discrete_cols=discrete_cols,
                          continuous_cols=continuous_cols)
```

### Correcting from pandas.DataFrame as Data
```python
from neuroCombat import neuroCombat

import pandas as pd
import numpy as np

data = pd.read_csv('data/bladder-expr.txt', delimiter='\t')
covars = pd.read_csv('data/bladder-pheno.txt', delimiter='\t')

discrete_cols = ['cancer']
continuous_cols = ['age']
batch_col = 'batch'

data_combat = neuroCombat(data=data,
                          covars=covars,
                          batch_col=batch_col,
                          discrete_cols=discrete_cols,
                          continuous_cols=continuous_cols)
```

## Performance
On the bladder example above, we report the following:

- Combat.py repository (not mine): ~0.6s
- R's SVA package : ~10s
- neuroCombat repository (mine): ~0.1 s