# neuroCombat

## Overview
This repository contains the ComBat algorithm for correcting batch effects in neuroimaging (or microarray) data. This code runs ~150 times faster than the R version, and is incredibly simplified in that you do NOT have to create any design matrices, etc. All you have to do is pass in TWO numpy arrays or pandas DataFrames (the dataset to correct, and the dataset containing the batch/confound/target variables).

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

	pheno = pd.read_table('examples/sva/bladder-pheno.txt', index_col=0) # Y (cognitive) data
	dat = pd.read_table('examples/sva/bladder-expr.txt', index_col=0) # X (imaging) data)
	dat = dat.T

	categorical_targets = ["cancer"]
	numerical_targets 	= ["age"]
	batch_var 			= "batch"

	result = neuroCombat(X=dat, Y=pheno, batch_var=batch_var,
		categorical_targets=categorical_targets, numerical_targets=numerical_targets)
```



### Correcting from Numpy Arrays
NOTE: If you read in the Y dataset as a numpy array, you MUST include `y_feature_labels`, which is a numpy array (or python list) of strings containing the Y column labels.

```python
	import numpy as np
	from neuroCombat.neuroCombat import neuroCombat

	X = np.load('examples/sva/bladder-expr.npy')
	X = X.T # neuroimaging people like our data to be shape = (samples, features)
	Y = np.load('examples/sva/bladder-pheno.npy') # shape = (samples, features)
	y_feature_labels = np.load('examples/sva/feature_labels.npy')

	categorical_targets = ["cancer"]
	numerical_targets 	= ["age"]
	batch_var 			= "batch"

	result1 = neuroCombat(X=X, Y=Y, batch_var=batch_var,
		categorical_targets=categorical_targets, numerical_targets=numerical_targets,
		y_feature_labels=y_feature_labels)
```


## Performance
On the example above, we report the following:

- Combat.py repository (not mine): ~0.6s
- R's SVA package : ~10s
- neuroCombat repository (mine): ~0.1 s

So, Combat.py repository reports ~ 24 times speed on the R version, and we report ~ 6 times speed speedup on the Combat.py version. This means ~ 150 time speedup on the R version. This margin gets larger as the dataset gets larger.