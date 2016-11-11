# neuroCombat

## Overview
This repository contains the ComBat algorithm for correcting batch effects in neuroimaging (or microarray) data. This code runs ~150 times faster than the R version, and is incredibly simplified in that you do NOT have to create any design matrices, etc. All you have to do is pass in TWO numpy arrays or pandas DataFrames (the dataset to correct, and the dataset containing the batch/confound/target variables).

## Example
```python
	X = np.load('bladder-expr.npy')
	Y = np.load('combat_py/bladder-pheno.npy')
	y_feature_labels = np.load('combat_py/feature_labels.npy')
	categorical_targets = ["cancer"]
	numerical_targets 	= ["age"]
	batch_var 			= "batch"

	result = neuroCombat(X=X, Y=Y, batch_var=batch_var,
		categorical_targets=categorical_targets, numerical_targets=numerical_targets,
		y_feature_labels=y_feature_labels)
```

## Performance
On the example above, we report the following:

- Combat.py repository (not mine): ~0.6s
- R's SVA package : ~10s
- neuroCombat repository (mine): ~0.1 s

So, Combat.py repository reports ~ 24 times speed on the R version, and we report ~ 6 times speed speedup on the Combat.py version. This means ~ 150 time speedup on the R version. This margin gets larger as the dataset gets larger.