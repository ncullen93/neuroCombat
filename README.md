# neuroCombat

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