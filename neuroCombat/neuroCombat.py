"""
ComBat for correcting batch effects in neuroimaging data
"""
#from __future__ import division
import pandas as pd
import numpy as np
import numpy.linalg as la


def neuroCombat(X, Y, batch_var, categorical_targets=None, numerical_targets=None,
	y_feature_labels=None):
	"""
	Run ComBat to correct for batch effects in neuroimaging data

	Arguments
	---------
	X : a pandas data frame or numpy array
		- neuroimaging data
		- shape = (features, samples)

	Y : a pandas data frame (or numpy array if included y_feature_labels)
		- demographic/phenotypic/behavioral/batch data 
		- shape = (samples, features)
	
	batch_vars : string
		- batch effect variable
		- e.g. scan site

	cat_target_vars : string or list of strings
		- variables which are categorical that you want to predict
		- e.g. binary depression or no depression

	num_target_vars : string or list of strings
		- variables which are numerical that you want to predict
		- e.g. depression sub-scores

	y_feature_labels : a list of strings
		- column labels for Y dataset (must includes batch variable, target variables, etc)
	"""
	##############################
	### CLEANING UP INPUT DATA ###
	##############################
	if not isinstance(categorical_targets, list):
		if categorical_targets is None:
			categorical_targets = []
		else:
			categorical_targets = [categorical_targets]
	if not isinstance(numerical_targets, list):
		if numerical_targets is None:
			numerical_targets = []
		else:
			numerical_targets = [numerical_targets]

	if isinstance(Y, np.ndarray):
		assert (y_feature_labels is not None), 'Must include y_feature_labels w/ numpy array input'
		Y = Y.astype('object')
		for i in range(Y.shape[-1]):
			try:
				Y[:,i] = Y[:,i].astype('float32')
			except:
				pass
	elif isinstance(Y, pd.DataFrame):
		y_feature_labels = np.array(Y.columns)
		Y = np.array(Y, dtype='object') 
		for i in range(Y.shape[-1]):
			try:
				Y[:,i] = Y[:,i].astype('float32')
			except:
				pass
	if isinstance(X, pd.DataFrame):
		X = np.array(X,dtype='float32')

	##############################

	# get column indices for relevant variables
	batch_col 	= np.where(y_feature_labels==batch_var)[0][0]
	cat_cols 	= [np.where(y_feature_labels==c_var)[0][0] for c_var in categorical_targets]
	num_cols 	= [np.where(y_feature_labels==n_var)[0][0] for n_var in numerical_targets]

	# create dictionary that stores batch info
	(batch_levels, sample_per_batch)= np.unique(Y[:,batch_col],return_counts=True)
	info_dict = {}
	info_dict['batch_levels']		= batch_levels.astype('int')
	info_dict['n_batch'] 			= len(batch_levels)
	info_dict['n_sample'] 			= int(Y.shape[0])
	info_dict['sample_per_batch']	= sample_per_batch.astype('int')
	info_dict['batch_info']			= [list(np.where(Y[:,batch_col]==idx)[0]) for idx in batch_levels]

	# create design matrix
	design 					= make_design_matrix(Y, batch_col, cat_cols, num_cols)
	
	# standardize data across features
	s_data, s_mean, v_pool 	= standardize_across_features(X, design, info_dict)
	
	# fit L/S models and find priors
	LS_dict 				= fit_LS_model_and_find_priors(s_data, design, info_dict)

	# find parametric adjustments
	gamma_star, delta_star 	= find_parametric_adjustments(s_data, LS_dict, info_dict)

	# adjust data
	bayes_data 				= adjust_data_final(s_data, design, gamma_star, delta_star, 
												s_mean, v_pool, info_dict)

	return np.array(bayes_data)

def make_design_matrix(Y, batch_col, cat_cols, num_cols):
	"""
	Return Matrix containing the following parts:
		- one-hot matrix of batch variable (full)
		- one-hot matrix for each categorical_targts (removing the first column)
		- column for each numerical_targets
	"""
	def to_categorical(y, nb_classes=None):
		if not nb_classes:
			nb_classes = np.max(y)+1
		Y = np.zeros((len(y), nb_classes))
		for i in range(len(y)):
			Y[i, y[i]] = 1.
		return Y
	
	hstack_list = []

	### batch one-hot ###
	batch = np.array(Y[:,batch_col], dtype='int') # batch_vars
	batch = batch - (np.min(batch) - 0) # min = zero
	batch_onehot = to_categorical(batch, len(np.unique(batch)))
	hstack_list.append(batch_onehot)

	### categorical one-hots ###
	for cat_col in cat_cols:
		cat = np.unique(np.array(Y[:,cat_col]),return_inverse=True)[1]
		cat_onehot = to_categorical(cat, len(np.unique(cat)))[:,1:]
		hstack_list.append(cat_onehot)

	### numerical vectors ###
	for num_col in num_cols:
		num = np.array(Y[:,num_col],dtype='float32')
		num = num.reshape(num.shape[0],1)
		hstack_list.append(num)

	design = np.hstack(hstack_list)
	return design

def standardize_across_features(X, design, info_dict):
	n_batch 			= info_dict['n_batch']
	n_sample 			= info_dict['n_sample']
	sample_per_batch	= info_dict['sample_per_batch']

	B_hat = np.dot(np.dot(la.inv(np.dot(design.T, design)), design.T), X.T)
	grand_mean = np.dot((sample_per_batch/ float(n_sample)).T, B_hat[:n_batch,:])
	var_pooled = np.dot(((X - np.dot(design, B_hat).T)**2), np.ones((n_sample, 1)) / float(n_sample))

	stand_mean = np.dot(grand_mean.T.reshape((len(grand_mean), 1)), np.ones((1, n_sample)))
	tmp = np.array(design.copy())
	tmp[:,:n_batch] = 0
	stand_mean  += np.dot(tmp, B_hat).T

	s_data = ((X- stand_mean) / np.dot(np.sqrt(var_pooled), np.ones((1, n_sample))))

	return s_data, stand_mean, var_pooled

def aprior(gamma_hat):
	m = np.mean(gamma_hat)
	s2 = np.var(gamma_hat,ddof=1)
	return (2 * s2 +m**2) / float(s2)

def bprior(gamma_hat):
	m = gamma_hat.mean()
	s2 = np.var(gamma_hat,ddof=1)
	return (m*s2+m**3)/s2

def postmean(g_hat, g_bar, n, d_star, t2):
	return (t2*n*g_hat+d_star * g_bar) / (t2*n+d_star)

def postvar(sum2, n, a, b):
	return (0.5 * sum2 + b) / (n / 2.0 + a - 1.0)

def fit_LS_model_and_find_priors(s_data, design, info_dict):
	n_batch 	= info_dict['n_batch']
	batch_info 	= info_dict['batch_info'] 
	
	batch_design = design[:,:n_batch]
	gamma_hat = np.dot(np.dot(la.inv(np.dot(batch_design.T, batch_design)), batch_design.T), s_data.T)

	delta_hat = []
	for i, batch_idxs in enumerate(batch_info):
		delta_hat.append(np.var(s_data[:,batch_idxs],axis=1,ddof=1))
	
	gamma_bar = np.mean(gamma_hat, axis=1) 
	t2 = np.var(gamma_hat,axis=1, ddof=1)

	a_prior = list(map(aprior, delta_hat))
	b_prior = list(map(bprior, delta_hat))

	LS_dict = {}
	LS_dict['gamma_hat'] 	= gamma_hat
	LS_dict['delta_hat'] 	= delta_hat
	LS_dict['gamma_bar']	= gamma_bar
	LS_dict['t2'] 			= t2
	LS_dict['a_prior'] 		= a_prior
	LS_dict['b_prior']		= b_prior
	return LS_dict

def it_sol(sdat, g_hat, d_hat, g_bar, t2, a, b, conv=0.0001):
	n = (1 - np.isnan(sdat)).sum(axis=1)
	g_old = g_hat.copy()
	d_old = d_hat.copy()

	change = 1
	count = 0
	while change > conv:
		g_new = postmean(g_hat, g_bar, n, d_old, t2)
		sum2 = ((sdat - np.dot(g_new.reshape((g_new.shape[0], 1)), np.ones((1, sdat.shape[1])))) ** 2).sum(axis=1)
		d_new = postvar(sum2, n, a, b)

		change = max((abs(g_new - g_old) / g_old).max(), (abs(d_new - d_old) / d_old).max())
		g_old = g_new #.copy()
		d_old = d_new #.copy()
		count = count + 1
	adjust = (g_new, d_new)
	return adjust 

def find_parametric_adjustments(s_data, LS, info_dict):
	batch_info 	= info_dict['batch_info'] 

	gamma_star, delta_star = [], []
	for i, batch_idxs in enumerate(batch_info):
		temp = it_sol(s_data[:,batch_idxs], LS['gamma_hat'][i],
					LS['delta_hat'][i], LS['gamma_bar'][i], LS['t2'][i], 
					LS['a_prior'][i], LS['b_prior'][i])

		gamma_star.append(temp[0])
		delta_star.append(temp[1])

	return np.array(gamma_star), np.array(delta_star)

def adjust_data_final(s_data, design, gamma_star, delta_star, stand_mean, var_pooled, info_dict):
	sample_per_batch 	= info_dict['sample_per_batch']
	n_batch 			= info_dict['n_batch']
	n_sample	 		= info_dict['n_sample']
	batch_info 			= info_dict['batch_info']

	batch_design = design[:,:n_batch]

	bayesdata = s_data
	gamma_star = np.array(gamma_star)
	delta_star = np.array(delta_star)

	for j, batch_idxs in enumerate(batch_info):
		dsq = np.sqrt(delta_star[j,:])
		dsq = dsq.reshape((len(dsq), 1))
		denom =  np.dot(dsq, np.ones((1, sample_per_batch[j])))
		numer = np.array(bayesdata[:,batch_idxs] - np.dot(batch_design[batch_idxs,:], gamma_star).T)

		bayesdata[:,batch_idxs] = numer / denom

	vpsq = np.sqrt(var_pooled).reshape((len(var_pooled), 1))
	bayesdata = bayesdata * np.dot(vpsq, np.ones((1, n_sample))) + stand_mean

	return bayesdata


if __name__=='__main__':
	X = np.load('../examples/sva/bladder-expr.npy')
	Y = np.load('../examples/sva/bladder-pheno.npy')
	y_feature_labels = np.load('../examples/sva/feature_labels.npy')
	categorical_targets = ["cancer"]
	numerical_targets 	= []
	batch_var 			= "batch"

	result = neuroCombat(X=X, Y=Y, batch_var=batch_var,
		categorical_targets=categorical_targets, numerical_targets=numerical_targets,
		y_feature_labels=y_feature_labels)