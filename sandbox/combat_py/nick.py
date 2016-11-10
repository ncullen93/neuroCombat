"""
How it works:
 	- need X data (imaging)
 	- need Y data (cognitive)
 	- pick out adjustment (confound) variables from Y (e.g age, sex)
 	- pick out variables of interest (target vars) from Y (e.g. depression score)
 	- pick out batch variables from Y (e.g. scan site)

Null Model = adjustment variables but not target variables and NOT batch variables
Full Model = adjustment vars AND target vars
"""


import pandas as pd
import patsy
import sys
import numpy.linalg as la
import numpy as np


# (57, 5)
pheno = pd.read_table('bladder-pheno.txt', index_col=0)
# (22283, 57)
data = pd.read_table('bladder-expr.txt', index_col=0)

# all these are from phenotype (i.e. demographic) data
#confound_vars 	= ["age"]
#batch_vars		= ["batch"]
target_vars 	= ["cancer","age"]

#ebat = combat(data, pheno, confound_vars, batch_vars, target_vars)


batch 		= pheno.batch
model_str	= '~ ' + ' + '.join(target_vars)
model 		= patsy.dmatrix(model_str, pheno, return_type="dataframe")

### INPUT TO COMBAT : data (genes/imaging), batch (vector), model (pheno/demog)
data  = data
batch = batch
model = model
numerical_covariates = ["age"]

#### COMBAT CODE ####

# add batch var to model
if model is not None and isinstance(model, pd.DataFrame):
    model["batch"] = list(batch)
else:
    model = pd.DataFrame({'batch': batch})

# get batch information
batch_items = model.groupby("batch").groups.items() # list of (batch_idx, list of samples in batch)
batch_levels = [k for k, v in batch_items] # list of batch idxs
batch_info = [v for k, v in batch_items] # list of list of samples in each batch by idx
n_batch = len(batch_info) # number batches
n_batches = np.array([len(v) for v in batch_info]) # number of samples per batch in list
n_array = float(sum(n_batches)) # total number of samples

# remove intercept column from model
drop_cols = [cname for cname, inter in  ((model == 1).all()).iterkv() if inter == True]
drop_idxs = [list(model.columns).index(cdrop) for cdrop in drop_cols]
model = model[[c for c in model.columns if not c in drop_cols]]

# get model column idxs of numerical covariates
numerical_covariates = [list(model.columns).index(c) if isinstance(c, str) else c
        for c in numerical_covariates if not c in drop_cols]


##### make new design matrix = one-hot batch matrix + one-hot target categorical cols + numerical cols
#design = design_mat(model, numerical_covariates, batch_levels)
mod = model.copy()
numerical_covariates = numerical_covariates
batch_levels = batch_levels

# shape = (samples, n_batches) .. essentially a ONE-HOT matrix
design = patsy.dmatrix("~ 0 + C(batch, levels=%s)" % str(batch_levels),
                                              mod, return_type="dataframe")
# drop batch from mod
mod = mod.drop(["batch"], axis=1)
# columns idxs that are numerical
numerical_covariates = list(numerical_covariates)
sys.stderr.write("found %i batches\n" % design.shape[1])
# columns names that are not numerical
other_cols = [c for i, c in enumerate(mod.columns)
              if not i in numerical_covariates]
# factor matrix - intercept + target variables
factor_matrix = mod[other_cols]
# add factor matrix to end of batch one-hot matrix
design = pd.concat((design, factor_matrix), axis=1)
# add numerical columns to factor matrix as well...
if numerical_covariates is not None:
    sys.stderr.write("found %i numerical covariates...\n"
                        % len(numerical_covariates))
    for i, nC in enumerate(numerical_covariates):
        cname = mod.columns[nC]
        sys.stderr.write("\t{0}\n".format(cname))
        design[cname] = mod[mod.columns[nC]]
sys.stderr.write("found %i categorical variables:" % len(other_cols))
sys.stderr.write("\t" + ", ".join(other_cols) + '\n')

##########################################################
#### standarize data across genes (regions, etc..) ######
########################################################
# shape of B_hat = (design.shape[1], data.shape[1]) = (design cols, num genes)
B_hat = np.dot(np.dot(la.inv(np.dot(design.T, design)), design.T), data.T)
# shape = data.shape[1] = (num genes, )
grand_mean = np.dot((n_batches / n_array).T, B_hat[:n_batch,:])
# shape = (data.shape[1],1) = (num_genes,1)
var_pooled = np.dot(((data - np.dot(design, B_hat).T)**2), np.ones((n_array, 1)) / n_array)
# shape = data.shape = (num genes, num samples)
stand_mean = np.dot(grand_mean.T.reshape((len(grand_mean), 1)), np.ones((1, n_array)))
# copy design matrix
tmp = np.array(design.copy())
tmp[:,:n_batch] = 0
stand_mean  += np.dot(tmp, B_hat).T
# standardized data
s_data = ((data - stand_mean) / np.dot(np.sqrt(var_pooled), np.ones((1, n_array))))

######################################
### Fit L/S model and find priors ####
######################################
# one-hot of batches
batch_design = design[design.columns[:n_batch]]
# shape = (num batches, num genes)
gamma_hat = np.dot(np.dot(la.inv(np.dot(batch_design.T, batch_design)), batch_design.T), s_data.T)

# calculate delta hats for each batch
delta_hat = [] # len = num batches
for i, batch_idxs in enumerate(batch_info):
	#batches = [list(model.columns).index(b) for b in batches]
	delta_hat.append(s_data[batch_idxs].var(axis=1))


gamma_bar = gamma_hat.mean(axis=1) 
t2 = gamma_hat.var(axis=1)

a_prior = list(map(aprior, delta_hat))
b_prior = list(map(bprior, delta_hat))


######################################
#### FIND PARAMETRIC ADJUSTMENTS #####
######################################
def aprior(gamma_hat):
	m = gamma_hat.mean()
	s2 = gamma_hat.var()
	return (2 * s2 +m**2) / s2

def bprior(gamma_hat):
	m = gamma_hat.mean()
	s2 = gamma_hat.var()
	return (m*s2+m**3)/s2

def postmean(g_hat, g_bar, n, d_star, t2):
	return (t2*n*g_hat+d_star * g_bar) / (t2*n+d_star)

def postvar(sum2, n, a, b):
	return (0.5 * sum2 + b) / (n / 2.0 + a - 1.0)

def it_sol(sdat, g_hat, d_hat, g_bar, t2, a, b, conv=0.0001):
	n = (1 - np.isnan(sdat)).sum(axis=1)
	g_old = g_hat.copy()
	d_old = d_hat.copy()

	change = 1
	count = 0
	while change > conv:
		#print g_hat.shape, g_bar.shape, t2.shape
		g_new = postmean(g_hat, g_bar, n, d_old, t2)
		sum2 = ((sdat - np.dot(g_new.reshape((g_new.shape[0], 1)), np.ones((1, sdat.shape[1])))) ** 2).sum(axis=1)
		d_new = postvar(sum2, n, a, b)

		change = max((abs(g_new - g_old) / g_old).max(), (abs(d_new - d_old) / d_old).max())
		g_old = g_new #.copy()
		d_old = d_new #.copy()
		count = count + 1
	adjust = (g_new, d_new)
	return adjust 

gamma_star, delta_star = [], []
for i, batch_idxs in enumerate(batch_info):
	#print '18 20 22 28 29 31 32 33 35 40 46'
	#print batch_info[batch_id]

	temp = it_sol(s_data[batch_idxs], gamma_hat[i],
				delta_hat[i], gamma_bar[i], t2[i], a_prior[i], b_prior[i])

	gamma_star.append(temp[0])
	delta_star.append(temp[1])

############################
### Actually adjust data ###
############################
bayesdata = s_data
gamma_star = np.array(gamma_star)
delta_star = np.array(delta_star)

for j, batch_idxs in enumerate(batch_info):
	dsq = np.sqrt(delta_star[j,:])
	dsq = dsq.reshape((len(dsq), 1))
	denom =  np.dot(dsq, np.ones((1, n_batches[j])))
	numer = np.array(bayesdata[batch_idxs] - np.dot(batch_design.ix[batch_idxs], gamma_star).T)

	bayesdata[batch_idxs] = numer / denom
   
vpsq = np.sqrt(var_pooled).reshape((len(var_pooled), 1))
bayesdata = bayesdata * np.dot(vpsq, np.ones((1, n_array))) + stand_mean


return bayesdata