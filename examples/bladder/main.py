
import os
os.chdir('/users/ncullen/desktop/projects/_third_party/neuroCombat/examples/bladder')
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

