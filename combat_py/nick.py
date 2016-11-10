"""
How it works:
 	- need X data (imaging)
 	- need Y data (cognitive)
 	- pick out confound variables from Y (e.g age, sex)
 	- pick out batch variables from Y (e.g. scan site)
 	- pick out target/predictor variables from Y (e.g. depression score)


"""


import pandas as pd
import patsy
import sys
import numpy.linalg as la
import numpy as np


# (57, 5)
pheno = pd.read_table('bladder-pheno.txt', index_col=0)
# (22283, 57)
dat = pd.read_table('bladder-expr.txt', index_col=0)
# (57,4)
mod = patsy.dmatrix("~ age + cancer", pheno, return_type="dataframe")
