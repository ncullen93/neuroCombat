



pheno = pd.read_table('bladder-pheno.txt', index_col=0) # Y (cognitive) data
dat = pd.read_table('bladder-expr.txt', index_col=0) # X (imaging) data

mod = patsy.dmatrix("~ 1", pheno, return_type="dataframe")

ebat = combat(dat, pheno.batch, mod, "age")