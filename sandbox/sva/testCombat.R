# imports
library(sva)
library(bladderbatch)
data(bladderdata)

# load data
pheno <- pData(bladderEset)
edata <- exprs(bladderEset)

# create models
batch <- pheno$batch
pheno$age = c(1:7, rep(1:10, 5))
modcombat <- model.matrix(~ age + as.factor(cancer),data=pheno)

c2 <- ComBat(dat=edata,batch=batch,mod=modcombat)
