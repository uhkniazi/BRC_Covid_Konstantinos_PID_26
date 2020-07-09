# File: 01_EDA.R
# Auth: umar.niazi@kcl.ac.uk
# DESC: EDA and some univariate analysis of the lite data
# Date: 7/7/2020


source('header.R')

dfData = read.csv('dataExternal/gstt_raw_clinical_data_v2_lite.csv', header=T,
                  stringsAsFactors = T)
colnames(dfData)
df = dfData[,-c(5,6)]
df = stack(df)
dim(dfData)
length(unique(df$ind))
library(lattice)
df$resp = dfData$hospital_outcome
df$resp2 = dfData$outcome_28d

xyplot(resp ~ values | ind, data=df, type=c('g', 'p'), pch=19, cex=0.6,
       par.strip.text=list(cex=0.7), 
       scales = list(x=list(rot=45, cex=0.4), relation='free'),
       ylab='Hospital Outcome')

xyplot(resp2 ~ values | ind, data=df, type=c('g', 'p'), pch=19, cex=0.6,
       par.strip.text=list(cex=0.7), 
       scales = list(x=list(rot=45, cex=0.4), relation='free'),
       ylab='28 Day outcome')

### check for overlap and balance
cn = colnames(dfData)[-c(5,6)]
f = dfData$hospital_outcome

pdf('temp/overlaps.pdf')
par(mfrow=c(2,2))
for (i in seq_along(cn)){
  x = as.numeric(dfData[,cn[i]])
  hist2(x[f == 'ALIVE'], x[f=='DIED'], main=cn[i], legends=c('A', 'D'))
}
dev.off(dev.cur())
################################################
### some quick variable selection and model fitting
### for prediction purposes
################################################
if(!require(downloader) || !require(methods)) stop('Library downloader and methods required')

url = 'https://raw.githubusercontent.com/uhkniazi/CCrossValidation/experimental/CCrossValidation.R'
download(url, 'CCrossValidation.R')

# load the required packages
source('CCrossValidation.R')
# delete the file after source
unlink('CCrossValidation.R')

## prepare data to use in the analysis
dfData.bk = dfData
# drop one of the age variables repeated twice in input sheet
grep('age', colnames(dfData))
dfData = dfData[,-18]

x = rep(NA, length=ncol(dfData))
for(i in 1:ncol(dfData)){
  x[i] = sum(is.na(dfData[,i]))
}
# remove variables with large NAs
data.frame(colnames(dfData), x)
which(x > 10)
dfData = dfData[,-c(40, 45, 46)]
table(complete.cases(dfData))

x = rep(NA, length=ncol(dfData))
for(i in 1:ncol(dfData)){
  x[i] = sum(is.na(dfData[,i]))
}
# remove variables with large NAs
data.frame(colnames(dfData), x)
x = which(x > 0)
cn = colnames(dfData)[x]

## impute the average values for missing data 
library(Hmisc)
for(i in 1:length(cn)){
  dfData[,cn[i]] = impute(dfData[,cn[i]], mean)
}
colnames(dfData)

## remove variables with 0 sd i.e. not changing 
s = apply(dfData, 2, sd)
s = which(s == 0)

lData.train = list(data=dfData[,-c(5,6, s)], covariates=dfData[,c(5,6)])

## drop certain variables
i = grep('days_ad', colnames(lData.train$data))
lData.train$data = lData.train$data[,-i]
########################## perform a random forest step
dfData = lData.train$data
fGroups = lData.train$covariates$hospital_outcome

oVar.r = CVariableSelection.RandomForest(dfData, fGroups, boot.num = 100)

plot.var.selection(oVar.r)

######################## Stan section for binomial regression approach
dfData = data.frame(lData.train$data)
dim(dfData)
dfData$fGroups = lData.train$covariates$hospital_outcome
rm(fGroups)
levels(dfData$fGroups)
dfData$sex = as.numeric(dfData$sex)-1
lData = list(resp=ifelse(dfData$fGroups == 'DIED', 1, 0), mModMatrix=model.matrix(fGroups ~ 1 + ., data=dfData))

library(rethinking)
library(rstan)
rstan_options(auto_write = TRUE)
options(mc.cores = parallel::detectCores())
setwd(gcswd)
stanDso = rstan::stan_model(file='binomialRegressionSharedCoeffVariance.stan')

lStanData = list(Ntotal=length(lData$resp), Ncol=ncol(lData$mModMatrix), X=lData$mModMatrix,
                 y=lData$resp)

# ## give initial values
# initf = function(chain_id = 1) {
#   list(betas=rep(0, times=ncol(lStanData$X)), tau=0.5)
# }


fit.stan = sampling(stanDso, data=lStanData, iter=2000, chains=4, pars=c('tau', 'betas2', 'log_lik'), cores=4,# init=initf,
                    control=list(adapt_delta=0.99, max_treedepth = 13))

#save(fit.stan, file='temp/fit.stan.binom_guess.rds')

print(fit.stan, c('betas2', 'tau'))
print(fit.stan, 'tau')
traceplot(fit.stan, 'tau')

## get the coefficient of interest
mCoef = extract(fit.stan)$betas2
dim(mCoef)
## get the intercept
iIntercept = mCoef[,1]
mCoef = mCoef[,-1]
colnames(mCoef) = colnames(lData$mModMatrix)[2:ncol(lData$mModMatrix)]

## coeftab object 
ct.1 = coeftab(fit.stan)
rn = rownames(ct.1@coefs)
i = grep('betas', rn)
rownames(ct.1@coefs)[i[-1]] = colnames(mCoef)
rownames(ct.1@se)[i[-1]] = colnames(mCoef)
plot(ct.1, pars=colnames(mCoef))

# ## find correlated variables
dim(dfData)
mData = as.matrix(dfData[,-36])
length(as.vector(mData))
mCor = cor((mData+runif(2170, 1e-4, 1e-3)), use="na.or.complete")
library(caret)
### find the columns that are correlated and should be removed
n = findCorrelation((mCor), cutoff = 0.8, names=T)
data.frame(n)
sapply(n, function(x) {
  (abs(mCor[,x]) >= 0.8)
})

n = sapply(n, function(x) {
  rownames(mCor)[(abs(mCor[,x]) >= 0.8)]
})

n = c(n[1,], n[2,])

pairs(mCoef[sample(1:nrow(mCoef), 500),n], col='grey', cex=0.5)

## to drop
i = c(grep('Creat', colnames(lData.train$data)),
      grep('haema', colnames(lData.train$data)),
      grep('Physio', colnames(lData.train$data)),
      grep('WCC|chronic_h|age_sc', colnames(lData.train$data)))
colnames(lData.train$data[,i])
lData.train$data = lData.train$data[,-i]

######### second stan model with a smaller number of variables
dfData = data.frame(lData.train$data)
dim(dfData)
dfData$fGroups = lData.train$covariates$hospital_outcome
rm(fGroups)
levels(dfData$fGroups)
dfData$sex = as.numeric(dfData$sex)-1
lData = list(resp=ifelse(dfData$fGroups == 'DIED', 1, 0), mModMatrix=model.matrix(fGroups ~ 1 + ., data=dfData))

lStanData = list(Ntotal=length(lData$resp), Ncol=ncol(lData$mModMatrix), X=lData$mModMatrix,
                 y=lData$resp)

# ## give initial values
# initf = function(chain_id = 1) {
#   list(betas=rep(0, times=ncol(lStanData$X)), tau=0.5)
# }


fit.stan.2 = sampling(stanDso, data=lStanData, iter=2000, chains=4, pars=c('tau', 'betas2', 'log_lik'), cores=4,# init=initf,
                    control=list(adapt_delta=0.99, max_treedepth = 13))

#save(fit.stan, file='temp/fit.stan.binom_guess.rds')

print(fit.stan.2, c('betas2', 'tau'))

## get the coefficient of interest
mCoef = extract(fit.stan.2)$betas2
dim(mCoef)
## get the intercept
iIntercept = mCoef[,1]
mCoef = mCoef[,-1]
colnames(mCoef) = colnames(lData$mModMatrix)[2:ncol(lData$mModMatrix)]

## coeftab object 
ct.1 = coeftab(fit.stan.2)
rn = rownames(ct.1@coefs)
i = grep('betas', rn)
rownames(ct.1@coefs)[i[-1]] = colnames(mCoef)
rownames(ct.1@se)[i[-1]] = colnames(mCoef)
plot(ct.1, pars=colnames(mCoef))

compare(fit.stan, fit.stan.2)
plot(compare(fit.stan, fit.stan.2))

### predictive performance of the model
## binomial prediction
mypred = function(theta, data){
  betas = theta # vector of betas i.e. regression coefficients for population
  ## data
  mModMatrix = data$mModMatrix
  # calculate fitted value
  iFitted = mModMatrix %*% betas
  # using logit link so use inverse logit
  #iFitted = plogis(iFitted)
  return(iFitted)
}


mCoef = extract(fit.stan.2)$betas2
dim(mCoef)
colnames(mCoef) = c('Intercept', colnames(lData$mModMatrix)[2:ncol(lData$mModMatrix)])
library(lattice)
## get the predicted values
## create model matrix
X = as.matrix(cbind(rep(1, times=nrow(dfData)), dfData[,colnames(mCoef)[-1]]))
colnames(X) = colnames(mCoef)
head(X)
ivPredict = plogis(mypred(colMeans(mCoef), list(mModMatrix=X))[,1])
xyplot(ivPredict ~ fGroups, xlab='Actual Group', 
       ylab='Predicted Probability of Died',
       data=dfData)
# outlier samples
i = which(ivPredict < 0.5 & dfData$fGroups == 'DIED')
ivOutliers = i




# ## top variables according to random forest and binomial models
# dfRF = CVariableSelection.RandomForest.getVariables(oVar.r)
# # select the top 30 variables
# cvTopVars = rownames(dfRF)[1:20]
# cvTopVars = cvTopVars[!cvTopVars %in% i]
# # top from binomial model
# m = colMeans(mCoef)
# m = abs(m)
# m = sort(m, decreasing = T)
# cvTopVars.binomial = names(m)[1:20] #names(m[m > 0.25])
# cvTopVars.binomial = cvTopVars.binomial[!cvTopVars.binomial %in% i]
# dfData = lData.train$data
# fGroups = lData.train$covariates$hospital_outcome
# dim(dfData)
# 
# oVar.sub = CVariableSelection.ReduceModel(dfData[,cvTopVars.binomial], fGroups, boot.num = 100)
# 
# # plot the number of variables vs average error rate
# plot.var.selection(oVar.sub)
# 
# # use the top 30 genes to find top combinations of genes
# dfData = data.frame(t(lData.train$data[cvTopGenes.binomial, ]))
# 
# oVar.sub2 = CVariableSelection.ReduceModel(dfData, fGroups, boot.num = 100)
# 
# # plot the number of variables vs average error rate
# plot.var.selection(oVar.sub2)