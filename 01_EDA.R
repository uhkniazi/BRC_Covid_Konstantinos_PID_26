# File: 01_EDA.R
# Auth: umar.niazi@kcl.ac.uk
# DESC: EDA and some analysis of the merged data gstt and kch (clinical)
# Date: 15/7/2020


source('header.R')

dfData = read.csv('dataExternal/combined_dataset_covid_predictive_analytics.csv', header=T,
                  stringsAsFactors = T)
colnames(dfData)
dfData = dfData[,-c(1)]
colnames(dfData)
str(dfData)
summary(dfData)
dfData$outcome_28d = factor(dfData$outcome_28d)
df = dfData[,-c(31)]
df = stack(df)
dim(dfData)
length(unique(df$ind))
library(lattice)
df$resp = dfData$outcome_28d

xyplot(resp ~ values | ind, data=df, type=c('g', 'p'), pch=19, cex=0.6,
       par.strip.text=list(cex=0.7), 
       scales = list(x=list(rot=45, cex=0.4), relation='free'),
       ylab='28 Day outcome')

### check for overlap and balance
cn = colnames(dfData)[-c(31)]
f = dfData$outcome_28d

pdf('temp/overlaps_d28.pdf')
par(mfrow=c(2,2))
for (i in seq_along(cn)){
  x = as.numeric(dfData[,cn[i]])
  hist2(x[f == '0'], x[f=='1'], main=cn[i], legends=c('A', 'D'))
}
dev.off(dev.cur())

xyplot(outcome_28d ~ Age, groups=Viremia, data=dfData, type=c('g', 'r', 'p'),
       auto.key=T)

xyplot(outcome_28d ~ FiO2, groups=Viremia, data=dfData, type=c('g', 'smooth', 'p'),
       auto.key=T)



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

table(complete.cases(dfData))
x = rep(NA, length=ncol(dfData))
for(i in 1:ncol(dfData)){
  x[i] = sum(is.na(dfData[,i]))
}
# remove variables with large NAs
data.frame(colnames(dfData), x)
i = which(x > 2)
tapply(dfData[,'Urea'], dfData$outcome_28d, function(x) sum(is.na(x)))

dfData = dfData[,-c(i)]
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
cvImputed = cn

## remove variables with 0 sd i.e. not changing 
s = apply(dfData, 2, sd)
s = which(s == 0)

as.data.frame(colnames(dfData))
lData.train = list(data=dfData[,-c(19, 20)], covariates=dfData[,c(19, 20)])

colnames(lData.train$data)

cn = colnames(dfData)[-c(2, 4, 12, 18:20)]

pdf('temp/xyplots.pdf')
for (i in seq_along(cn)){
  print(xyplot(outcome_28d ~ dfData[,cn[i]], groups=Viremia, data=dfData, type=c('g', 'r', 'p'),
               auto.key=T, xlab=cn[i]))
  
  print(xyplot(outcome_28d ~ dfData[,cn[i]], groups=Viremia, data=dfData, type=c('g', 'smooth', 'p'),
               auto.key=T, xlab=cn[i]))
  
}
dev.off(dev.cur())


########################## perform a random forest step
dfData = lData.train$data
fGroups = lData.train$covariates$outcome_28d

oVar.r = CVariableSelection.RandomForest(dfData, fGroups, boot.num = 100)

plot.var.selection(oVar.r)

######################## Stan section for binomial regression approach
dfData = data.frame(lData.train$data)
dim(dfData)
dfData$fGroups = lData.train$covariates$outcome_28d
table(dfData$fGroups)
rm(fGroups)
levels(dfData$fGroups)
#dfData$sex = as.numeric(dfData$sex)-1
lData = list(resp=ifelse(dfData$fGroups == '1', 1, 0), mModMatrix=model.matrix(fGroups ~ 1 + ., data=dfData))

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


fit.stan = sampling(stanDso, data=lStanData, iter=1000, chains=4, pars=c('tau', 'betas2', 'log_lik'), cores=4,# init=initf,
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

m = abs(colMeans(mCoef))
m = sort(m, decreasing = T)

l2 = barplot(m, 
             las=2, xaxt='n', col='grey', main='Top Variables')
axis(1, at = l2, labels = names(m), tick = F, las=2, cex.axis=0.7 )

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


mCoef = extract(fit.stan)$betas2
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


# ## find correlated variables
dim(dfData)
mData = as.matrix(dfData[,-19])
length(as.vector(mData))
mCor = cor((mData+runif(1425, 1e-4, 1e-3)), use="na.or.complete")
library(caret)
image(mCor)
### find the columns that are correlated and should be removed
n = findCorrelation((mCor), cutoff = 0.7, names=T)
data.frame(n)
# sapply(n, function(x) {
#   (abs(mCor[,x]) >= 0.7)
# })
# 
n = sapply(n, function(x) {
  rownames(mCor)[(abs(mCor[,x]) >= 0.7)]
})
# 
# n = c(n[1,], n[2,])
# 
# pairs(mCoef[sample(1:nrow(mCoef), 500),n], col='grey', cex=0.5)
# 
# ## to drop
# i = c(grep('Age_score', colnames(lData.train$data)))
#       # grep('age_score', colnames(lData.train$data)))
#       # grep('Physio', colnames(lData.train$data)),
#       # grep('WCC|age_sc|miR.451a.5p|miR.192.5p', colnames(lData.train$data)))
# colnames(lData.train$data[,i])
# lData.train$data = lData.train$data[,-i]

cvTopVariables.rf = rownames(CVariableSelection.RandomForest.getVariables(oVar.r))[1:8]
cvTopVariables.bin = names(m)[1:6]
cvTopVariables = unique(c(cvTopVariables.rf, cvTopVariables.bin))

lData.train$data$Age_vir = lData.train$data$Age * lData.train$data$Viremia
lData.train$data$Age_vir2 = lData.train$data$Age_vir^2
cvTopVariables = cvTopVariables[-4]
lData.train$data = lData.train$data[,cvTopVariables]
########################## perform a second random forest step
dfData = lData.train$data
fGroups = lData.train$covariates$outcome_28d

oVar.r = CVariableSelection.RandomForest(dfData, fGroups, boot.num = 100)

plot.var.selection(oVar.r)

######### second stan model with a smaller number of variables
dfData = data.frame(lData.train$data)
dim(dfData)
dfData$fGroups = lData.train$covariates$outcome_28d
rm(fGroups)
levels(dfData$fGroups)

table(dfData$fGroups)
# i = which(dfData$fGroups == 'DIED')
# i2 = which(dfData$fGroups == 'ALIVE')
# i = sample(i, 10)
# i2 = sample(i2, 10)
# dfData = dfData[c(i, i2),]
levels(dfData$fGroups)

lData = list(resp=ifelse(dfData$fGroups == '1', 1, 0), mModMatrix=model.matrix(fGroups ~ 1 + ., data=dfData))

lStanData = list(Ntotal=length(lData$resp), Ncol=ncol(lData$mModMatrix), X=lData$mModMatrix,
                 y=lData$resp)
stanDso = rstan::stan_model(file='binomialRegressionSharedCoeffVariance.stan')

## give initial values
initf = function(chain_id = 1) {
  list(betas=rep(0, times=ncol(lStanData$X)), tau=0.5)
}


fit.stan.2 = sampling(stanDso, data=lStanData, iter=1000, chains=4, pars=c('tau', 'betas2', 'log_lik'), cores=4, init=initf,
                    control=list(adapt_delta=0.99, max_treedepth = 13))

#save(fit.stan, file='temp/fit.stan.binom_guess.rds')

print(fit.stan.2, c('betas2', 'tau'))
traceplot(fit.stan.2, 'tau')
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

m = abs(colMeans(mCoef))
m = sort(m, decreasing = T)

l2 = barplot(m, 
             las=2, xaxt='n', col='grey', main='Top Variables')
axis(1, at = l2, labels = names(m), tick = F, las=2, cex.axis=0.7 )


# compare(fit.stan, fit.stan.2)
# plot(compare(fit.stan, fit.stan.2))

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


## draw a ROC curve first for calibration performance test
ivTruth = dfData$fGroups == '1'
p = prediction(ivPredict, ivTruth)
perf.alive = performance(p, 'tpr', 'fpr')
dfPerf.alive = data.frame(c=perf.alive@alpha.values, t=perf.alive@y.values[[1]], f=perf.alive@x.values[[1]], 
                          r=perf.alive@y.values[[1]]/perf.alive@x.values[[1]])
colnames(dfPerf.alive) = c('c', 't', 'f', 'r')
plot(perf.alive)

cPredict = rep('1', times=length(ivPredict))
cPredict[ivPredict < 0.2] = '0'
table(dfData$fGroups, cPredict)


##############################
### due to presence of outliers try the model with a guess parameter
dfData = data.frame(lData.train$data)
dim(dfData)
dfData$fGroups = lData.train$covariates$outcome_28d
rm(fGroups)
levels(dfData$fGroups)

table(dfData$fGroups)
# i = which(dfData$fGroups == 'DIED')
# i2 = which(dfData$fGroups == 'ALIVE')
# i = sample(i, 10)
# i2 = sample(i2, 10)
# dfData = dfData[c(i, i2),]
# table(dfData$fGroups)

#dfData$sex = as.numeric(dfData$sex)-1
lData = list(resp=ifelse(dfData$fGroups == '1', 1, 0), mModMatrix=model.matrix(fGroups ~ 1 + ., data=dfData))

stanDso = rstan::stan_model(file='binomialGuessMixtureRegressionSharedCoeffVariance.stan')

lStanData = list(Ntotal=length(lData$resp), Ncol=ncol(lData$mModMatrix), X=lData$mModMatrix,
                 y=lData$resp)

## give initial values
initf = function(chain_id = 1) {
  list(betas=rep(0, times=ncol(lStanData$X)), tau=0.5)
}


fit.stan.3 = sampling(stanDso, data=lStanData, iter=1000, chains=4, pars=c('tau', 'betas2', 'log_lik'), cores=4, init=initf,
                    control=list(adapt_delta=0.99, max_treedepth = 13))

#save(fit.stan, file='temp/fit.stan.binom_guess.rds')

print(fit.stan.3, c('betas2', 'tau'))
print(fit.stan.3, 'tau')
traceplot(fit.stan.3, 'tau')

## get the coefficient of interest
mCoef = extract(fit.stan.3)$betas2
dim(mCoef)
## get the intercept
iIntercept = mCoef[,1]
mCoef = mCoef[,-1]
colnames(mCoef) = colnames(lData$mModMatrix)[2:ncol(lData$mModMatrix)]

## coeftab object 
ct.1 = coeftab(fit.stan.3)
rn = rownames(ct.1@coefs)
i = grep('betas', rn)
rownames(ct.1@coefs)[i[-1]] = colnames(mCoef)
rownames(ct.1@se)[i[-1]] = colnames(mCoef)
plot(ct.1, pars=colnames(mCoef))

compare(fit.stan, fit.stan.2, fit.stan.3)
plot(compare(fit.stan, fit.stan.2, fit.stan.3))

mCoef = extract(fit.stan.3)$betas2
dim(mCoef)
colnames(mCoef) = c('Intercept', colnames(lData$mModMatrix)[2:ncol(lData$mModMatrix)])
library(lattice)
## get the predicted values
## create model matrix
dfData = data.frame(lData.train$data)
dim(dfData)
dfData$fGroups = lData.train$covariates$outcome_28d

X = as.matrix(cbind(rep(1, times=nrow(dfData)), dfData[,colnames(mCoef)[-1]]))
colnames(X) = colnames(mCoef)
head(X)
ivPredict = plogis(mypred(colMeans(mCoef), list(mModMatrix=X))[,1])
xyplot(ivPredict ~ fGroups, xlab='Actual Group', 
       ylab='Predicted Probability of Died',
       data=dfData)
# # outlier samples
i = which(ivPredict < 0.5 & dfData$fGroups == '1')
# p = rep('ALIVE', times=length(dfData$fGroups))
# p[ivPredict > 0.5] = 'DIED'
# table(dfData$fGroups, p)

############################################################################
######## addressing imbalance
############################################################################
stanDso = rstan::stan_model(file='binomialRegressionSharedCoeffVariance.stan')

dfData = data.frame(lData.train$data)
dim(dfData)
dfData$fGroups = lData.train$covariates$outcome_28d
rm(fGroups)
levels(dfData$fGroups)

table(dfData$fGroups)

lFits = vector(mode = 'list', length=10)
for(b in 1:10){
  dfData = data.frame(lData.train$data)
  dfData$fGroups = lData.train$covariates$outcome_28d
  i = which(dfData$fGroups == '1')
  i2 = which(dfData$fGroups == '0')
  i = sample(i, 15)
  i2 = sample(i2, 30)
  dfData = dfData[c(i, i2),]
  
  
  lData = list(resp=ifelse(dfData$fGroups == '1', 1, 0), mModMatrix=model.matrix(fGroups ~ 1 + ., data=dfData))
  
  lStanData = list(Ntotal=length(lData$resp), Ncol=ncol(lData$mModMatrix), X=lData$mModMatrix,
                   y=lData$resp)
  
  lFits[[b]] = sampling(stanDso, data=lStanData, iter=1000, chains=4, pars=c('tau', 'betas2', 'log_lik'), cores=4, init=initf,
                        control=list(adapt_delta=0.99, max_treedepth = 13))
  
  ## get the coefficient of interest
  mCoef = extract(lFits[[b]])$betas2
  dim(mCoef)
  ## get the intercept
  iIntercept = mCoef[,1]
  mCoef = mCoef[,-1]
  colnames(mCoef) = colnames(lData$mModMatrix)[2:ncol(lData$mModMatrix)]
  
  ## coeftab object 
  ct.1 = coeftab(lFits[[b]])
  rn = rownames(ct.1@coefs)
  i = grep('betas', rn)
  rownames(ct.1@coefs)[i[-1]] = colnames(mCoef)
  rownames(ct.1@se)[i[-1]] = colnames(mCoef)
  plot(ct.1, pars=colnames(mCoef))
}

ct.1 = coeftab(lFits[[1]], lFits[[2]], lFits[[3]], 
               lFits[[4]], lFits[[5]], lFits[[6]],
               lFits[[7]], lFits[[8]], lFits[[9]],
               lFits[[10]])
rn = rownames(ct.1@coefs)
i = grep('betas', rn)
rownames(ct.1@coefs)[i[-1]] = colnames(mCoef)
rownames(ct.1@se)[i[-1]] = colnames(mCoef)
plot(ct.1, pars=c(colnames(mCoef)))

## select variables to drop
# drop some variables that show consistent zero coefficients or opposite under 
# sub-sampling - see the loop above where this was done to select
# the variables after fitting multiple sub-sampled models
i = grep('Crea|Fi|Hema', cvTopVariables)
cvTopVariables[i]
cvTopVariables = cvTopVariables[-i]

dfData = lData.train$data
fGroups = lData.train$covariates$outcome_28d
dim(dfData)

oVar.sub = CVariableSelection.ReduceModel(dfData, fGroups, boot.num = 200)

# plot the number of variables vs average error rate
plot.var.selection(oVar.sub)

cvVar = CVariableSelection.ReduceModel.getMinModel(oVar.sub, size = 2)


oCV.lda = CCrossValidation.LDA(dfData[,cvVar], dfData[,cvVar], fGroups, fGroups, level.predict = '1',
                               boot.num = 100, k.fold = 10) 

plot.cv.performance(oCV.lda)

rm(stanDso)
url = 'https://raw.githubusercontent.com/uhkniazi/CCrossValidation/experimental/bernoulli.stan'
download(url, 'bernoulli.stan')

oCV.s = CCrossValidation.StanBern(dfData[,cvVar], dfData[,cvVar], fGroups, fGroups, level.predict = '1',
                                  boot.num = 15, k.fold = 10, ncores = 2, nchains = 2) 

save(oCV.s, file='temp/oCV.stan.rds')

plot.cv.performance(oCV.s)
