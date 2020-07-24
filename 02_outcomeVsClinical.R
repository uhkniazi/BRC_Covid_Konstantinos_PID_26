# File: 02_outcomeVsClinical.R
# Auth: umar.niazi@kcl.ac.uk
# DESC: Model to show the performance of the outcome vs selected clinical variables
# Date: 20/7/2020

source('header.R')

load('temp/lData.train.rds')

library(rethinking)
library(rstan)
rstan_options(auto_write = TRUE)
options(mc.cores = parallel::detectCores())

library(lattice)

colnames(lData.train$data)
dfData = lData.train$data[,c(1, 3, 5, 7)]
dfData$outcome_28d = lData.train$covariates$outcome_28d

xyplot(outcome_28d ~ Age | factor(Viremia), data=dfData, type=c('g', 'r', 'p'),
       auto.key=T)

xyplot(outcome_28d ~ BMI | factor(Viremia), data=dfData, type=c('g', 'r', 'p'),
       auto.key=T)

xyplot(outcome_28d ~ days_2Onset | factor(Viremia), data=dfData, type=c('g', 'r', 'p'),
       auto.key=T)


xyplot(outcome_28d ~ Age | factor(Viremia), data=dfData, type=c('g', 'smooth', 'p'),
       auto.key=T)

xyplot(outcome_28d ~ Age, data=dfData, type=c('g', 'r', 'p'),
       auto.key=T)

xyplot(outcome_28d ~ Age, data=dfData, type=c('g', 'smooth', 'p'),
       auto.key=T)

## see page 332 rethinking book for priors and model description
#dfData$Viremia = factor(dfData$Viremia)
fit.1 = glm(outcome_28d ~ Age*Viremia, data=dfData, family='binomial')
fit.2 = glm(outcome_28d ~ Age + Viremia, data=dfData, family='binomial')

anova(fit.1, fit.2)

m1 = model.matrix(outcome_28d ~ ., data=dfData)
#m2 = model.matrix(outcome_28d ~ Age:Viremia -1, data=dfData)
lData = list(resp=ifelse(dfData$outcome_28d == '1', 1, 0), 
             mModMatrix=m1)
stanDso = rstan::stan_model(file='binomial.stan')

lStanData = list(Ntotal=length(lData$resp), Ncol=ncol(lData$mModMatrix), X=lData$mModMatrix,
                 y=lData$resp)

fit.stan.1 = sampling(stanDso, data=lStanData, iter=5000, chains=4, pars=c('betas', 'log_lik'), cores=4)

print(fit.stan.1, pars=c('betas'))
traceplot(fit.stan.1, pars=c('betas'))

# plot coefficients
ct.1 = coeftab(fit.stan.1)
rn = rownames(ct.1@coefs)
i = grep('betas', rn)
rownames(ct.1@coefs)[i[-1]] = colnames(lData$mModMatrix)[2:ncol(lData$mModMatrix)]
rownames(ct.1@se)[i[-1]] = colnames(lData$mModMatrix)[2:ncol(lData$mModMatrix)]
plot(ct.1, pars=colnames(lData$mModMatrix)[2:ncol(lData$mModMatrix)])
ct.1.1 = ct.1

## second formulation with only viremia
m1 = model.matrix(outcome_28d ~ Age -1, data=dfData)
m2 = model.matrix(outcome_28d ~ Viremia -1, data=dfData)
lData = list(resp=ifelse(dfData$outcome_28d == '1', 1, 0), 
             mModMatrix=cbind(1, m1, m2))

lStanData = list(Ntotal=length(lData$resp), Ncol=ncol(lData$mModMatrix), X=lData$mModMatrix,
                 y=lData$resp)

fit.stan.2 = sampling(stanDso, data=lStanData, iter=5000, chains=4, pars=c('betas', 'log_lik'), cores=4)

print(fit.stan.2, pars=c('betas'))

# plot coefficients
ct.1 = coeftab(fit.stan.2)
rn = rownames(ct.1@coefs)
i = grep('betas', rn)
rownames(ct.1@coefs)[i[-1]] = colnames(lData$mModMatrix)[2:ncol(lData$mModMatrix)]
rownames(ct.1@se)[i[-1]] = colnames(lData$mModMatrix)[2:ncol(lData$mModMatrix)]
plot(ct.1, pars=colnames(lData$mModMatrix)[2:ncol(lData$mModMatrix)])

compare(fit.stan.1, fit.stan.2)
plot(compare(fit.stan.1, fit.stan.2))

############# performance of model 1 vs model 2
if(!require(downloader) || !require(methods)) stop('Library downloader and methods required')

url = 'https://raw.githubusercontent.com/uhkniazi/CCrossValidation/experimental/CCrossValidation.R'
download(url, 'CCrossValidation.R')

# load the required packages
source('CCrossValidation.R')
# delete the file after source
unlink('CCrossValidation.R')

rm(stanDso)
colnames(dfData)

oCV.1 = CCrossValidation.StanBern(dfData[,-5], dfData[,-5], dfData$outcome_28d,
                                  dfData$outcome_28d, level.predict = '1',
                                  boot.num = 10, k.fold = 10, ncores = 2, nchains = 2) 

plot.cv.performance(oCV.1)

## try second model
oCV.2 = CCrossValidation.StanBern(dfData[,-c(3:5)], dfData[,-c(3:5)], dfData$outcome_28d,
                                  dfData$outcome_28d, level.predict = '1',
                                  boot.num = 10, k.fold = 10, ncores = 2, nchains = 2) 

plot.cv.performance(oCV.2)


############### work with model 2 as it is simpler and similar to model 1
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

## setup the data
dfData = lData.train$data[,c(1, 3)]
dfData$outcome_28d = lData.train$covariates$outcome_28d
str(dfData)
m1 = model.matrix(outcome_28d ~ Age -1, data=dfData)
m2 = model.matrix(outcome_28d ~ Viremia -1, data=dfData)
lData = list(resp=ifelse(dfData$outcome_28d == '1', 1, 0), 
             mModMatrix=cbind(1, m1, m2))


mCoef = extract(fit.stan.2)$betas
dim(mCoef)
colnames(mCoef) = c('Intercept', colnames(lData$mModMatrix)[2:ncol(lData$mModMatrix)])
head(mCoef)
library(lattice)
## get the predicted values
## create model matrix
X = model.matrix(outcome_28d ~ ., data=dfData)
head(X)

ivPredict = plogis(mypred(colMeans(mCoef), list(mModMatrix=X))[,1])
xyplot(ivPredict ~ dfData$outcome_28d, xlab='Actual Group', 
       ylab='Predicted Probability of Died',
       data=dfData)

library(ROCR)
## draw a ROC curve first for calibration performance test
ivTruth = dfData$outcome_28d == '1'
p = prediction(ivPredict, ivTruth)
perf.alive = performance(p, 'tpr', 'fpr')
dfPerf.alive = data.frame(c=perf.alive@alpha.values, t=perf.alive@y.values[[1]], f=perf.alive@x.values[[1]], 
                          r=perf.alive@y.values[[1]]/perf.alive@x.values[[1]])
colnames(dfPerf.alive) = c('c', 't', 'f', 'r')
plot(perf.alive)

## show prediction performance based on cutoff
cPredict = rep('1', times=length(ivPredict))
cPredict[ivPredict < 0.25] = '0'
table(truth=dfData$outcome_28d, predicted=cPredict)

## draw some counterfactual plots
## create the plots for regression with each predictor (not input) fixed at its average
## see Data Analysis ... Regression & Multilevel M [Gelman] for jitter.binary function
jitter.binary = function(a, jitt=.05){
  ifelse (a==0, runif (length(a), 0, jitt), runif (length(a), 1-jitt, 1))
}

outcome.jitter = jitter.binary(dfData$outcome_28d)

plot(dfData$Age, outcome.jitter, pch=20, 
     xlab='Affect of increasing Age', ylab='Probability of Death Class',
     main='Prediction of Death class vs Age')
x = seq(min(dfData$Age), max(dfData$Age), length.out = 100)
m = cbind(1, x, 0)
c = colMeans(mCoef)
lines(x, plogis(m %*% c), col='black')
m = cbind(1, x, 1)
lines(x, plogis(m %*% c), col='red')

legend('topleft', legend = c('Vir=0', 'Vir=1'), fill=c('black', 'red'))

####### repeat with model 1
dfData = lData.train$data[,c(1, 3, 5, 7)]
dfData$outcome_28d = lData.train$covariates$outcome_28d
str(dfData)

m1 = model.matrix(outcome_28d ~ ., data=dfData)
lData = list(resp=ifelse(dfData$outcome_28d == '1', 1, 0), 
             mModMatrix=m1)

mCoef = extract(fit.stan.1)$betas
dim(mCoef)
colnames(mCoef) = c('Intercept', colnames(lData$mModMatrix)[2:ncol(lData$mModMatrix)])
head(mCoef)

outcome.jitter = jitter.binary(dfData$outcome_28d)

plot(dfData$Age, outcome.jitter, pch=20, 
     xlab='Affect of increasing Age', ylab='Probability of Death Class',
     main='Prediction of Death class vs Age')
x = seq(min(dfData$Age), max(dfData$Age), length.out = 100)
colnames(mCoef)
m = cbind(1, x, 0, mean(dfData$BMI), mean(dfData$days_2Onset))
c = colMeans(mCoef)
lines(x, plogis(m %*% c), col='black')
m = cbind(1, x, 1, mean(dfData$BMI), mean(dfData$days_2Onset))
lines(x, plogis(m %*% c), col='red')

legend('topleft', legend = c('Vir=0', 'Vir=1'), fill=c('black', 'red'))

### BMI
plot(dfData$BMI, outcome.jitter, pch=20, 
     xlab='Affect of increasing BMI', ylab='Probability of Death Class',
     main='Prediction of Death class vs BMI')
x = seq(min(dfData$BMI), max(dfData$BMI), length.out = 100)
colnames(mCoef)
m = cbind(1, mean(dfData$Age), 0, x, mean(dfData$days_2Onset))
c = colMeans(mCoef)
lines(x, plogis(m %*% c), col='black')
m = cbind(1, mean(dfData$Age), 1, x, mean(dfData$days_2Onset))
lines(x, plogis(m %*% c), col='red')

legend('topleft', legend = c('Vir=0', 'Vir=1'), fill=c('black', 'red'))


### days 2 onset
plot(dfData$days_2Onset, outcome.jitter, pch=20, 
     xlab='Affect of increasing days2Onset', ylab='Probability of Death Class',
     main='Prediction of Death class vs days2Onset')
x = seq(min(dfData$days_2Onset), max(dfData$days_2Onset), length.out = 100)
colnames(mCoef)
m = cbind(1, mean(dfData$Age), 0, mean(dfData$BMI), x)
c = colMeans(mCoef)
lines(x, plogis(m %*% c), col='black')
m = cbind(1, mean(dfData$Age), 1, mean(dfData$BMI), x)
lines(x, plogis(m %*% c), col='red')

legend('bottomright', legend = c('Vir=0', 'Vir=1'), fill=c('black', 'red'))
