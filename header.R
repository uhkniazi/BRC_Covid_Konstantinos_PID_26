# File: header.R
# Auth: umar.niazi@kcl.ac.uk
# DESC: global variables
# Date: 07/17/2020


## variables
g_pid = 26
#g_did = c(36, 37, 40)
gcswd = getwd()

p.old = par()

###### utility functions

f_LoadObject = function(r.obj.file)
{
  # temp environment to load object
  env <- new.env()
  # read object
  nm <- load(r.obj.file, env)[1]
  return(env[[nm]])
}


# utility function to calculate gamma distribution used as a prior for scale
gammaShRaFromModeSD = function( mode , sd ) {
  # function changed a little to return jeffery non-informative prior
  if ( mode <=0 || sd <=0 ) return( list( shape=0.5 , rate=0.0001 ) ) 
  rate = ( mode + sqrt( mode^2 + 4 * sd^2 ) ) / ( 2 * sd^2 )
  shape = 1 + mode * rate
  return( list( shape=shape , rate=rate ) )
}

f_plotVolcano = function(dfGenes, main, p.adj.cut = 0.1, fc.lim = c(-3, 3)){
  p.val = -1 * log10(dfGenes$P.Value)
  fc = dfGenes$logFC
  # cutoff for p.value y.axis
  y.cut = -1 * log10(0.01)
  col = rep('lightgrey', times=length(p.val))
  c = which(dfGenes$adj.P.Val < p.adj.cut)
  col[c] = 'red'
  plot(fc, p.val, pch=20, xlab='Fold Change', ylab='-log10 P.Value', col=col, main=main, xlim=fc.lim)
  abline(v = 0, col='grey', lty=2)
  abline(h = y.cut, col='red', lty=2)
  # second cutoff for adjusted p-values
  y.cut = quantile(p.val[c], probs=0.95)
  abline(h = y.cut, col='red')
  # identify these genes
  g = which(p.val > y.cut)
  lab = dfGenes[g, 'SYMBOL']
  text(dfGenes$logFC[g], y = p.val[g], labels = lab, pos=2, cex=0.6)
}

plotMeanFC = function(m, dat, p.cut, title){
  col = rep('grey', length.out=(nrow(dat)))
  col[dat$adj.P.Val < p.cut] = 'red'
  #mh = cut(m, breaks=quantile(m, 0:50/50), include.lowest = T)
  plot(m, dat$logFC, col=col, pch=20, main=title, xlab='log Mean', ylab='logFC', ylim=c(-2, 2), cex=0.6)
}

## utility function for plotting
getms = function(f){
  m = mean(f)
  se = sd(f)
  m.up = m+1.96*se
  m.down = m-1.96*se
  ret= c(m, m.up, m.down)
  names(ret) = c('m', 'm.up', 'm.down')
  return(ret)
}


## shannon diversity
## code copied from
## https://rdrr.io/github/microbiome/microbiome/src/R/diversities.R
shannon <- function(x) {
  
  # Ignore zeroes
  x <- x[x > 0]
  
  # Species richness (number of species)
  S <- length(x)
  
  # Relative abundances
  p <- x/sum(x)
  
  # Shannon index
  (-sum(p * log(p)))
  
}


### binning the variable into categories
f_binData = function(x, m=0, m2=0.3, b=4, lab=c(0, 1, 2)){
  br = c(m, seq(m2, max(x), length.out = b-1))
  return(cut(x, breaks = br, right = T, include.lowest = T, labels = lab))
}


## hist2
## create an overlapping histogram by splitting data vector
## into 2 groups, typical usage would be to check overlap and
## balance of covariates between two experimental groups
## code idea from https://stackoverflow.com/questions/3541713/how-to-plot-two-histograms-together-in-r
hist2 = function(x, y, main='', xlab='Values', ylab='Density', legends=c('Value 1', 'Value 2'), legend.pos='topright'){
  p1 = hist(x,plot=FALSE)
  p2 = hist(y,plot=FALSE)
  ## calculate the range of the graph
  xlim = range(p1$breaks,p2$breaks)
  ylim = range(0,p1$density,
               p2$density)
  plot(p1,xlim = xlim, ylim = ylim,
       col = rgb(1,0,0,0.4),xlab = xlab, ylab=ylab,
       freq = FALSE, ## relative, not absolute frequency
       main = main)
  plot(p2,xlim = xlim, ylim = ylim,
       xaxt = 'n', yaxt = 'n', ## don't add axes
       col = rgb(0,0,1,0.4), add = TRUE,
       freq = FALSE)
  legend(legend.pos, legends,
         fill = rgb(1:0,0,0:1,0.4), bty = 'n',
         border = NA)
}
