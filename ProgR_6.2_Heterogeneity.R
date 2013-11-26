################################################################
################################################################
# Script: ProgR_6.2_Heterogeneity.R
# Author: Ingrid Lonnstedt
# Date:  30/08/2013
# R version: R 3.0.2
# Details: Our own absolute CNs and heterogeneity
################################################################
################################################################



################################################################
################################################################
#
# File paths, functions and date
#
################################################################
################################################################
setwd(paste(getwd(), '/RESPONSIFY', sep=''))
#setwd('/Users/lonnstedt/Documents/RESPONSIFY')

date = format(Sys.Date())
source(file.path(getwd(), 'programs', 'ProgR_0.1_Functions.R'))

################################################################
################################################################
#
# Check which HAPSEG solution looks best
#
################################################################
################################################################
info = read.delim(paste(prefix.ann, '../../Infokey.txt', sep = ''), stringsAsFactors = F)
info$sample.name = format(info$sample.name, nsmall = 4, trim = T)
library(PSCBS)
source(file.path(getwd(), 'programs', 'make_HAPSEG_seg_obj.R'))


for (i in 1:56){
  d = file.path(getwd(), 'ABSOLUTE', 'SampleData')
  d = paste(d,'/', info$sample.name[i], sep = '')
  dir.create(d, recursive = T, showWarnings = F)
  pdf(file.path(d, paste(info$sample.name[i], '_TAPS_after_HAPSEGS.pdf', sep = '')), 
      paper = 'a4r', height = 6, width = 12)  
  par(mfrow = c(1,3), mar = c(4, 4, 4, .1), oma = c(0, 0, 1, 0))
  for (h in 1:3){
    dh = paste('HAPSEG', h, sep = '')
    d = file.path(getwd(), 'ABSOLUTE', dh)
    d = paste(d,'/', info$sample.name[i], sep = '')
    #load(file.path(d, paste(info$sample.name[i], '_', info$arrname[i], '_segdat.RData', sep = '')))
    seg.dat.fn = file.path(d, paste(info$sample.name[i], '_', info$arrname[i], '_segdat.RData', sep = ''))
    seg.dat = HAPSEGMakeSegObj(seg.dat.fn, gender = 'Female', filter_segs=FALSE, min_probes=10, max_sd=100, 
                               verbose=TRUE)  
    as.seg.dat = seg.dat$as.seg.dat
    nr = as.integer(nrow(as.seg.dat))
    mat = as.seg.dat[1:(nr/2),]
    names(mat)[6] = 'a1'
    mat$a2 = as.seg.dat$copy_num[(nr/2+1):(nr)]
    mat$W = mat$W*2
    
    x = mat$a1 + mat$a2
    y = 2*(mat$a2/(x)-.5)
    xlim = c(quantile(x, probs = c(.025,.975), na.rm = T))
    ylim = c(quantile(y, probs = c(0, 1), na.rm = T))
    plot(x, y, pch = 16, col = col.transp('grey', .3), 
         xlim = xlim, ylim = ylim, xlab = 'Total CN intensity', 
         ylab = 'Allelic Imbalance 2*|BAF - 0.5|', main = dh)
    title(main = info$sample.name[i], outer = T, line = 0)
  }
  dev.off()
}

################################################################
################################################################
#
# Process SNP array data for CN and purity estimation
#
################################################################
################################################################
info = read.delim(file.path(getwd(), 'AROMA','annotationData','Infokey.txt'), stringsAsFactors = F)
info$sample.name = format(info$sample.name, nsmall = 4, trim = T)


mysamples = c(24, 30, 56, 36, 3, 21)

###Arguments
i = 24


####################
###Read PSCBS data
####################
library(PSCBS)
d = file.path(getwd(), 'AROMA', 'reports', 'SNP_segmentation', 'CBSfitsWithNormals')
da = loadObject(paste(d, '/fit', i, '.RData', sep = ''))
mat = da$output
names(mat)[names(mat) == 'c1Mean'] = 'a1'
names(mat)[names(mat) == 'c2Mean'] = 'a2'
names(mat)[names(mat) == 'tcnStart'] = 'Start'
names(mat)[names(mat) == 'tcnEnd'] = 'End'
names(mat)[names(mat) == 'chromosome'] = 'Chromosome'
mat$length = mat$End - mat$Start
mat = subset(mat, !is.na(length) & Chromosome<23 & !is.na(a2) & !is.na(a1))
mat$W = mat$length/sum(mat$length)

#Figure
x = mat$a1 + mat$a2
y = 2*(mat$a2/(x)-.5)
xlim = c(quantile(x, probs = c(0,.975), na.rm = T))
ylim = c(quantile(y, probs = c(0, 1), na.rm = T))
plot(x, y, pch = 16, col = col.transp('grey', .3), 
     xlim = xlim, ylim = ylim, xlab = 'Total CN intensity', 
     ylab = 'Allelic Imbalance 2*|BAF - 0.5|', main = 'TAPS')
abline(v = 1, lty = 'dotted')
d = file.path(getwd(), 'HET', info$sample.name[i])
dir.create(d, recursive = T, showWarnings = F)
tx = ifelse((info$ind[i] + 100) %in% info$ind & i !=54, 'Sample has germline SNP array',
            'Sample does not have germline SNP array')
mtext(tx)
dev.print(png, file=file.path(d, paste(info$sample.name[i],'_Standardized_TAPS.png', 
                                       sep = '')),  width=640, height=640)

####################
###Read HAPSEG output data
####################
source(file.path(getwd(), 'programs', 'make_HAPSEG_seg_obj.R'))
d = file.path(getwd(), 'ABSOLUTE', 'WithNormals','HAPSEG2')
d = paste(d,'/', info$sample.name[i], sep = '')
#load(file.path(d, paste(info$sample.name[i], '_', info$arrname[i], '_segdat.RData', sep = '')))
seg.dat.fn = file.path(d, paste(info$sample.name[i], '_', info$arrname[i], '_segdat.RData', sep = ''))
seg.dat = HAPSEGMakeSegObj(seg.dat.fn, gender = 'Female', filter_segs=FALSE, min_probes=10, max_sd=100, 
                           verbose=TRUE)  
as.seg.dat = seg.dat$as.seg.dat
nr = as.integer(nrow(as.seg.dat))
mat = as.seg.dat[1:(nr/2),]
names(mat)[6] = 'a1'
mat$a2 = as.seg.dat$copy_num[(nr/2+1):(nr)]
mat$W = mat$W*2


####################
###Starting values
####################
mystart = start.alphamax.f(mat)
alpha = mystart$alphamax
F = mystart$F
plot.transformed(mat, alpha, 
                 F, xlim = c(0,10), ylim = c(0,6))

####################
###Evaluate rho function and weights
####################

#Input
w = 1-exp(-mat$length/500000) #Sample 24, 30, 36, 3, 21
w = 1-exp(-mat$length/5000000) #Sample 56
rhofunc = 'tukey'
nmad = 3 #Sample 24, 30, 56, 36


#Change weight constant until 1'000'000 length segs have almost top weight
#Change nmod until Contrib to rho  gets constant at good Distance to lattice point
evaluate.weights(mat, mystart$alphamax, mystart$F, nmad = nmad, w=w)

evaluate.weights.2D(mat, alpha = mystart$alphamax, F = mystart$F, 
                    nmad = 3, w=w)
evaluate.weights.2D = function(mat, alpha, F, rhofunc = 'tukey', 
                               nmad = 1, w = 1, maxCN = 24, nu = 3){
  
  #Rotate first
  Fm1 = solve(F)
  c1 = x[,1]*Fm1[1,1] + x[,2] * Fm1[1,2]
  c2 = x[,1]*Fm1[2,1] + x[,2] * Fm1[2,2]
  
  #Assign observations to closest lattice point
  A2 = 0:maxCN
  A1 = matrix(A2, ncol = maxCN + 1, nrow = maxCN + 1)
  A2 = t(A1)
  E1 = c(alphamax*A1 + (1-alphamax))
  E2 = c(alphamax*A2 + (1-alphamax))
  alldist = apply(cbind(c1, c2), 1, dist, y = cbind(E1, E2))
  closest = apply(alldist, 2, which.min)
  
  #Calculate residuals
  e = cbind(E1, E2)[closest,]
  xcoord = mat$a2-e[,2]
  ycoord = mat$a1-e[,1]
  X = cbind(ycoord, xcoord)
  
  c = cov.trob(X, center = c(0,0), nu = nu)$cov
  D = log(mahalanobis(X, center = c(0,0), cov = c) )
  #Did not work without tukey (even with weight)
  rhoargs = list(c = nmad)
  r = D*w
  eachrho = do.call(rhofunc, args = list(x = r, args = rhoargs))
  
  #Plot
  par(mfrow=c(2,2), mar = c(4, 4, 3, .1))
  plot(r, eachrho, log = 'xy', xlab = 'r = Distance * Weight', 
       ylab = 'Contribution to rho')
  plot(D, eachrho, log = 'xy', xlab = 'Distance to lattice point', 
       ylab = 'Contribution to rho')
  plot(mat$length, w, log = 'xy', xlab = 'Segment length', 
       ylab = 'Weight')
  plot(D, r, log = 'xy', ylab = 'r = Distance * Weight', 
       xlab = 'Distance to lattice point')
  title(main = 'Evaluation of rho function and weights', line = -1,
        outer = T)
}


####################
###Optimization
####################

#w, rhofunc and nmad from above!
di = getdist(mat, mystart$alphamax, mystart$F) 
rhoargs = list(k = mad(di), c = mad(di)*nmad, nu = 3, maxCN = 24)
rhoargs = list(k = mad(di), c = nmad, nu = 3, maxCN = 24)
x = cbind(mat$a1, mat$a2)

#pars = solnp(pars = c(mystart$alphamax, mystart$F), fun = rho.alphamax.f, x = x,
#             w = w, maxCN = 24, rhofunc = rhofunc, 
#             rhoargs = rhoargs, UB = c(1, mystart$UBf), 
#             LB = c(0,mystart$LBf))$pars
pars = solnp(pars = c(mystart$alphamax, mystart$F), fun = rho.alphamax.f,
             x = x,
             w = w, maxCN = 24, rhofunc = rhofunc, 
             rhoargs = rhoargs)$pars
pars = solnp(pars = c(mystart$alphamax, mystart$F), fun = rho.2Dt, y = x,
             w = w, ineqfun = rho.2Dt.ineq, ineqLB = 0, ineqUB = Inf,
             rhoargs = rhoargs, LB = c(0, 0, -Inf, -Inf, 0),
             UB = c(1, Inf, Inf, Inf, Inf))$pars
alpha = pars[1]
F = matrix(pars[-1], ncol = 2, nrow = 2)

plot.transformed.2D(mat, alpha, F, w=w, rhoargs=rhoargs, 
                    xlim = c(0, 10), ylim = c(0,6))
plot.transformed(mat, alpha, F, w=w, rhoargs=rhoargs, rhofunc = rhofunc, 
                 xlim = c(0,10), ylim = c(0,6))

#All possible alphas and Fs
maxCN = 24
alpha = pars[1]
F = matrix(pars[-1], ncol = 2, nrow = 2)
EEs = alpha*c(0:maxCN) + (1-alpha) 
alphas = NULL
fs = NULL
for (j in 1:length(EEs)){
  alphas = c(alphas, alpha/(EEs[j]+alpha))
  if (j == 1) fs = 1 else {
    fs = c(fs, fs[j-1]*alphas[j-1]/alphas[j])
  }
}

###For subclonal CN estimation
#Estimate cov without using weights to get maximum cov
#Estimate each Mahalanobis distance with cov and incorporated weight in some way

###For rotation with t
#Count each distance weighted down


### Try with 2D: did not move the start positions
rho.2Dt.ineq = function(pars, y, w, rhoargs = list(nu = 3, maxCN = 24)){
  alphamax = pars[1]
  F = matrix(pars[-1], ncol = 2, nrow = 2)
  
  #Rotate first (observed points' scale did not work)
  Fm1 = solve(F)
  c1 = x[,1]*Fm1[1,1] + x[,2] * Fm1[1,2]
  c2 = x[,1]*Fm1[2,1] + x[,2] * Fm1[2,2]
  
  #Assign observations to closest lattice point
  A2 = 0:maxCN
  A1 = matrix(A2, ncol = maxCN + 1, nrow = maxCN + 1)
  A2 = t(A1)
  E1 = c(alphamax*A1 + (1-alphamax))
  E2 = c(alphamax*A2 + (1-alphamax))
  alldist = apply(cbind(c1, c2), 1, dist, y = cbind(E1, E2))
  closest = apply(alldist, 2, which.min)
  
  #Calculate residuals
  e = cbind(E1, E2)[closest,]
  xcoord = y[,2]-e[,2]
  ycoord = y[,1]-e[,1]
  X = cbind(ycoord, xcoord)
  
  c = cov.trob(X, center = c(0,0), nu = rhoargs$nu)$cov
  det(c)
}

minfunc = function(pars, y, w=1, 
               rhoargs = list(nu = 3, nmad = 3,
               rhofunc = 'tukey'), 
               minfunc = minfunc.1D){
  individual.rho = do.call(minfunc, list(y=y, w=w, rhoargs = rhoargs,
          alpha = pars[1], F = matrix(pars[-1], ncol = 2, nrow = 2)))
  sum(individual.rho)
}

minfunc.2D = function(y, w=1, alpha, F, rhoargs = NULL){
  
}

rho.2Dt = function(pars, y, w, rhoargs = list(nu = 3, maxCN = 24, nmad = 3)){
  alphamax = pars[1]
  F = matrix(pars[-1], ncol = 2, nrow = 2)
    
  #Rotate first (observed points' scale did not work)
  Fm1 = solve(F)
  c1 = x[,1]*Fm1[1,1] + x[,2] * Fm1[1,2]
  c2 = x[,1]*Fm1[2,1] + x[,2] * Fm1[2,2]
  
  #Assign observations to closest lattice point
  A2 = 0:maxCN
  A1 = matrix(A2, ncol = maxCN + 1, nrow = maxCN + 1)
  A2 = t(A1)
  E1 = c(alphamax*A1 + (1-alphamax))
  E2 = c(alphamax*A2 + (1-alphamax))
  alldist = apply(cbind(c1, c2), 1, dist, y = cbind(E1, E2))
  closest = apply(alldist, 2, which.min)
  
  #Calculate residuals
  e = cbind(E1, E2)[closest,]
  xcoord = y[,2]-e[,2]
  ycoord = y[,1]-e[,1]
  X = cbind(ycoord, xcoord)
  
  c = cov.trob(X, center = c(0,0), nu = rhoargs$nu)$cov
  D = sqrt(mahalanobis(X, center = c(0,0), cov = c))
      #Did not work without tukey (even with weight)
  rhoargs = list(k = nmad*mad(di), c = nmad*mad(di))
  rhoargs = list(k = nmad*mad(di), c = nmad)
  r = D*w
  sum(do.call(rhofunc, args = list(x = r, args = rhoargs)))
}

di = dist(x, e)
xcoord = x[,2]-e[,2]
ycoord = x[,1]-e[,1]
wei = w/sum(w)
scaleup = 1/min(w)
wei*scaleup


rho.alphamax.f2 = function(pars, x, w, maxCN=24, rhofunc = 'huber', 
                          rhoargs = NULL){
#This does not work in practice on 24
  alphamax = pars[1]
  F = matrix(pars[-1], ncol = 2, nrow = 2)
  
  #Rotate first
  Fm1 = solve(F)
  c1 = x[,1]*Fm1[1,1] + x[,2] * Fm1[1,2]
  c2 = x[,1]*Fm1[2,1] + x[,2] * Fm1[2,2]
  
  #Assign observations to closest lattice point
  A2 = 0:maxCN
  A1 = matrix(A2, ncol = maxCN + 1, nrow = maxCN + 1)
  A2 = t(A1)
  E1 = c(alphamax*A1 + (1-alphamax))
  E2 = c(alphamax*A2 + (1-alphamax))
  alldist = apply(cbind(c1, c2), 1, dist, y = cbind(E1, E2))
  closest = apply(alldist, 2, which.min)
  
  #Calculate residuals
  e = cbind(E1, E2)[closest,]
  r = dist(c(c1, c2), e)*w
  sum(do.call(rhofunc, args = list(x = r, args = rhoargs)))
}


####################
### Choise of alpha with help from mutations
####################

maf.fn = paste(getwd(), 
               '/AROMA/rawData/responsify/Exome/MAF/Data3AllVariants_mutect_nonsilent_mindepth20_withAlleleFreqs', '_', 
               info$sample.name[i], '.maf', sep = '')
maf.fn = paste(getwd(), 
               '/Datasets/Exome/MAF/Data3AllVariants_mutect_nonsilent_mindepth20_withAlleleFreqs', '_', 
               info$sample.name[i], '.maf', sep = '')

maf = read.delim(maf.fn, row.names = NULL, stringsAsFactors = FALSE, 
                 check.names = FALSE, na.strings = c("NA", "---"), 
                 blank.lines.skip = TRUE, comment.char = "#")
maf$depth = maf$t_alt_count + maf$t_ref_count

#Assign mutations to segments
maf$segment = NA
for (seg in 1:nrow(mat)){
  ix = which(maf$Start_position >= mat$Start.bp[seg] & 
               maf$Start_position <= mat$End.bp[seg]
             & as.character(maf$Chromosome) == as.character(mat$Chromosome[seg]))
  maf$segment[ix] = seg
} 

muts = unique(subset(maf, select = c('Chromosome','Start_position','End_position',
                                     't_alt_count','t_ref_count','segment')))
muts$af = muts$t_alt_count/(muts$t_alt_count+muts$t_ref_count)
muts = subset(muts, !is.na(segment) )

#Plot observed variant frequencies 
#NB: investigate n = 1:4 if surely gridrow 4 must correspond to at least 0 copies
par(mfrow = c(2,2))
for (n in 1:4){
  F = matrix(pars[-1], ncol = 2, nrow = 2)
  alpha = alphas[n]
  F = F*fs[n]
  Fm1 = solve(F)
  c1 = mat$a1*Fm1[1,1] + mat$a2 * Fm1[1,2]
  c2 = mat$a1*Fm1[2,1] + mat$a2 * Fm1[2,2]
  c1 = (c1 - (1 - alpha))/alpha
  c2 = (c2 - (1 - alpha))/alpha
  
  
  minor = alpha*c1/(alpha*(c1+c2) + (1-alpha)*2)
  major = alpha*c2/(alpha*(c1+c2) + (1-alpha)*2)
  both = alpha*(c1+c2)/(alpha*(c1+c2) + (1-alpha)*2)
    
  #plot(both, ylim = c(0,1))
  #points(major, col = 'blue')
  #points(minor, col = 'red')
  #points(muts$segment, muts$af, col = 'green')
  
  ord = order(both, major)
  plot(both[ord], ylim = c(0,1), ylab = 'Allele fraction', xlab = 'CN segment', 
       main = paste('CN = 0 at gridline',n, 'from bottom/left'))
  points(major[ord], col = 'blue')
  points(minor[ord], col = 'red')
  muts$index = NA
  for (j in 1:nrow(muts)){
    if (!is.na(muts$segment[j])) muts$index[j] = 
      which(ord == muts$segment[j])
  }
  points(muts$index, muts$af, col = 'green', pch = 16)
  
}

#Choose n
n = 3 #Sample 24
n = 2 #Sample 21
n = 6 #Sample 30
n = 9 #Sample 3
n = 7 #Sample 56
alpha = alphas[n]
F = matrix(pars[-1], ncol = 2, nrow = 2)
plot.transformed(mat, alpha, F = F*fs[n], w=w, rhoargs=rhoargs, rhofunc = rhofunc, 
                 xlim = c(0,4), ylim = c(0,2.5))

#Set transformed CNs
Fm1 = solve(F*fs[n])
mat$a1t = mat$a1*Fm1[1,1] + mat$a2 * Fm1[1,2]
mat$a2t = mat$a1*Fm1[2,1] + mat$a2 * Fm1[2,2]

plot.lattice(x = mat$a2t, y = mat$a1t, alpha = alphas[n],
             xlim = c(0,2.5), ylim = c(0,2), 
             main = info$sample.name[i])
d = file.path(getwd(), 'HET', info$sample.name[i])
dev.print(png, file=file.path(d, paste(info$sample.name[i],'_latticeplot_transformed.png', 
                                       sep = '')),  width=740, height=600)



####################
### Nonlinearity fix
####################
#Not possible on samples 56, 36, 3, 21
alpha = alphas[n]
dummy = nonlinAdj.asympOrig(mat$a1t, alpha = alpha, maxCN = 6)
dummy = nonlinAdj.asympOrig(mat$a2t, alpha = alpha, maxCN = 6)
dummy = nonlinAdj.asympOrig(c(mat$a1t, mat$a2t), alpha = alpha, maxCN = 6)
dummy = nonlinAdj.logis(mat$a1t, alpha = alpha, maxCN = 6)
dummy = nonlinAdj.logis(mat$a2t, alpha = alpha, maxCN = 6)#24
dummy = nonlinAdj.logis(c(mat$a1t, mat$a2t), alpha = alpha, maxCN = 6)

mat$a2ta = nonlinAdj.logis(y=mat$a2t, alpha = alpha, maxCN = 6)
mat$a1ta = mat$a1t

plot.lattice(x = mat$a2ta, y = mat$a1ta, alpha = alphas[n],
             xlim = c(0,6), ylim = c(0,2.5))
d = file.path(getwd(), 'HET', info$sample.name[i])
dev.print(png, file=file.path(d, paste(info$sample.name[i],'_latticeplot_transformed_reatten.png', 
                                       sep = '')),  width=740, height=600)



####################
### Assign decimal copy numbers 
####################

#Calculate lattice points from alpha
alpha = alphas[n]
mat$c1 = (mat$a1t - (1 - alpha))/alpha
mat$c2 = (mat$a2t - (1 - alpha))/alpha
plot.lattice(x = mat$c2, y = mat$c1, alpha = 1,
             xlim = c(0,7), ylim = c(0,5), 
             main = info$sample.name[i])


####################
### Overlay lattice points: subclonal CN segments
####################

#Calculate lattice points from alpha
alpha = alphas[n]
maxCN = 24
A2 = 0:maxCN
A1 = matrix(A2, ncol = maxCN + 1, nrow = maxCN + 1)
A2 = t(A1)
E1 = c(alpha*A1 + (1-alpha))
E2 = c(alpha*A2 + (1-alpha))


### OVERLAY PLOT Without reattenuation
par(mfrow = c(1,2))
### Each segment's distance to it's lattice point
x = cbind(mat$a1t, mat$a2t)
alldist = apply(x, 1, dist, y = cbind(E1, E2))
closest = apply(alldist, 2, which.min)
e = cbind(E1, E2)[closest,]
di = dist(x, e)
xcoord = x[,2]-e[,2]
ycoord = x[,1]-e[,1]


mod = rlm(ycoord~xcoord)
smoothScatter(xcoord, ycoord, xlim = c(-alpha, alpha), 
              ylim = c(-alpha, alpha), main = 'Without reattenuation')
abline(v=0, h=0)
rect(xleft = -alpha/2, ybottom = -alpha/2, xright = alpha/2, 
     ytop = alpha/2, lty = 'dotted')

### OVERLAY PLOT With reattenuation
### Each segment's distance to it's lattice point
x = cbind(mat$a1ta, mat$a2ta)
alldist = apply(x, 1, dist, y = cbind(E1, E2))
closest = apply(alldist, 2, which.min)
e = cbind(E1, E2)[closest,]
di = dist(x, e)
xcoord = x[,2]-e[,2]
ycoord = x[,1]-e[,1]

smoothScatter(xcoord, ycoord, xlim = c(-alpha, alpha),
              ylim = c(-alpha, alpha), main = 'With reattenuation')
abline(v=0, h=0)
rect(xleft = -alpha/2, ybottom = -alpha/2, xright = alpha/2, 
     ytop = alpha/2, lty = 'dotted')

d = file.path(getwd(), 'HET', info$sample.name[i])
dev.print(png, file=file.path(d, 
          paste(info$sample.name[i],'_overlay_transformed_reatten.png', 
                                       sep = '')),  width=640, height=400)


### Distributional distiances
x = cbind(mat$a1t, mat$a2t)
alldist = apply(x, 1, dist, y = cbind(E1, E2))
closest = apply(alldist, 2, which.min)
e = cbind(E1, E2)[closest,]
di = dist(x, e)
xcoord = x[,2]-e[,2]
ycoord = x[,1]-e[,1]
X = cbind(ycoord, xcoord)
wei = w/sum(w)
#X = cbind(ycoord, xcoord)[1:10,]
#wei = w[1:10]/sum(w[1:10])
scaleup = 1/min(w)
wei*scaleup

#Multivariate t robust covariance
c = cov(X)
c5 = cov.trob(X, center = c(0,0), nu = 5)$cov
c4 = cov.trob(X, center = c(0,0), nu = 4)$cov
c3 = cov.trob(X, center = c(0,0), nu = 3)$cov #more robust
cn = cov.trob(X, center = c(0,0), nu = nrow(X)-2)$cov #less robust

#Or with weights
c = cov(X)
c5 = cov.trob(X, center = c(0,0), nu = 5, 
              wt = round(w*scaleup))$cov
c4 = cov.trob(X, center = c(0,0), nu = 4, wt = wei)$cov
c3 = cov.trob(X, center = c(0,0), nu = 3, wt = wei)$cov #more robust
cn = cov.trob(X, center = c(0,0), nu = nrow(X)-2, wt = wei)$cov #less robust

D = mahalanobis(X, center = c(0,0), cov = c)
D5 = mahalanobis(X, center = c(0,0), cov = c5)
D4 = mahalanobis(X, center = c(0,0), cov = c4)
D3 = mahalanobis(X, center = c(0,0), cov = c3)
Dn = mahalanobis(X, center = c(0,0), cov = cn)

par(mfrow = c(2,2))
plot(density(D))
plot(density(D5))
plot(density(D4))
plot(density(D3))
plot(density(Dn))

par(mfrow = c(2,2))
hist(-2*log(D))
hist(-2*log(Dn))
hist(-2*log(D5))
#hist(-2*log(D4))
hist(-2*log(D3))

par(mfrow = c(2,2))
qqplot(qexp(ppoints(length(D)), rate = .5), D,
       main = 'Normal', ylim = NULL)
abline(0, 1, col = 'gray')
qqplot(qexp(ppoints(length(D)), rate = .5), Dn,
       main = 't_n', ylim = NULL)
abline(0, 1, col = 'gray')
qqplot(qexp(ppoints(length(D)), rate = .5), D5,
       main = 't_5', ylim = NULL)
abline(0, 1, col = 'gray')
qqplot(qexp(ppoints(length(D)), rate = .5), D3,
       main = 't_3', ylim = NULL)
abline(0, 1, col = 'gray')



pal = brewer.pal(9, 'YlOrRd')
cols = pal[-c(1:2)]
cols = colorRampPalette(cols)(100)

par(mfrow = c(1,2))
colby = log(Dn)
colbreaks = seq(min(colby), max(colby), length = length(cols))
delta = colbreaks[2] - colbreaks[1]
colind = (colby-colbreaks[1]) %/% delta + 1
colind[colind>length(cols)] = length(cols)
mycols = col.transp(cols[colind], .5)
plot(xcoord, ycoord, xlim = c(-alpha, alpha), col = mycols, pch = 16,
              ylim = c(-alpha, alpha), main = 'Colour by D^2_n')
abline(v=0, h=0)
rect(xleft = -alpha/2, ybottom = -alpha/2, xright = alpha/2, 
     ytop = alpha/2, lty = 'dotted')
index = Dn>2
points(xcoord[index], ycoord[index], pch = 16)

colby = log(D3)
colbreaks = seq(min(colby), max(colby), length = length(cols))
delta = colbreaks[2] - colbreaks[1]
colind = (colby-colbreaks[1]) %/% delta + 1
colind[colind>length(cols)] = length(cols)
mycols = col.transp(cols[colind], .5)
plot(xcoord, ycoord, xlim = c(-alpha, alpha), col = mycols, pch = 16,
     ylim = c(-alpha, alpha), main = 'Colour by D^2_3')
abline(v=0, h=0)
rect(xleft = -alpha/2, ybottom = -alpha/2, xright = alpha/2, 
     ytop = alpha/2, lty = 'dotted')
index = D3>3
points(xcoord[index], ycoord[index], pch = 16)

##################
#Junk
hist(di, xlim = c(0,1), breaks = 1000)
hist(exp(-di**2), xlim = c(0,1), breaks = 1000)
plot.lattice(x = mat$a2t, y = mat$a1t, alpha = alphas[n],
             xlim = c(0,2.5), ylim = c(0,2), 
             main = info$sample.name[i], colby = log(D3))
d = file.path(getwd(), 'HET', info$sample.name[i])
dev.print(png, file=file.path(d, 
  paste(info$sample.name[i],'_latticeplot_transformed_Dcol.png', 
                sep = '')),  width=740, height=600)


qqplot(qchisq(ppoints(length(D)), df = 3), D,
       main = expression("Q-Q plot of Mahalanobis" * ~D^2 *
                           " vs. quantiles of" * ~ chi[3]^2))
abline(0, 1, col = 'gray')


#with mu
mu = matrix(wei %*% X, nrow = 2, ncol = nrow(X))
dv = t(X) - mu
sig = matrix(0, nrow = 2, ncol = 2)
for (j in 1:length(wei)){
  sig = sig + (wei[j]*(dv[,j] %*% t(dv[,j])))
}
sig = sig/(1-sum(wei**2))

#without mu
dv = t(X)
sig = matrix(0, nrow = 2, ncol = 2)
for (j in 1:length(wei)){
  sig = sig + (wei[j]*(dv[,j] %*% t(dv[,j])))
}
sig = sig/(1-sum(wei**2))

#Just ordinary covariance
R = cov2cor(cov(X))
Rm1 = solve(R)
t = numeric(nrow(mat))
apply(X, 2, sd)
x = cbind(X[])
for (j in 1:nrow(mat)){
  xj = X[j,]
  t[j] = t(X[j,]) %*% Rm1 %*% X[j,]
}


####################
### Output HAPSEG transformed version dataset
####################

#Without reattenuation output for sample 24, 56, 21
seg.dat.new = hapseg.format(seg.dat, mat, xname = 'a2t', 
                            yname = 'a1t')

#Without reattenuation output for sample 30
seg.dat.new = hapseg.format(seg.dat, mat, xname = 'a2ta', 
                            yname = 'a1ta')

#Save workspace
d = file.path(getwd(), 'HET', info$sample.name[i])
save.image(file = file.path(d, paste(info$sample.name[i], '.RData', sep = '')))


#Print the new object
d = file.path(getwd(), 'ABSOLUTE', 'WithNormals', 'HAPSEG2')
d = paste(d,'/', info$sample.name[i], sep = '')
seg.dat = seg.dat.new
save(seg.dat, file = file.path(d, paste(info$sample.name[i], '_', info$arrname[i], 
                                              '_scaled_segdat.RData', sep = '')))


#Load
d = file.path(getwd(), 'HET', info$sample.name[i])
load(file = file.path(d, paste(info$sample.name[i], '.RData', sep = '')))



#######################
### Subclonal mutations
#######################

n = 3
alpha = alphas[n]
c1 = pmax(mat$c1, 0)
c2 = pmax(mat$c2, 0)

minor = alpha*c1/(alpha*(c1+c2) + (1-alpha)*2)
major = alpha*c2/(alpha*(c1+c2) + (1-alpha)*2)
both = alpha*(c1+c2)/(alpha*(c1+c2) + (1-alpha)*2)


#Plot variant fractions
par(mfrow = c(1,1))
ord = order(both, major)
plot(both[ord], ylim = c(0,1), ylab = 'Variant fraction', xlab = 'CN segment')
points(major[ord], col = 'blue')
points(minor[ord], col = 'red')
muts$index = NA
for (j in 1:nrow(muts)){
  if (!is.na(muts$segment[j])) muts$index[j] = 
    which(ord == muts$segment[j])
}
points(muts$index, muts$af, col = 'green', pch = 16)

#Highlight 'subclonal' mutations
sign.level = .01
muts$major = muts$t_alt_count >= qbinom(sign.level/2, size = muts$t_alt_count + muts$t_ref_count,
                                        prob = major[muts$segment]) &
  muts$t_alt_count <= qbinom(1 - sign.level/2, size = muts$t_alt_count + muts$t_ref_count,
                             prob = major[muts$segment])
muts$minor = muts$t_alt_count >= qbinom(sign.level/2, size = muts$t_alt_count + muts$t_ref_count,
                                        prob = minor[muts$segment]) &
  muts$t_alt_count <= qbinom(1 - sign.level/2, size = muts$t_alt_count + muts$t_ref_count,
                             prob = minor[muts$segment])
muts$both = muts$t_alt_count >= qbinom(sign.level/2, size = muts$t_alt_count + muts$t_ref_count,
                                       prob = both[muts$segment]) &
  muts$t_alt_count <= qbinom(1 - sign.level/2, size = muts$t_alt_count + muts$t_ref_count,
                             prob = both[muts$segment])
muts$subclonal = !muts$both & !muts$major & !muts$minor


points(muts$index[muts$subclonal], muts$af[muts$subclonal], col = 'darkgreen', pch = 16)

#Junk
par(mfrow = c(3,1))
hist(na.omit(muts$af-both[muts$segment]), breaks = 100, main = '',
     xlab = 'Variant fraction - E(fraction if in both)')
hist(na.omit(muts$af-major[muts$segment]), breaks = 100,  main = '',
     xlab = 'Variant fraction - E(fraction if in major)')
hist(na.omit(muts$af-minor[muts$segment]), breaks = 100,  main = '',
     xlab = 'Variant fraction - E(fraction if in minor)')
#hist(muts$af, breaks = 50)

#Mutations on top of CN grid
plot.het(x = mat$a2t, y = mat$a1t, alpha = alphas[n],
         xlim = c(0,3), ylim = c(0,2.5), xlab = 'Major CN', ylab = 'Minor CN')
points(mat$a2t[maf$segment], mat$a1t[maf$segment], pch = '.', cex = 2)

points(mat$a2t[(muts$segment)[submuts]], mat$a1t[(muts$segment)[submuts]], 
       pch = 4, cex = 1)





#Plot
plot.het = function(x, y, alpha, main = NULL, xlim = NULL, xlab = NULL, ylab = NULL,
                        ylim = NULL, maxCN = 24){
  require('RColorBrewer')
  require('graphics')
  
  #Calculate lattice points from alpha
  A2 = 0:maxCN
  A1 = matrix(A2, ncol = maxCN + 1, nrow = maxCN + 1)
  A2 = t(A1)
  
  E1 = c(alpha*A1 + (1-alpha))
  E2 = c(alpha*A2 + (1-alpha))
  
  ### Each segment's distance to it's lattice point -> weight
  alldist = apply(cbind(y,x), 1, dist, y = cbind(c(E1), c(E2)))
  closest = apply(alldist, 2, which.min)
  e = cbind(c(E1), c(E2))[closest,]
  r = dist(cbind(y,x), e)
  
  #Colour palette
  pal = brewer.pal(9, 'YlOrRd')
  cols = pal[-c(1:2)]
  cols = colorRampPalette(cols)(100)
    
  #Plot frame
  if (is.null(xlim)) xlim = c(0, max(c(x, E2), na.rm = T))
  if (is.null(ylim)) ylim = c(0, max(y, na.rm = T)*1.2)
  plot(x, y, pch = 16, col = col.transp('grey', .3), type = 'n',
       xlim = xlim, ylim = ylim, xaxt = 'n', yaxt = 'n', xlab = xlab, ylab = ylab)
  
  #Points
  colby = r
  colbreaks = seq(min(colby), max(colby[x<xlim[2]]), length = length(cols))
  delta = colbreaks[2] - colbreaks[1]
  colind = (colby-colbreaks[1]) %/% delta + 1
  colind[colind>length(cols)] = length(cols)
  mycols = col.transp(cols[colind], .5)
  points(x, y, pch = 16, col = mycols, cex = 1)
  
  #Colour scale
  tmp = c(xlim[2]/2, xlim[2])
  tmp = seq(tmp[1], tmp[2], length=length(cols))
  points(tmp, rep(0, length(tmp)), col = cols, pch = 16)
  
  
  #Gridlines
  EEs = unique(E1)
  for (j in 1:length(EEs)){
    lines(range(EEs), rep(EEs[j], 2), lty = 'dotted', col = 'grey20')
    lines(rep(EEs[j], 2), range(EEs), lty = 'dotted', col = 'grey20')
  }
  abline(0,1)
  
  title(main = main, outer = T, line = -1)
}

#######################
### Refined plots Subclonal mutations
#######################

n = 3
alpha = alphas[n]


#Express vaf as number of copies per cell (cellular multiplicity)
muts$mult = muts$af/alpha*
  (alpha*((mat$c1 + mat$c2)[muts$segment])+(1-alpha)*2)

#Plot Copy numbers and multiplicities
plotmat = mat
ylim = c(-1, 10)
plotmat$c1 = pmin(plotmat$c1, ylim[2])
plotmat$c2 = pmin(plotmat$c2, ylim[2])
plotmat$tot = plotmat$c1 + plotmat$c2
plotmat$tot = pmin(plotmat$tot, ylim[2])
chroms = as.character(unique(plotmat$Chromosome))
par(mfrow = c(1,1), mar = c(4, 4, 1, .1))
for (j in 1:length(chroms)){
  chrom = chroms[j]
  xlim = c(0, max(c(plotmat$End.bp[as.character(plotmat$Chromosome) == chrom], 
            muts$End_position[as.character(muts$Chromosome) == chrom])))
  plot(1, 1, type = 'n', ylim = ylim, 
       xlim = xlim, ylab = 'Copies per cell', 
       xlab = paste('Position along chromosome', chrom))
  for (k in 0:round(ylim[2])){
    abline(h = k, lty = 'dotted', col = 'grey20')
  }
  
  ix.plotmat = which(as.character(plotmat$Chromosome) == chrom)
  for (k in ix.plotmat){
    lines(c(plotmat$Start.bp[k], plotmat$End.bp[k]), rep(plotmat$c1[k], 2), 
          lty  = 'dotted', lwd = 3, col = col.transp('blue'))
    lines(c(plotmat$Start.bp[k], plotmat$End.bp[k]), rep(plotmat$c2[k], 2), 
          lty  = 'dotted', lwd = 3, col = col.transp('red'))
    lines(c(plotmat$Start.bp[k], plotmat$End.bp[k]), rep(plotmat$c1[k]+plotmat$c2[k], 2), 
          lty  = 'solid', lwd = 3, col = col.transp('darkgreen'))
  }
  ix = which(as.character(muts$Chromosome) == chrom)
  points(muts$Start_position[ix], pmin(muts$mult[ix], ylim[2]), pch = 2)
  
}





####################
### WEHI retreat talk plot
####################
#I did this for sample 24

maf.fn = paste(getwd(), 
               '/AROMA/rawData/responsify/Exome/MAF/Data3AllVariants_mutect_nonsilent_mindepth20_withAlleleFreqs', '_', 
               info$sample.name[i], '.maf', sep = '')

maf = read.delim(maf.fn, row.names = NULL, stringsAsFactors = FALSE, 
                 check.names = FALSE, na.strings = c("NA", "---"), 
                 blank.lines.skip = TRUE, comment.char = "#")



#Assign mutations to segments
maf$segment = NA
for (seg in 1:nrow(mat)){
  ix = which(maf$Start_position >= mat$Start.bp[seg] & 
               maf$Start_position <= mat$End.bp[seg]
             & as.character(maf$Chromosome) == as.character(mat$Chromosome[seg]))
  maf$segment[ix] = seg
} 


myplot.transformed(mat, alpha = alphas[n], F = F, w=w, rhoargs=rhoargs, 
                 rhofunc = rhofunc, 
                 xlim = c(0,4), ylim = c(0,2.5))

myplot.transformed = function(mat, alpha, F, w=NULL, rhofunc = 'leastsq',
                            rhoargs = NULL, main = NULL, xlim = NULL, ylim = NULL,
                            maxCN = 24, color.by = 'eachrho'){
  require('RColorBrewer')
  require('graphics')
  
  #Recalculate transformed lattice points from F and alpha
  Fm1 = solve(F)  
  A2 = 0:maxCN
  A1 = matrix(A2, ncol = maxCN + 1, nrow = maxCN + 1)
  A2 = t(A1)
  
  E1 = c(alpha*A1 + (1-alpha))
  E2 = c(alpha*A2 + (1-alpha))
  
  ### Positions of lattice points in observed scale and
  ### each segment's distance to it's lattice point -> weight
  x = cbind(mat$a1, mat$a2)
  e1 = F[1,1] * E1 + F[1,2]*E2
  e2 = F[2,1] * E1 + F[2,2]*E2
  alldist = apply(x, 1, dist, y = cbind(e1, e2))
  closest = apply(alldist, 2, which.min)
  e = cbind(e1, e2)[closest,]
  if (is.null(w)) w = 1
  r = dist(x, e)/w
  eachweight = do.call(rhofunc, args = list(x = r, args = rhoargs))
  
  #Colour palette
  #display.brewer.all()
  pal = brewer.pal(9, 'YlOrRd')
  #display.brewer.pal(8, 'Set1')
  cols = pal[-c(1:2)]
  cols = colorRampPalette(cols)(100)
  
  #Transformed positions of data
  a1p = x[,1]*Fm1[1,1] + x[,2] * Fm1[1,2]
  a2p = x[,1]*Fm1[2,1] + x[,2] * Fm1[2,2]
  
  try(close.screen(all.screens = T), silent = T)
  split.screen(rbind(c(0.2, 0.98, 0.2, 0.9),
                     c(0.1,0.2,0.2, 0.9), 
                     c(0.2, 0.98, 0.1, 0.2),
                     c(0.1, 0.2, 0.1, 0.2)))
  screen(1)
  par(mar = c(0,0,0,0))
  
  
  #Plot frame
  if (is.null(xlim)) xlim = c(0, max(c(a2p, E2), na.rm = T))
  if (is.null(ylim)) ylim = c(0, max(a1p, na.rm = T)*1.2)
  plot(a2p, a1p, pch = 16, col = col.transp('grey', .3), 
       type = 'n', xlim = xlim, ylim = ylim, xaxt = 'n', yaxt = 'n')
  
  #Points
  colby = eachweight
  if (color.by[1] == 'log') colby = log(eachweight+1)
  if (is.numeric(color.by)) colby = color.by
  colbreaks = seq(min(colby), max(colby[a2p<xlim[2]]), length = length(cols))
  delta = colbreaks[2] - colbreaks[1]
  colind = (colby-colbreaks[1]) %/% delta + 1
  colind[colind>length(cols)] = length(cols)
  mycols = col.transp(cols[colind], .5)
  if (color.by == 'nocol') mycols = col.transp('grey', .3)
  points(a2p, a1p, pch = 16, col = mycols, cex = 1)
  
  #Gridlines
  EEs = unique(E1)
  for (j in 1:length(EEs)){
    lines(range(EEs), rep(EEs[j], 2), lty = 'dotted', col = 'grey20')
    lines(rep(EEs[j], 2), range(EEs), lty = 'dotted', col = 'grey20')
  }
  abline(0,1)
  
  #Colour scale
  if (color.by != 'nocol'){
    tmp = c(xlim[2]/2, xlim[2])
    tmp = seq(tmp[1], tmp[2], length=length(cols))
    points(tmp, rep(0, length(tmp)), col = cols, pch = 16)
  }
  
  #Mutations
  points(a2p[maf$segment], a1p[maf$segment], pch = '.', cex = 2)
#  Highlight mutations in subclonal segments
#  mat$outlier = 0
#  mat$outlier[colby > max(colby)*3/4] = 1
#  mat$outlier[maf$segment]==1
#  points((a2p[maf$segment])[mat$outlier[maf$segment]==1], 
#         (a1p[maf$segment])[mat$outlier[maf$segment]==1], pch = 4)
  
  
  #Single scale density y
  screen(2)
  par(mar = c(0,0,0,0))
  yden = density(a1p,kernel="gaussian", bw = .03)
  roughy = round(sum(diff(yden$y)**2)/diff(yden$x[1:2]), 3)
  yden$y = log(yden$y+1)
  plot(1, 1, ylim = ylim, xlim = range(yden$y), type = 'n', xaxt = 'n')
  polygon(yden$y, yden$x, col = col.transp(1, .5), border = NA)
  title(ylab = 'Transformed a1', line = -1, outer = T)
  
  #Single scale density x
  screen(3)
  par(mar = c(0,0,0,0))
  xden = density(a2p,kernel="gaussian", bw = .03)
  roughx = round(sum(diff(xden$y)**2)/diff(xden$x[1:2]), 3)
  xden$y = log(xden$y+1)
  plot(1, 1, xlim = xlim, ylim = range(xden$y), type = 'n', yaxt = 'n')
  polygon(xden$x, xden$y, col = col.transp(1, .5), border = NA)
  title(xlab = 'Transformed a2', line = -1, outer = T)
  
  screen(4)
  par(mar = c(0,0,0,0))
  plot(1,1, type = 'n', xaxt = 'n', yaxt = 'n', xlim = c(0,1), ylim = c(0,1), bty = 'n')
  text(0, 1, roughx, adj = c(0,1))
  text(1, 0, roughy, adj = c(1,0))
  
  title(main = main, outer = T, line = -1)
  close.screen(all.screens = TRUE)
  c(roughx = roughx, roughy=roughy)
  
}
