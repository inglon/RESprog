################################################################
################################################################
# Script: ProgR_6.1_Heterogeneity.R
# Author: Ingrid Lonnstedt
# Date:  30/08/2013
# R version: R 3.0.0
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
# Transform segment dataset
#
################################################################
################################################################
info = read.delim(paste(prefix.ann, '../../Infokey.txt', sep = ''), stringsAsFactors = F)
info$sample.name = format(info$sample.name, nsmall = 4, trim = T)
library(PSCBS)
library(Rsolnp)

###Arguments
alpha = .64
maxCN = 24

###Standard setup
A2 = 0:maxCN
A1 = matrix(A2, ncol = maxCN + 1, nrow = maxCN + 1)
A2 = t(A1)

E1 = alpha*A1 + (1-alpha)
E2 = alpha*A2 + (1-alpha)
EE = cbind(E1, E2)

###Sample
i = 15

####################
###Read PSCBS data
####################
d = file.path(getwd(), 'AROMA', 'reports', 'SNP_segmentation', 'CBSfits')
da = loadObject(paste(d, '/fit', i, '.RData', sep = ''))
mat = da$output
names(mat)[names(mat) == 'c1Mean'] = 'a1'
names(mat)[names(mat) == 'c2Mean'] = 'a2'
names(mat)[names(mat) == 'tcnStart'] = 'Start'
names(mat)[names(mat) == 'tcnEnd'] = 'End'
names(mat)[names(mat) == 'chromosome'] = 'Chromosome'
mat$length = mat$End - mat$Start
mat = subset(mat, !is.na(length) & Chromosome<23 & !is.na(a2) & !is.na(a1))

####################
###Or read HAPSEG output data
####################
source(file.path(getwd(), 'programs', 'make_HAPSEG_seg_obj.R'))
d = file.path(getwd(), 'ABSOLUTE', 'HAPSEG2')
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


###Plot TAPS for talk
x = mat$a1 + mat$a2
y = 2*(mat$a2/(x)-.5)
xlim = c(quantile(x, probs = c(0,.975), na.rm = T))
ylim = c(quantile(y, probs = c(0, 1), na.rm = T))
plot(x, y, pch = 16, col = col.transp('grey', .3), 
     xlim = xlim, ylim = ylim, xlab = 'Total CN intensity', 
     ylab = 'Allelic Imbalance 2*|BAF - 0.5|', main = 'TAPS')

xden = density(mat$a2+mat$a1,kernel="gaussian", bw = .02)
plot(xden$x, xden$y, type = 'n', xlim = xlim, 
     xlab = 'Total CN array intensity', yaxt = 'n')
polygon(xden$x, xden$y, col = col.transp(1, .5), border = NA)

####################
###Start values
####################
par(mfrow = c(1,2))
x = mat$a1 + mat$a2
y = 2*(mat$a2/(x)-.5)
xlim = c(quantile(x, probs = c(.025,.975), na.rm = T))
ylim = c(quantile(y, probs = c(0, 1), na.rm = T))
plot(x, y, pch = 16, col = col.transp('grey', .3), 
     xlim = xlim, ylim = ylim, xlab = 'Total CN intensity', 
     ylab = 'Allelic Imbalance 2*|BAF - 0.5|', main = 'TAPS')
title(main = info$sample.name[i], outer = T, line = -1)
#xlim = c(0, max(mat$a2, na.rm = T))
#ylim = xlim
xlim = c(0, quantile(mat$a2, probs = c(.9)))
ylim = c(0, quantile(mat$a1, probs = c(.975)))
plot(mat$a2, mat$a1, pch = 16, col = col.transp('grey', .3), 
     xlim = xlim, ylim = ylim, xlab = 'a2', ylab = 'a1',
     main = 'Minor and major CN intensity')
abline(0,1, col = 'red')

print('Mark two subsequent CN clusters along a "vertical" line:')
pos = locator(2)
dy = abs(diff(pos$y))
print('Mark two distant points in one "vertical" line:')
pos = locator(2)
l = lm(x~y, data = pos)$coef[2]
print('Mark two subsequent CN clusters along a horizontal line:')
pos = locator(2)
dx = abs(diff(pos$x))
print('Mark two distant points in one "horizontal" line:')
pos = locator(2)
k = lm(y~x, data = pos)$coef[2]

f11 = dy
f22 = dx

F = matrix(c(f11, l*f11, k*f22, f22), ncol = 2, nrow = 2)

### Starting value positions of lattice points
e1 = F[1,1] * E1 + F[1,2]*E2
e2 = F[2,1] * E1 + F[2,2]*E2
ee1 = c(e1)[!is.na(c(e1))]
ee2 = c(e2)[!is.na(c(e2))]

###Assign observations to closest lattice point
x = cbind(mat$a1, mat$a2)
alldist = apply(x, 1, dist, y = cbind(ee1, ee2))
closest = apply(alldist, 2, which.min)
#points(mat$a2, mat$a1, pch = 16, col = col.transp(closest, .3), cex = .7)
points(ee2, ee1)


####################
###Optimize rotation fit
####################

###Optimization starts here
w = 1#log(mat$W*1000)
rhofunc = 'leastsq'
rhoargs = list(k = 1, c = 1)
x = cbind(mat$a1, mat$a2)
pars = solnp(pars = c(F), fun = rho, x = x,
             w = w, EE = EE, rhofunc = rhofunc, 
             rhoargs = rhoargs)$pars
Fout = matrix(pars, ncol = 2, nrow = 2)

###Plot results in original scale
xlim = c(quantile(mat$a2, probs = c(.025,.975)))
ylim = c(quantile(mat$a1, probs = c(.025,.975)))

plot(mat$a2, mat$a1, pch = 16, col = col.transp('grey', .3), 
     xlim = xlim, ylim = ylim, xlab = 'a2', ylab = 'a1')
abline(0,1, col = 'red')
title(main = i)
ee1out = Fout[1,1] * EE[,1] + Fout[1,2]*EE[,2]
ee2out = Fout[2,1] * EE[,1] + Fout[2,2]*EE[,2]
alldist = apply(x, 1, dist, y = cbind(ee1out, ee2out))
closest = apply(alldist, 2, which.min)
#points(mat$a2, mat$a1, pch = 16, col = col.transp(closest, .3), cex = .7)
#points(ee2, ee1)
points(ee2out, ee1out, col = 'blue')


###########################
#Display results in theoretical scales
###########################


###Plot with colouring according to cluster ee
Fm1 = solve(Fout)
a1p = mat$a1*Fm1[1,1] + mat$a2 * Fm1[1,2]
a2p = mat$a1*Fm1[2,1] + mat$a2 * Fm1[2,2]
xlim = c(0, max(c(a2p, EE[,2])*1.2, na.rm = T))
ylim = c(0, max(a1p, na.rm = T)*1.2)
#xlim = c(quantile(a2p, probs = c(.01,.99)))
#ylim = c(quantile(a1p, probs = c(.01,.99)))
plot(a2p, a1p, pch = 16, col = col.transp('grey', .3), 
     xlim = xlim, ylim = ylim, xlab = 'a2', ylab = 'a1')
abline(0,1, col = 'red')
title(main = i)

alldist = apply(cbind(a1p, a2p), 1, dist, y = cbind(EE[,1], EE[,2]))
closest = apply(alldist, 2, which.min)
points(a2p, a1p, pch = 16, col = col.transp(closest, .3), cex = .7)
points(EE[,2], EE[,1])

###Plot in theoretical scales with weight colouring

#Recalculate optimized fit weights
ee1 = Fout[1,1] * EE[,1] + Fout[1,2]*EE[,2]
ee2 = Fout[2,1] * EE[,1] + Fout[2,2]*EE[,2]
alldist = apply(x, 1, dist, y = cbind(ee1, ee2))
closest = apply(alldist, 2, which.min)
e = cbind(ee1, ee2)[closest,]
r = dist(x, e)*w
di = dist(x, e)
eachweight = do.call(rhofunc, args = list(x = r, args = rhoargs))

#Transform
Fm1 = solve(Fout)
a1p = mat$a1*Fm1[1,1] + mat$a2 * Fm1[1,2]
a2p = mat$a1*Fm1[2,1] + mat$a2 * Fm1[2,2]

#Colour palette
#display.brewer.all()
pal = brewer.pal(9, 'YlOrRd')
#display.brewer.pal(8, 'Set1')
cols = pal[-c(1:2)]
cols = colorRampPalette(cols)(100)

#Plot frame
xlim = c(0, max(c(a2p, EE[,2]), na.rm = T))
ylim = c(0, max(a1p, na.rm = T)*1.2)
plot(a2p, a1p, pch = 16, col = col.transp('grey', .3), 
     xlim = xlim, ylim = ylim, xlab = 'Transformed a2', 
     ylab = 'Transformed a1', type = 'n')
title(main = paste('Sample', i, rhofunc, 'minfunction'))
mtext('Weight = log(constant*segment length)')

colby = eachweight
colbreaks = seq(min(colby), max(colby), length = length(cols))
delta = colbreaks[2] - colbreaks[1]
colind = (colby-colbreaks[1]) %/% delta
points(a2p, a1p, pch = 16, col = col.transp(cols[colind], .5), 
       cex = 1)

#Gridlines
EEs = unique(EE[,1])
for (j in 1:length(EEs)){
  lines(range(EEs), rep(EEs[j], 2), lty = 'dotted', col = 'grey20')
  lines(rep(EEs[j], 2), range(EEs), lty = 'dotted', col = 'grey20')
}
abline(0,1)

#Colour scale
tmp = c(xlim[2]/2, xlim[2])
tmp = seq(tmp[1], tmp[2], length=length(cols))
points(tmp, rep(0, length(tmp)), col = cols, pch = 16)
text(xlim[2],(ylim[1]*97/100 + ylim[2]*3/100),
     'Contribution to minfunction', adj = 1)

par(mfrow=c(2,2))
plot(w, eachweight, log = 'xy')
plot(r, eachweight, log = 'xy')
plot(di, eachweight, log = 'xy', xlab = 'Distance to lattice point')
plot(di, r, log = 'xy', xlab = 'Distance to lattice point')
title(main = 'weight = log(constant * segment length)', outer = T, line = -1)



################################################################
################################################################
#
# Display different transformation versions
#
################################################################
################################################################
info = read.delim(paste(prefix.ann, '../../Infokey.txt', sep = ''), stringsAsFactors = F)
info$sample.name = format(info$sample.name, nsmall = 4, trim = T)
library(PSCBS)
library(Rsolnp)

###Arguments
alpha = .85
maxCN = 24

###Standard setup
A2 = 0:maxCN
A1 = matrix(A2, ncol = maxCN + 1, nrow = maxCN + 1)
A2 = t(A1)
ix.na = A1>A2
A1[ix.na] = NA
A2[ix.na] = NA

E1 = alpha*A1 + (1-alpha)
E2 = alpha*A2 + (1-alpha)
EE = cbind(EE1 = c(E1)[!is.na(c(E1))], EE2 = c(E2)[!is.na(c(E2))])


###Sample
i = 21

####################
###Read HAPSEG output data
####################
source(file.path(getwd(), 'programs', 'make_HAPSEG_seg_obj.R'))
d = file.path(getwd(), 'ABSOLUTE', 'HAPSEG2')
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

###################
###Optimization
###################
F = transform.startval(a1 = mat$a1, a2 = mat$a2, alpha = alpha)
w = 1#log(5*mat$W + 1)
rhofunc = 'leastsq'
rhoargs = list(k = 1, c = 0.1)
x = cbind(mat$a1, mat$a2)
pars = solnp(pars = c(F), fun = rho, x = x,
             w = w, EE = EE, rhofunc = rhofunc, 
             rhoargs = rhoargs, LB = 0.5*c(F[1,1], -Inf, -Inf, F[2,2]),
             UB = 1.5*c(F[1,1],Inf, Inf, F[2,2]))$pars
Fout = matrix(pars, ncol = 2, nrow = 2)
plot.transformed(a1=mat$a1, a2=mat$a2, w=w, Fout=Fout, EE=EE, rhofunc=rhofunc, 
                            rhoargs=rhoargs, color.by = 'log')

mtext('Weight = 1', adj = 0, cex = .8)
mtext('Colour by contribution to minfunc', adj = 1, cex = .8)
mtext('Weight = segment length proportion', adj = 0, cex = .8)
mtext('Weight = log(constant*segment length)', adj = 0, cex = .8)

plot.transformed(a1=mat$a1, a2=mat$a2, w=w, Fout=Fout, EE=EE, rhofunc=rhofunc, 
                 rhoargs=rhoargs, color.by = 'nocol',
                 xlim = c(0,3), ylim = c(0, 3))

plot.transformed(mat$a1, mat$a2, w=w, Fout=Fout, EE=EE, rhofunc=rhofunc, 
                 rhoargs=rhoargs, color.by = mat$W)
mtext('Colour by segment length', adj = 1, cex = .8)

#Recalculate individual influence of optimized solution
ee1 = Fout[1,1] * EE[,1] + Fout[1,2]*EE[,2]
ee2 = Fout[2,1] * EE[,1] + Fout[2,2]*EE[,2]
alldist = apply(x, 1, dist, y = cbind(ee1, ee2))
closest = apply(alldist, 2, which.min)
e = cbind(ee1, ee2)[closest,]
r = dist(x, e)*w
di = dist(x, e)
eachweight = do.call(rhofunc, args = list(x = r, args = rhoargs))

#Plot weights
par(mfrow=c(2,2))
plot(w, eachweight, log = 'xy')
plot(r, eachweight, log = 'xy')
plot(di, eachweight, log = 'xy', xlab = 'Distance to lattice point')
plot(di, r, log = 'xy', xlab = 'Distance to lattice point')
title(main = 'weight = log(constant * segment length)', outer = T, line = -1)

library(numDeriv)
grad
################################################################
################################################################
#
# Transformation figures
#
################################################################
################################################################

  
  x = c(0, 1, 2, 3, 4, 5, 0, 1, 2, 3, 4, 5, 0, 1, 2, 3, 4, 5, 0, 1, 2, 3, 4, 5, 0, 1, 2, 3, 4, 5)
  y = c(0,0,0,0,0,0,1,1,1,1,1,1,2,2,2,2,2,2,3,3,3,3,3,3,4,4,4,4,4,4)
  a = 2.02 
  b = -.3
  c = -.7
  d = 0.5
  xp = a*x + b*y
  yp = c*x + d*y
  xlim = range(c(x, xp))
  ylim = range(c(y, yp))
  plot(x, y, col = 1, xlim = xlim, ylim = ylim)

  points(xp, yp, col = 4, pch = 16)
  
  
  cr = 0:10
  x = cr*(1+2)/(1+cr*2)
  plot(cr, x)
  
  cr = 0:10
  x = cr*(1+0)/(1+cr*0)
  plot(cr, x)
  abline(0,1)
  #

################################################################
################################################################
#
# Test
#
################################################################
################################################################


F = matrix(c(1, -.2,
             0.1, 1), byrow = T, ncol = 2)
x = rep(0:5, each = 6)
y = rep(0:5, 6)
plot(x, y, xlim = c(-2, 8), ylim = c(-2, 8), type = 'n')
for (j in length(x)) abline(v = x, col = 'grey', lty = 'dotted')
for (j in length(y)) abline(h = x, col = 'grey', lty = 'dotted')
abline(0,1, col = 'red')
ex = F[1,1]*x + F[1,2]*y
ey = F[2,1]*x + F[2,2]*y
points(ex, ey, col = 'blue')

#Do diagonal clusters violate the grid? 
#Does PSCBS look better than HAPSEG?

################################################################
################################################################
#
# Test
#
################################################################
################################################################
###Arguments
alpha = .39
maxCN = 24

###Standard setup
A2 = 0:maxCN
A1 = matrix(A2, ncol = maxCN + 1, nrow = maxCN + 1)
A2 = t(A1)

E1 = alpha*A1 + (1-alpha)
E2 = alpha*A2 + (1-alpha)

plot(E2+E1, 2*(E2/(E1+E2)-.5), xlim = c(0, 8), ylim = c(0, 1), xlab = 
       'Total copy number intensity',xaxt = 'n',
     ylab = 'Allelic Imbalance 2*|BAF-0.5|', main = 'TAPS')
plot(E2, E1, xlim = c(0, 3), ylim = c(0, 3), xlab = 'Major copy number intensity',
     ylab = 'Minor copy number intensity', 
     main = 'Minor and major CN intensity',
     yaxt = 'n', xaxt = 'n')


################################################################
################################################################
#
# Without alpha
#
################################################################
################################################################
info = read.delim(file.path(getwd(), 'AROMA','annotationData','Infokey.txt'), stringsAsFactors = F)
info$sample.name = format(info$sample.name, nsmall = 4, trim = T)
library(Rsolnp)


###Arguments
i = 21


####################
###Read HAPSEG output data
####################
source(file.path(getwd(), 'programs', 'make_HAPSEG_seg_obj.R'))
d = file.path(getwd(), 'ABSOLUTE', 'HAPSEG2')
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
###Read PSCBS data
####################
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


####################
###Start values
####################
maxlines = 30
maxCN = 24

par(mfrow = c(1,2))
x = mat$a1 + mat$a2
y = 2*(mat$a2/(x)-.5)
xlim = c(quantile(x, probs = c(.025,.975), na.rm = T))
ylim = c(quantile(y, probs = c(0, 1), na.rm = T))
xlim = c(0,2)
plot(x, y, pch = 16, col = col.transp('grey', .3), 
     xlim = xlim, ylim = ylim, xlab = 'Total CN intensity', 
     ylab = 'Allelic Imbalance 2*|BAF - 0.5|', main = 'TAPS')
title(main = info$sample.name[i], outer = T, line = -1)
#xlim = c(0, max(mat$a2, na.rm = T))
#ylim = xlim
xlim = c(0, quantile(mat$a2, probs = c(.9)))
ylim = c(0, quantile(mat$a1, probs = c(.975)))
plot(mat$a2, mat$a1, pch = 16, col = col.transp('grey', .3), 
     xlim = xlim, ylim = ylim, xlab = 'a2', ylab = 'a1',
     main = 'Minor and major CN intensity')
abline(0,1, col = 'red')

####################
###Optimize fit
####################
mystart = start.alphamax.f(mat)


###Optimization starts here
library(Rsolnp)
w = 1#log(mat$W*1000)
x = cbind(mat$a1, mat$a2)
rhofunc = 'tukey'
rhoargs = list(k = 1, c = 1)
pars = solnp(pars = c(mystart$alphamax, mystart$F), fun = rho.alphamax.f, x = x,
             w = w, maxCN = 24, rhofunc = rhofunc, 
             rhoargs = rhoargs, UB = c(1, mystart$UBf), 
             LB = c(0,mystart$LBf))$pars
alpha = pars[1]
F = matrix(pars[-1], ncol = 2, nrow = 2)


####################
### Look at fit
####################

xlim = NULL
ylim = NULL
maxCN = 24
color.by = 'eachweight'


  require('RColorBrewer')
  
  #Recalculate transformed lattice points from F and alpha
  Fm1 = solve(F)  
  A2 = 0:maxCN
  A1 = matrix(A2, ncol = maxCN + 1, nrow = maxCN + 1)
  A2 = t(A1)
  
  E1 = c(alpha*A1 + (1-alpha))
  E2 = c(alpha*A2 + (1-alpha))
  
  ### Starting value positions of lattice points in observed scale
  e1 = F[1,1] * E1 + F[1,2]*E2
  e2 = F[2,1] * E1 + F[2,2]*E2
  alldist = apply(x, 1, dist, y = cbind(e1, e2))
  closest = apply(alldist, 2, which.min)
  e = cbind(e1, e2)[closest,]
  r = dist(x, e)/w
  eachweight = do.call(rhofunc, args = list(x = r, args = rhoargs))
  
  #Transform
  a1p = x[,1]*Fm1[1,1] + x[,2] * Fm1[1,2]
  a2p = x[,1]*Fm1[2,1] + x[,2] * Fm1[2,2]
  
  #Colour palette
  #display.brewer.all()
  pal = brewer.pal(9, 'YlOrRd')
  #display.brewer.pal(8, 'Set1')
  cols = pal[-c(1:2)]
  cols = colorRampPalette(cols)(100)
  

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
colbreaks = seq(min(colby), max(colby), length = length(cols))
delta = colbreaks[2] - colbreaks[1]
colind = (colby-colbreaks[1]) %/% delta + 1
mycols = col.transp(cols[colind], .5)
if (color.by == 'nocol') mycols = col.transp('grey', .3)
points(a2p, a1p, pch = 16, col = mycols, cex = 1)

#Gridlines
for (j in 1:length(EEs)){
  lines(range(EEs), rep(EEs[j], 2), lty = 'dotted', col = 'grey20')
  lines(rep(EEs[j], 2), range(EEs), lty = 'dotted', col = 'grey20')
}
abline(0,1)

#Single scale density y
screen(2)
par(mar = c(0,0,0,0))
EEs = unique(E1)
yden = density(a1p,kernel="gaussian", bw = .03)
roughy = round(sum(diff(xden$y)**2)/diff(xden$x[1:2]), 3)
yden$y = log(yden$y+1)
plot(1, 1, ylim = ylim, xlim = range(yden$y), type = 'n', xaxt = 'n')
polygon(yden$y, yden$x, col = col.transp(1, .5), border = NA)

#Single scale density x
screen(3)
par(mar = c(0,0,0,0))
xden = density(a2p,kernel="gaussian", bw = .03)
roughx = round(sum(diff(xden$y)**2)/diff(xden$x[1:2]), 3)
xden$y = log(xden$y+1)
plot(1, 1, xlim = xlim, ylim = range(yden$y), type = 'n', yaxt = 'n')
polygon(xden$x, xden$y, col = col.transp(1, .5), border = NA)

#I am here trying to type the roughness in screen 4...
screen(4)
par(mar = c(0,0,0,0))
plot(1,1, type = 'n', xaxt = 'n', yaxt = 'n', xlim = c(0,1), ylim = c(0,1), bty = 'n')
text(0, 1, roughx, pos = 1, offset = 0)
text(0, sum(range(axTicks(2))*c(.9,.1)), roughy, srt = 90, pos = 3, offset = 0)


       xlab = 'Transformed a2', 
       ylab = 'Transformed a1', )

plot(1:100, rnorm(100), xaxs = "i", ylim = c(-3, 3), yaxt = "n", col = "red")

screen(3)

par(mar = c(0, 0, 0, 0))
plot(1:100, rnorm(100), xaxs = "i", ylim = c(-3, 3), yaxt = "n", col = "red")
title(main = paste('Sample', i, rhofunc, 'minfunction'), outer = T, line = -1)
close.screen(all.screens = TRUE)

  #Single scale densities
  
  #diff(xden$y)/diff(xden$x[1:2])
  
  
  
  
  #Colour scale
  if (color.by != 'nocol'){
    tmp = c(xlim[2]/2, xlim[2])
    tmp = seq(tmp[1], tmp[2], length=length(cols))
    points(tmp, rep(-.05, length(tmp)), col = cols, pch = 16)
  }
  



################################################################
################################################################
#
# Looking at datasets
#
################################################################
################################################################
info = read.delim(file.path(getwd(), 'AROMA','annotationData','Infokey.txt'), stringsAsFactors = F)
info$sample.name = format(info$sample.name, nsmall = 4, trim = T)
library(Rsolnp)


par(ask = T, mfrow = c(1,2))
for ( i in 2:56){

####################
###Read HAPSEG output data
####################
source(file.path(getwd(), 'programs', 'make_HAPSEG_seg_obj.R'))
d = file.path(getwd(), 'ABSOLUTE', 'HAPSEG2')
d = file.path(getwd(), 'ABSOLUTE', 'WithNormals', 'HAPSEGbirdseed')
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
     ylab = 'Allelic Imbalance 2*|BAF - 0.5|', main = 'TAPS')
title(main = info$sample.name[i], outer = T, line = -1)
#xlim = c(0, max(mat$a2, na.rm = T))
#ylim = xlim
xlim = c(0, quantile(mat$a2, probs = c(.9)))
ylim = c(0, quantile(mat$a1, probs = c(.975)))
plot(mat$a2, mat$a1, pch = 16, col = col.transp('grey', .3), 
     xlim = xlim, ylim = ylim, xlab = 'a2', ylab = 'a1',
     main = 'Minor and major CN intensity')
abline(0,1, col = 'red')
}





### Look at HAPSEGbirdseed
#
i = 15
source(file.path(getwd(), 'programs', 'make_HAPSEG_seg_obj.R'))
d = file.path(getwd(), 'ABSOLUTE', 'WithNormals', 'HAPSEGbirdseed')
arrname = info$arrname[i]
d = paste(d,'/', arrname, sep = '')
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
     ylab = 'Allelic Imbalance 2*|BAF - 0.5|', main = 'TAPS')
title(main = info$sample.name[i], outer = T, line = -1)
#xlim = c(0, max(mat$a2, na.rm = T))
#ylim = xlim
xlim = c(0, quantile(mat$a2, probs = c(.9)))
ylim = c(0, quantile(mat$a1, probs = c(.975)))
plot(mat$a2, mat$a1, pch = 16, col = col.transp('grey', .3), 
     xlim = xlim, ylim = ylim, xlab = 'a2', ylab = 'a1',
     main = 'Minor and major CN intensity')
abline(0,1, col = 'red')




