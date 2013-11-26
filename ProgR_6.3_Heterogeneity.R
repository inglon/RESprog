################################################################
################################################################
# Script: ProgR_6.3_Heterogeneity.R
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
info = read.delim(file.path(getwd(), 'AROMA','annotationData','Infokey.txt'), stringsAsFactors = F)
info$sample.name = format(info$sample.name, nsmall = 4, trim = T)
library(PSCBS)
source(file.path(getwd(), 'programs', 'make_HAPSEG_seg_obj.R'))


d = file.path(getwd(), 'ABSOLUTE', 'WithNormals', 'SampleData')
pdf(file.path(d, 'TAPS_after_HAPSEGS.pdf'), 
    paper = 'a4r', height = 6, width = 12)  
for (i in 1:56){
  d = file.path(getwd(), 'ABSOLUTE', 'WithNormals', 'SampleData')
  d = paste(d,'/', info$sample.name[i], sep = '')
  dir.create(d, recursive = T, showWarnings = F)
  #  pdf(file.path(d, paste(info$sample.name[i], '_TAPS_after_HAPSEGS.pdf', sep = '')), 
  #      paper = 'a4r', height = 6, width = 12)  
  par(mfrow = c(1,2), mar = c(4, 4, 4, .1), oma = c(0, 0, 1, 0))
  for (h in 1:2){
    dh = paste('HAPSEG', h, sep = '')
    d = file.path(getwd(), 'ABSOLUTE', 'WithNormals', dh)
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
  #  dev.off()
}
dev.off()


################################################################
################################################################
#
# Copy number estimation
#
################################################################
################################################################
mysamples = c(24, 30, 56, 36, 3, 21)
#Load a sample
info = read.delim(file.path(getwd(), 'AROMA','annotationData','Infokey.txt'), stringsAsFactors = F)
info$sample.name = format(info$sample.name, nsmall = 4, trim = T)
i = 24
d = file.path(getwd(), 'HET', info$sample.name[i])
load(file = file.path(d, paste(info$sample.name[i], '.RData', sep = '')))


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
dev.print(png, file=file.path(d, paste(info$sample.name[i],'_1_Standardized_TAPS.png', 
                                       sep = '')),  width=640, height=640)

####################
###Read HAPSEG output data
####################
h = 2
if (i %in% c(1:3, 6, 10, 22, 28, 35, 42, 46, 48)) h = 1


source(file.path(getwd(), 'programs', 'make_HAPSEG_seg_obj.R'))
d = file.path(getwd(), 'ABSOLUTE', 'WithNormals','HAPSEG')
d = paste(d, h, sep = '')
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
d = file.path(getwd(), 'HET', info$sample.name[i], 
              paste(info$sample.name[i], '_2_Start_values', sep = ''))
mystart = start.alphamax.f(mat, printfile = d, colbychrom = T)
alpha = mystart$alphamax
F = mystart$F
plot.transformed(mat, alpha, 
                 F, xlim = c(-5,20), ylim = c(-5,8))
d = file.path(getwd(), 'HET', info$sample.name[i], 
              paste(info$sample.name[i], '_2_Start_values_4.png', sep = ''))
dev.print(png, file = d, width = 640, height = 500)

plot.transformed(mat, alpha, 
                 F, xlim = c(0,7), ylim = c(0,4))
d = file.path(getwd(), 'HET', info$sample.name[i], 
              paste(info$sample.name[i], '_2_Start_values_5.png', sep = ''))
dev.print(png, file = d, width = 640, height = 500)

####################
###Evaluate rho function and weights
####################

#Input
w = 1-exp(-mat$length/500000) #Sample 24, 30, 36, 21
w = 1-exp(-mat$length/5000000) #Sample 56, 3
rhofunc = 'tukey'
nmad = 3 #Sample 24, 30, 56, 36


#Change weight constant until 1'000'000 length segs have almost top weight
#Change nmod until Contrib to rho  gets constant at good Distance to lattice point
evaluate.weights(mat, mystart$alphamax, mystart$F, nmad = nmad, w=w)
d = file.path(getwd(), 'HET', info$sample.name[i])
dev.print(png, file=file.path(d, paste(info$sample.name[i],'_3.1_evaluate_weights.png', 
                                       sep = '')),  width=900, height=900)


####################
###Optimization
####################

#w, rhofunc and nmad from above!
di = getdist(mat, mystart$alphamax, mystart$F) 
rhoargs = list(k = mad(di), c = mad(di)*nmad)
x = cbind(mat$a1, mat$a2)

pars = solnp(pars = c(mystart$alphamax, mystart$F), fun = rho.alphamax.f, x = x,
             w = w, maxCN = 24, rhofunc = rhofunc, 
             rhoargs = rhoargs)$pars

alpha = pars[1]
F = matrix(pars[-1], ncol = 2, nrow = 2)
plot.transformed(mat, alpha, F, w=w, rhoargs=rhoargs, rhofunc = rhofunc, 
                 xlim = c(0,7), ylim = c(0,4))
d = file.path(getwd(), 'HET', info$sample.name[i])
dev.print(png, file=file.path(d, paste(info$sample.name[i],'_3.2_rotated_grid.png', 
                                       sep = '')),  width=640, height=500)

#All possible alphas and Fs
maxCN = 24
alpha = pars[1]
F = matrix(pars[-1], ncol = 2, nrow = 2)

#All possible alphas
EEs = alpha*c(0:maxCN) + (1-alpha) 
alphas = NULL
fs = NULL
for (j in 1:length(EEs)){
  alphas = c(alphas, alpha/(EEs[j]+alpha))
  if (j == 1) fs = 1 else {
    fs = c(fs, fs[j-1]*alphas[j-1]/alphas[j])
  }
}

####################
### Choose alpha with help from mutations
####################


###Read mutations
maf.fn = paste(getwd(), 
               '/AROMA/rawData/responsify/Exome/MAF/AllVarintsMutectOncotatorPlusMutsigCategNo5214_mindepth10_withvaf', '_', 
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
nrow(maf)
muts = unique(subset(maf, select = c('Chromosome','Start_position','End_position',
                                     't_alt_count','t_ref_count','segment', 'depth',
                                     'effect', 'Hugo_Symbol')))
nrow(muts)
#rm(maf)
muts$af = muts$t_alt_count/(muts$t_alt_count+muts$t_ref_count)
muts = subset(muts, !is.na(segment) )
muts = subset(muts, muts$depth>=20)
nrow(muts)

#Plot observed variant frequencies 
#NB: investigate n = 1:4 if surely gridrow 4 must correspond to at least 0 copies
par(mfrow = c(2,2))
for (n in 3:6){
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
    
  ord = order(both, major)
  plot(both[ord], ylim = c(-0.2,1), ylab = 'Allele fraction', xlab = 'CN segment', 
       main = paste('CN = 0 at gridline',n, 'from bottom/left'))
  points(major[ord], col = 'blue')
  points(minor[ord], col = 'red')
  muts$index = NA
  for (j in 1:nrow(muts)){
    if (!is.na(muts$segment[j])) muts$index[j] = 
      which(ord == muts$segment[j])
  }
  points(muts$index, muts$af, col = 'green', pch = 16)
  abline(h = 0)
}
d = file.path(getwd(), 'HET', info$sample.name[i])
dev.print(png, file=file.path(d, paste(info$sample.name[i],'_4.1_minor_major_af.png', 
                                       sep = '')),  width=900, height=900)


###Plot by genome position
par(ask = FALSE, mfrow = c(2,1))
for (n in 5:6){
  F = matrix(pars[-1], ncol = 2, nrow = 2)
  alpha = alphas[n]
  F = F*fs[n]
  Fm1 = solve(F)
  c1 = mat$a1*Fm1[1,1] + mat$a2 * Fm1[1,2]
  c2 = mat$a1*Fm1[2,1] + mat$a2 * Fm1[2,2]
  c1 = (c1 - (1 - alpha))/alpha
  c2 = (c2 - (1 - alpha))/alpha 
  mult = muts$af*(alpha*((c1 + c2)[muts$segment])+(1-alpha)*2)/alpha
   plot.CNmut(plotmat = mat, plotmuts = muts, c1 = c1, c2 = c2, mult = mult,
             main = paste('n =', n, ', purity', round(alpha*100), '%'),
              restrict = TRUE)
}
d = file.path(getwd(), 'HET', info$sample.name[i])
dev.print(png, file=file.path(d, paste(info$sample.name[i],'_4.2_CNmuts_preliminary.png', 
                                       sep = '')),  width=900, height=900)


#Choose n
n = 3 #Sample 24
n = 4 #Sample 21
n = 5 #Sample 30
n = 9 #Sample 3
n = 7 #Sample 56
alpha = alphas[n]
F = matrix(pars[-1], ncol = 2, nrow = 2)
plot.transformed(mat, alpha, F = F*fs[n], w=w, rhoargs=rhoargs, rhofunc = rhofunc, 
                 xlim = c(0,3), ylim = c(0,1.5))
mtext(paste('Purity', round(alpha*100), '%'),adj = 0.1, outer = T, line = -2)
d = file.path(getwd(), 'HET', info$sample.name[i])
dev.print(png, file=file.path(d, paste(info$sample.name[i],'_4.3_rotated_grid_chosen.png', 
                                       sep = '')),  width=640, height=500)


#Set transformed CNs
alpha = alphas[n]
F = matrix(pars[-1], ncol = 2, nrow = 2)
Fm1 = solve(F*fs[n])
mat$a1t = mat$a1*Fm1[1,1] + mat$a2 * Fm1[1,2]
mat$a2t = mat$a1*Fm1[2,1] + mat$a2 * Fm1[2,2]
mat$purity = alpha

#Save dataset details
details = list(sample = i, n=n, alpha=alpha, F=F, Fm1=Fm1)

####################
### Nonlinearity fix
####################

alpha = alphas[n]
dummy = nonlinAdj.asympOrig(mat$a1t, alpha = alpha, maxCN = 6)
dummy = nonlinAdj.asympOrig(mat$a2t, alpha = alpha, maxCN = 6)
dummy = nonlinAdj.logis(mat$a1t, alpha = alpha, maxCN = 6)#24, 30, 56
dummy = nonlinAdj.logis(mat$a2t, alpha = alpha, maxCN = 6)

mat$a1ta = nonlinAdj.logis(mat$a1t, alpha = alpha, maxCN = 6)#24, 30, 56
title(main = 'Minor CN logistic regression')
mat$a1ta = mat$a1t
mat$a2ta = nonlinAdj.asympOrig(mat$a2t, alpha = alpha, maxCN = 6) 
mat$a2ta = nonlinAdj.logis(mat$a2t, alpha = alpha, maxCN = 6)#24, 56
title(main = 'Major CN logistic regression')
mat$a2ta = mat$a2t 

d = file.path(getwd(), 'HET', info$sample.name[i])
dev.print(png, file=file.path(d, paste(info$sample.name[i],'_5.1_nonlinearlity_fix.png', 
                                       sep = '')),  width=900, height=900)


####Plots of results
####

#Calculate lattice points from alpha
alpha = alphas[n]
maxCN = 24
A2 = 0:maxCN
A1 = matrix(A2, ncol = maxCN + 1, nrow = maxCN + 1)
A2 = t(A1)
E1 = c(alpha*A1 + (1-alpha))
E2 = c(alpha*A2 + (1-alpha))

#Plot CN lattice
plot.cn(x = mat$a2ta, y = mat$a1ta, 
             xlim = c(0,4), ylim = c(0,2), grid = FALSE,
        xlab = 'Transformed a2 after non-linearity fix',
        ylab = 'Transformed a1 after non-linearity fix')
EEs = unique(E1)
for (j in 1:length(EEs)){
  lines(range(EEs), rep(EEs[j], 2), lty = 'dotted', col = 'grey20')
  lines(rep(EEs[j], 2), range(EEs), lty = 'dotted', col = 'grey20')
}

d = file.path(getwd(), 'HET', info$sample.name[i])
dev.print(png, file=file.path(d, paste(info$sample.name[i],'_5.2_rotated_grid_linearized.png', 
                                       sep = '')),  width=640, height=500)




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


smoothScatter(xcoord, ycoord, xlim = c(-alpha*3/4, alpha*3/4),
              ylim = c(-alpha*3/4, alpha*3/4), main = 'Not linearized')
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

smoothScatter(xcoord, ycoord, xlim = c(-alpha*3/4, alpha*3/4),
              ylim = c(-alpha*3/4, alpha*3/4), main = 'Linearized')
abline(v=0, h=0)
rect(xleft = -alpha/2, ybottom = -alpha/2, xright = alpha/2, 
     ytop = alpha/2, lty = 'dotted')

d = file.path(getwd(), 'HET', info$sample.name[i])
dev.print(png, file=file.path(d, 
                              paste(info$sample.name[i],'_5.3_overlay_lattice_points.png', 
                                    sep = '')),  width=640, height=400)

####################
### Assign decimal copy numbers, ploidy
####################

#Used a1t for sample 24, 30, 
#Used a1ta for sample 56, 21

#Calculate lattice points from alpha
alpha = alphas[n]
mat$c1 = (mat$a1t - (1 - alpha))/alpha
mat$c2 = (mat$a2t - (1 - alpha))/alpha
mat$alpha = alphas[n]

par(mfrow = c(1,1), mar = c(3, 3, 2, .1))
plot.cn(x = mat$c2, y = mat$c1, mycols = 'blue',
        xlim = c(0,6), ylim = c(-.5,4), xlab = 'Major CN', ylab = 'Minor CN')

#Calculate ploidy
mat$ploi = sum((mat$c1 + mat$c2)*mat$W)
mtext(paste('Ploidy =', format(round(mat$ploi[1],1),nsmall = 1)), adj = 0)

d = file.path(getwd(), 'HET', info$sample.name[i])
dev.print(png, file=file.path(d, 
                              paste(info$sample.name[i],'_5.4_CNs.png', 
                                    sep = '')),  width=400, height=400)

#Save details
details$ploidy = mat$ploi[1]



################################################################
################################################################
#
# Subclonal Mutations
#
################################################################
################################################################

#Plot variant fractions
par(mfrow = c(1,1))
alpha = alphas[n]
c1 = pmax(0, mat$c1)
c2 = pmax(0, mat$c2)

minor = alpha*c1/(alpha*(c1+c2) + (1-alpha)*2)
major = alpha*c2/(alpha*(c1+c2) + (1-alpha)*2)
both = alpha*(c1+c2)/(alpha*(c1+c2) + (1-alpha)*2)

ord = order(both, major)
plot(both[ord], ylim = c(0,1), ylab = 'Variant fraction', xlab = 'CN segment')
points(major[ord], col = 'blue')
points(minor[ord], col = 'red')

#Add mutations
muts$index = NA
for (j in 1:nrow(muts)){
  if (!is.na(muts$segment[j])) muts$index[j] = 
    which(ord == muts$segment[j])
}
points(muts$index, muts$af, col = 'green', pch = 16)

#Highlight subclonal mutations
muts$subclonal = subclonal.muts(muts, mat)
points(muts$index[muts$subclonal], muts$af[muts$subclonal], 
       col = 'darkgreen', pch = 16)
d = file.path(getwd(), 'HET', info$sample.name[i])
dev.print(png, file=file.path(d, 
                              paste(info$sample.name[i],'_7.1_minor_major_af.png', 
                                    sep = '')),  width=640, height=400)




################################################################
################################################################
#
# Subclonal CNs
#
################################################################
################################################################

#Overlay lattice points into X
maxCN = 24
E2 = 0:maxCN
E1 = matrix(E2, ncol = maxCN + 1, nrow = maxCN + 1)
E2 = t(E1)
E1 = c(E1)
E2 = c(E2)

x = cbind(mat$c1, mat$c2)
alldist = apply(x, 1, dist, y = cbind(E1, E2))
closest = apply(alldist, 2, which.min)
e = cbind(E1, E2)[closest,]
di = dist(x, e)
xcoord = x[,2]-e[,2]
ycoord = x[,1]-e[,1]
X = cbind(ycoord, xcoord)

#Colour scale for figures
pal = brewer.pal(11, 'Spectral')
pal = rev(pal[c(1:4, 9:11)])
#display.brewer.pal(9, 'RdBu')
#display.brewer.all()

par(mfrow = c(2,1))
#Multivariate t robust covariance qq=plots
nu = 2
D = mahalanobis.w(X, nu, w)  
cols = (colors.by(D, pal, span = c(0,10))$mycols)[order(D)]
qqplot(qexp(ppoints(length(D)), rate = .5), D, col = cols,
       main = 'Exponential QQ-plot',ylim = c(0,20),
       ylab = 'Mahalanobis squared distance',
       xlab = 'Exponential quantiles')
abline(0, 1, col = 'gray')
#Overlay lattice points
col1 = colors.by(D, pal, span = c(0,10))
plot(X[,2], X[,1], xlim = c(-.75, .75), col = col1$mycols, pch = 16,
     ylim = c(-.75, .75), main = '', cex.main = 1,
     xlab = 'Distance to closest major integer CN',
     ylab = 'Distance to closest minor integer CN')
mtext('Colouring by squared Mahalanobis distance to lattice point', cex = .8)
abline(v=0, h=0)
rect(xleft = -.5, ybottom = -.5, xright = .5, 
     ytop = .5, lty = 'dotted')
#Colour scale
tmp = c(0.1, .75)
tmp = seq(tmp[1], tmp[2], length=length(col1$cols))
points(tmp, rep(-.75, length(tmp)), col = col1$cols, pch = 16)
abline(v=0, h=0)
rect(xleft = -.5, ybottom = -.5, xright = .5, 
     ytop = .5, lty = 'dotted')
d = file.path(getwd(), 'HET', info$sample.name[i])
dev.print(png, file=file.path(d, 
                              paste(info$sample.name[i],'_6.1_CN_distribution.png', 
                                    sep = '')),  width=400, height=700)

mat$D = D


#Subclonal CN segments given cutoff with 1-pexp(cut, .5) sign.level
cut = 3 #Sample 24 and others
cols = rep('blue', nrow(mat))
cols[mat$D > cut] = 'red'
par(mfrow = c(1,1))
plot(X[,2], X[,1], xlim = c(-.75, .75), col = cols, pch = 16,
     ylim = c(-.75, .75), main = paste('Subclonal CNs with cutoff =', 
                                       cut), cex.main = 1,
     xlab = 'Distance to closest major integer CN',
     ylab = 'Distance to closest minor integer CN')
abline(v=0, h=0)
rect(xleft = -.5, ybottom = -.5, xright = .5, 
     ytop = .5, lty = 'dotted')
legend('topleft', bty = 'n', pch = 16, col = c('blue','red'), 
       legend = c('Segment with integer CN','Segment with subclonal CN'))
d = file.path(getwd(), 'HET', info$sample.name[i])
dev.print(png, file=file.path(d, 
                              paste(info$sample.name[i],'_6.2_subclonal_CNs_overlaid.png', 
                                    sep = '')),  width=400, height=400)

#Choose how extreme CNs can be judged to be subclonal: CN subclonality defined here
plot.cn(x = mat$c2, y = mat$c1, mycols = cols, xlim = c(0,10), ylim = c(-2,4))
cn.cut = 5.5 #21 Choose this!
cn.cut = 4.5 #24, 56 Choose this!
cn.cut = 3.5 #30 Choose this!
cn.cut0 = -1
mat$subclonal = FALSE
mat$subclonal[mat$D>cut] = TRUE
mat$subclonal[mat$c2>cn.cut | mat$c1>cn.cut] = FALSE
mat$subclonal[mat$c2<cn.cut0 | mat$c1<cn.cut0] = FALSE

#Subclonal segments wrt mutations
mat$subclonal.muts = subclonal.segmuts(muts = muts,
                                      mat = mat, cn.cut = cn.cut, cn.cut0 = cn.cut0)


#% ITH wrt CNs
cols = rep('blue', nrow(mat))
cols[mat$subclonal.muts] = 'orange'
cols[mat$subclonal] = 'red'
par(mfrow = c(1,1), mar = c(3, 3, 2, .1))
plot.cn(x = mat$c2, y = mat$c1, mycols = cols, xlim = c(0,8), ylim = c(-.5,4),
        xlab = 'Major CN', ylab = 'Minor CN')
mat$het = sum(mat$length[mat$subclonal | mat$subclonal.muts])/sum(mat$length)
mtext(paste(round(mat$het[1]*100), '% heterogeneity'), adj = 1)
mtext(paste('Ploidy =', format(round(mat$ploi[1],1), nsmall = 1)), adj = 0)
legend('topleft', bty = 'n', pch = 16, col = c('blue','orange', 'red'), 
       legend = c('Segment with integer CN',
                  'Segment with subclonal mutations',
                  'Segment with subclonal CN'))
d = file.path(getwd(), 'HET', info$sample.name[i])
dev.print(png, file=file.path(d, 
                              paste(info$sample.name[i],'_6.3_subclonal_CNs.png', 
                                    sep = '')),  width=400, height=400)


details$het = mat$het[1]
details$cn.cut = cn.cut
details$cn.cut0 = cn.cut0
details$cut = cut
################################################################
################################################################
#
# Plot CNs and mutations
#
################################################################
################################################################

#Plot CN grid with mutations
cols = rep('blue', nrow(mat))
cols[mat$subclonal.muts] = 'orange'
cols[mat$subclonal] = 'red'
cols = col.transp(cols, .5)
plot.cn(x = mat$c2, y = mat$c1, mycols = cols, xlim = c(0,8), ylim = c(-.5,4),
        xlab = 'Major CN', ylab = 'Minor CN')
mtext(paste(round(mat$het[1]*100), '% heterogeneity'), adj = 1)
mtext(paste('Ploidy =', format(round(mat$ploi[1],1), nsmall = 1)), adj = 0)
points(mat$c2[maf$segment], mat$c1[maf$segment], pch = '.', cex = 4)
points(mat$c2[(muts$segment)[muts$subclonal]], mat$c1[(muts$segment)[muts$subclonal]], 
       pch = 4, cex = 1.5)
legend('topleft', bty = 'n', pch = c(16, 16, 16, 4, 20), col = c('blue','orange', 'red', 'black', 'black'), 
       legend = c('Segment with integer CN', 'Segment with subclonal mutations', 'Segment with subclonal CN',
                  'Clonal mutation', 'Subclonal mutation'))
d = file.path(getwd(), 'HET', info$sample.name[i])
dev.print(png, file=file.path(d, 
                              paste(info$sample.name[i],'_7.2_subclonal_CNs_withMuts.png', 
                                    sep = '')),  width=400, height=400)





#Express vaf as number of copies per cell (cellular multiplicity)
alpha = alphas[n]
muts$mult = muts$af/alpha*
  (alpha*((mat$c1 + mat$c2)[muts$segment])+(1-alpha)*2)


#Colour settings
mutcols = c('green3','orange')
mycols = rep('blue', nrow(mat))
mycols[mat$subclonal.muts] = 'orange'
mycols[mat$subclonal] = 'red'


#Plot Copy numbers and multiplicities whole genome
plot.CNmut(plotmat = mat, plotmuts = muts[muts$effect == 'nonsilent',], mycols = mycols,
           mutcols = mutcols)
mtext(paste(nrow(muts[muts$effect == 'nonsilent',]), 'nonsilent mutations'))
mtext(paste(round(mat$het[1]*100), '% heterogeneity'), adj = 1)
mtext(paste('Ploidy =', format(round(mat$ploi[1],1), nsmall = 1)), adj = 0)

d = file.path(getwd(), 'HET', info$sample.name[i])
dev.print(png, file=file.path(d, 
              paste(info$sample.name[i],'_7.3_Genome_CNs_nonsilentMutations.png', 
                  sep = '')),  width=940, height=600)

plot.CNmut(plotmat = mat, plotmuts = muts, mycols = mycols, mutcols = mutcols)
mtext(paste(nrow(muts), 'mutations'))
mtext(paste(round(mat$het[1]*100), '% heterogeneity'), adj = 1)
mtext(paste('Ploidy =', format(round(mat$ploi[1],1), nsmall = 1)), adj = 0)
dev.print(png, file=file.path(d, 
                              paste(info$sample.name[i],'_7.4_Genome_CNs_allMutations.png', 
                                    sep = '')),  width=940, height=600)



#Plot Copy numbers and multiplicities chromosome by chromosome
plotmat = mat

d = file.path(getwd(), 'HET', info$sample.name[i])
pdf(paste(d, '/', info$sample.name[i], '_7.5_CNs_and_Mutations_byChromosome.pdf', 
          sep = ''), 
    paper = 'a4r', height = 6, width = 12)  
ylim = c(-1, 8)
plotmat$c1 = pmax(plotmat$c1, 0)
plotmat$c2 = pmax(plotmat$c2, 0)
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
    lines(c(plotmat$Start.bp[k], plotmat$End.bp[k]), rep(plotmat$tot[k], 2), 
          lty  = 'solid', lwd = 7, col = 'grey', lend = 'butt')
    lines(c(plotmat$Start.bp[k], plotmat$End.bp[k]), rep(plotmat$c1[k], 2), 
          lty  = 'solid', lwd = 7, col = mycols[k], lend = 'butt')
    lines(c(plotmat$Start.bp[k], plotmat$End.bp[k]), rep(plotmat$c2[k], 2), 
          lty  = 'solid', lwd = 7, col = mycols[k], lend = 'butt')
  }
  ix = which(as.character(muts$Chromosome) == chrom &
               muts$effect != 'nonsilent')
  if (length(ix)>0){
    cols = rep(mutcols[1], length(ix))
    cols[muts$subclonal[ix]] = mutcols[2]
    points(muts$Start_position[ix], pmin(muts$mult[ix], ylim[2]), pch = 2,
           col = cols)
  }
  ix = which(as.character(muts$Chromosome) == chrom &
               muts$effect == 'nonsilent')
  if (length(ix)>0){
    cols = rep(mutcols[1], length(ix))
    cols[muts$subclonal[ix]] = mutcols[2]
    text(muts$Start_position[ix], pmin(muts$mult[ix], ylim[2]),
         col = cols, labels = muts$Hugo_Symbol[ix], cex = .7)
  }
  
}
dev.off()
par(ask = FALSE)


################################################################
################################################################
#
# Plot CNs and mutations
#
################################################################
################################################################


####################
### Output dataset details
####################

out = data.frame(ind = details$sample, sample = info$sample[i],
                 purity = details$alpha)

####################
### Output workspace
####################

#Save workspace
d = file.path(getwd(), 'HET', info$sample.name[i])
save.image(file = file.path(d, paste(info$sample.name[i], '.RData', sep = '')))

#Load
#d = file.path(getwd(), 'HET', info$sample.name[i])
#load(file = file.path(d, paste(info$sample.name[i], '.RData', sep = '')))

####################
### Output HAPSEG transformed version dataset
####################

seg.dat.new = hapseg.format(seg.dat, mat, xname = 'c2', 
                            yname = 'c1')

#Print the new object
d = file.path(getwd(), 'ABSOLUTE', 'WithNormals', 'HAPSEG2')
d = paste(d,'/', info$sample.name[i], sep = '')
seg.dat = seg.dat.new
save(seg.dat, file = file.path(d, paste(info$sample.name[i], '_', info$arrname[i], 
                                        '_scaled_segdat.RData', sep = '')))



################################################################
################################################################
#End of file
