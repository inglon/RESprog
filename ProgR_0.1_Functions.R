library(IRanges)
library(Rsolnp)
library(MASS)
require('RColorBrewer')


##################################################################
##################################################################
##
## Load or read sample data
##
##################################################################
##################################################################
loadSample = function(i = NULL, info = NULL){
  d = file.path(getwd(), 'HET', info$sample.name[i])
  load(file = file.path(d, paste(info$sample.name[i], '_obj.RData', sep = '')))
  obj
}

readCBS = function(obj, printfolder = NULL){
  library(PSCBS)
  d = file.path(getwd(), 'AROMA', 'reports', 'SNP_segmentation', 'CBSfitsWithNormals')
  da = loadObject(paste(d, '/fit', obj$i, '.RData', sep = ''))
  mat = da$output
  names(mat)[names(mat) == 'c1Mean'] = 'a1'
  names(mat)[names(mat) == 'c2Mean'] = 'a2'
  names(mat)[names(mat) == 'tcnStart'] = 'Start'
  names(mat)[names(mat) == 'tcnEnd'] = 'End'
  names(mat)[names(mat) == 'chromosome'] = 'Chromosome'
  mat$length = mat$End - mat$Start
  mat = subset(mat, !is.na(length) & Chromosome<23 & !is.na(a2) & !is.na(a1))
  mat$W = mat$length/sum(mat$length)
  obj$matCBS = mat
  
  
  #Figure
  x = pmax(mat$a1,0) + mat$a2
  y = 2*(mat$a2/(x)-.5)
  xlim = c(quantile(x, probs = c(0,.975), na.rm = T)*c(1, 1.2))
  ylim = c(quantile(y, probs = c(0, 1), na.rm = T))
  plot(x, y, pch = 16, col = col.transp('blue', .3), 
       xlim = xlim, ylim = ylim, xlab = 'Total CN intensity', 
       ylab = 'Allelic Imbalance 2*|BAF - 0.5|', main = 'TAPS')
  abline(v = 1, lty = 'dotted')
  points(x[mat$a1<0], y[mat$a1<0], col = col.transp('red', .5))
  tx = ifelse((obj$info$ind[i] + 100) %in% obj$info$ind & obj$i !=54, 'Sample has germline SNP array',
              'Sample does not have germline SNP array')
  mtext(tx)
  if (!is.null(printfolder)) {
    dir.create(printfolder, recursive = T, showWarnings = FALSE)
    dev.print(png, file=file.path(printfolder, paste(obj$info$sample.name[i],'_1.2_Standardized_TAPS.png', 
                            sep = '')),  width=640, height=640)
  }
  obj
}

readHAPSEG = function(obj, printfolder = NULL){
  info = obj$info
  i = obj$i
  source(file.path(getwd(), 'programs', 'make_HAPSEG_seg_obj.R'))
  d = file.path(getwd(), 'ABSOLUTE', 'WithNormals','HAPSEG')
  d = paste(d, obj$h, sep = '')
  d = paste(d,'/', info$sample.name[i], sep = '')
  seg.dat.fn = file.path(d, paste(info$sample.name[i], '_', info$arrname[i], '_segdat.RData', sep = ''))
  seg.dat = HAPSEGMakeSegObj(seg.dat.fn, gender = 'Female', filter_segs=FALSE, min_probes=10, max_sd=100, 
                             verbose=TRUE)  
  as.seg.dat = seg.dat$as.seg.dat
  nr = as.integer(nrow(as.seg.dat))
  mat = as.seg.dat[1:(nr/2),]
  names(mat)[6] = 'a1'
  mat$a2 = as.seg.dat$copy_num[(nr/2+1):(nr)]
  mat$W = mat$W*2
  obj$mat = mat
  obj$seg.dat = seg.dat
  
  x = mat$a1 + mat$a2
  y = 2*(mat$a2/(x)-.5)
  xlim = c(quantile(x, probs = c(.025,.975), na.rm = T)*c(.8, 1.2))
  ylim = c(quantile(y, probs = c(0, 1), na.rm = T))
  plot(x, y, pch = 16, col = col.transp('blue', .3), 
       xlim = xlim, ylim = ylim, xlab = 'Total CN intensity', 
       ylab = 'Allelic Imbalance 2*|BAF - 0.5|', 
       main = paste('HAPSEG',obj$h, sep = ''))
  points(x[mat$a1<0], y[mat$a1<0], col = col.transp('red', .5))
  title(main = info$sample.name[i], outer = T, line = 0)
  if (!is.null(printfolder)) {  
    dir.create(printfolder, recursive = T, showWarnings = FALSE)
    dev.print(png, 
              file=file.path(printfolder, paste(info$sample.name[i],'_1.1_HAPSEG_TAPS.png', 
                                                sep = '')),  width=640, height=640)
  }
  obj
}

##################################################################
##################################################################
##
## General
##
##################################################################
##################################################################


col.transp <- function(col1, alpha=0.5){
  out = NULL
  for (j in 1:length(col1)){
    tmp=as.list(t(col2rgb(col1[j])))
    tmp$alpha=alpha*255
    tmp$maxColorValue=255
    out = c(out, do.call('rgb',tmp))
  }
  out
}
colors.by = function(x, pal = brewer.pal(9,'YlOrRd'), span = range(na.omit(x)),
                     alpha = .5){
  
  cols = colorRampPalette(pal)(100)
  colbreaks = seq(min(span), max(span), length = length(cols))
  delta = colbreaks[2] - colbreaks[1]
  colind = pmax((x-colbreaks[1]) %/% delta + 1, 1)
  colind[colind>length(cols)] = length(cols)
  mycols = col.transp(cols[colind], alpha = alpha)
  list(mycols=mycols, cols = cols)
}

dist = function(x, y){
  if (is.null(dim(x))) x = matrix(x, ncol = 2)
  if (is.null(dim(y))) y = matrix(y, ncol = 2)
  sqrt((x[,1] - y[,1])**2 + (x[,2] - y[,2])**2)
}



##################################################################
##################################################################









##################################################################
##################################################################
##
## Rotation of CNs
##
##################################################################
##################################################################
startRotation = function(obj){
  i = obj$i
  d = file.path(obj$d, 
                paste(obj$info$sample.name[i], '_2_Start_values', sep = ''))
  mystart = start.alphamax.f(obj$mat, printfile = d, colbychrom = T)
  alpha = mystart$alphamax
  F = mystart$F
  plot.transformed(obj$mat, alpha, 
                   F, xlim = NULL, ylim = NULL)
  d = file.path(obj$d, 
                paste(info$sample.name[i], '_2_Start_values_4.png', sep = ''))
  dev.print(png, file = d, width = 640, height = 500)
  
  Fm1 = solve(F)
  a1p = obj$mat$a1*Fm1[1,1] + obj$mat$a2 * Fm1[1,2]
  a2p = obj$mat$a1*Fm1[2,1] + obj$mat$a2 * Fm1[2,2]
  
  xlim = c(0, quantile(a2p, probs = .95, na.rm = T)*1.2)
  ylim = c(0, quantile(a1p, probs = .95, na.rm = T)*1.2)
  plot.transformed(obj$mat, alpha, 
                   F, xlim = xlim, ylim = ylim)
  d = file.path(obj$d, 
                paste(info$sample.name[i], '_2_Start_values_5.png', sep = ''))
  dev.print(png, file = d, width = 640, height = 500) 
  mystart
}

start.alphamax.f = function(mat, maxCN = 24, maxlines = 30,
                            colbychrom = FALSE, printfile = FALSE){
  #  par(mfrow = c(1,2))
  #x = mat$a1 + mat$a2
  #y = 2*(mat$a2/(x)-.5)
  #xlim = c(quantile(x, probs = c(.025,.975), na.rm = T))
  #ylim = c(quantile(y, probs = c(0, 1), na.rm = T))
  #plot(x, y, pch = 16, col = col.transp('grey', .3), 
  #     xlim = xlim, ylim = ylim, xlab = 'Total CN intensity', 
  #     ylab = 'Allelic Imbalance 2*|BAF - 0.5|', main = 'TAPS')
  xlim = c(0, quantile(mat$a2, probs = c(.9)))
  ylim = c(0, quantile(mat$a1, probs = c(.975)))
  cols = col.transp('grey', .3)
  if (colbychrom) cols = col.transp(as.numeric(mat$Chromosome), .3)
  plot(mat$a2, mat$a1, pch = 16, col = cols, 
       xlim = xlim, ylim = ylim, xlab = 'Major intensity ratio a2', 
       ylab = 'Minor intensity ratio a1',
       main = 'Minor and major CN intensity')
  abline(0,1, col = 'red')
  if (printfile != FALSE) dev.print(png, file = paste(printfile, '_1.png', 
                                                      sep = ''), width = 640, height = 500)
  
  print('Mark two subsequent CN clusters along a "vertical" line:')
  posy = locator(2)
  print('Mark two distant points in one "vertical" line:')
  posl = locator(2)
  print('Mark two subsequent CN clusters along a horizontal line:')
  posx = locator(2)
  print('Mark two distant points in one "horizontal" line:')
  posk = locator(2)
  print('Mark the allelic balance cluster with maximal CN:')
  posb = locator(1)
  if (!is.null(posy)) dy = abs(diff(posy$y))
  if (!is.null(posx)) dx = abs(diff(posx$x))
  l = lm(x~y, data = posl)$coef[2]
  k = lm(y~x, data = posk)$coef[2]
  lines(posy$x, posy$y)
  lines(posl$x, posl$y, lty = 'dotted')
  lines(posx$x, posx$y)
  lines(posk$x, posk$y, lty = 'dotted')
  points(posb$x, posb$y, pch = 4)
  dx = (dx+dy)/2
  dy = dx
  f11 = dy#(dx+dy)/2
  f22 = dx#(dx+dy)/2
  if (printfile != FALSE) dev.print(png, file = paste(printfile, '_2.png', 
                                                      sep = ''), width = 640, height = 500)
  
  F = matrix(c(f11, l*f11, k*f22, f22), ncol = 2, nrow = 2)
  Fm1 = solve(F)
  
  ###Boundaries for F values
  LBf = c(0.5*f11, -Inf, -Inf, 0.5*f22)
  UBf = c(1.5*f11, Inf, Inf, 1.5*f22)
  #UBf[2]
  ix = which.max(posl$y)
  postmp = posl
  postmp$x[ix] = postmp$x[ix]+dx
  UBf[2] = lm(x~y, data = postmp)$coef[2]*UBf[1]
  if (l<0) UBf[2] = min(0, UBf[2])
  #UBf[3]
  ix = which.max(posk$x)
  postmp = posk
  postmp$y[ix] = postmp$y[ix]+dy
  UBf[3] = lm(y~x, data = postmp)$coef[2]*UBf[4]
  if (k<0) UBf[3] = min(0, UBf[3])
  #LBf[2]
  ix = which.max(posl$y)
  postmp = posl
  postmp$x[ix] = postmp$x[ix]-dx
  LBf[2] = lm(x~y, data = postmp)$coef[2]*LBf[1]  
  if (l>0) LBf[2] = max(0, LBf[2])
  #LBf[3]
  ix = which.max(posk$x)
  postmp = posk
  postmp$y[ix] = postmp$y[ix]-dy
  LBf[3] = lm(y~x, data = postmp)$coef[2]*LBf[4]  
  if (k>0) LBf[3] = max(0, LBf[3])
  
  ### Starting value positions of lattice points in E scale
  E1 = Fm1[1,1]*posy$y + Fm1[1,2]*posy$x
  deltay = max(E1) - min(E1)
  b = seq(max(E1), by = -deltay, to = 0)
  by = b[length(b)]
  ny = length(b)
  E1 = seq(by, by = deltay, length.out = maxlines)
  
  E2 = Fm1[2,1]*posx$y + Fm1[2,2]*posx$x
  deltax = max(E2) - min(E2)
  b = seq(max(E2), by = -deltax, to = 0)
  bx = b[length(b)]
  nx = length(b)
  E2 = seq(bx, by = deltax, length.out = maxlines)
  
  E1 = rep(E1, maxlines)
  E2 = rep(E2, each = maxlines)
  
  ### Starting value positions of lattice points in observed scale
  ### Before adjusted to common alpha
  e1 = F[1,1] * E1 + F[1,2]*E2
  e2 = F[2,1] * E1 + F[2,2]*E2
  points(e2, e1)
  
  ### Estimate startvalue of maximum alpha which is <1
  balance.x = Fm1[2,1]*posb$y + Fm1[2,2]*posb$x
  balance.y = Fm1[1,1]*posb$y + Fm1[1,2]*posb$x
  dists = dist(x = c(balance.y, balance.x), y = cbind(E1, E2))
  closest = which.min(dists)
  bbx = seq(E2[closest], by = -deltax, to = 0)
  bby = seq(E1[closest], by = -deltay, to = 0)
  n = min(length(bbx), length(bby))
  nxs = which(bbx>=0)[(length(bbx)-n+1):length(bbx)] 
  nys = which(bby>=0)[(length(bby)-n+1):length(bby)] 
  
  alphax = (1)/(bx+nxs)
  alphay = 1/(by+nys)
  print(paste('Purity estimates 1:', alphay[1]))
  print(paste('Purity estimates 2:', alphax[1]))
  alphamax = mean(c(alphax[1], alphay[1]))
  
  ### Re-estimate starting value positions of Es and F from alpha
  f11 = dy/alphamax#(dx+dy)/2/alphamax#
  f22 = dx/alphamax#(dx+dy)/2/alphamax#
  F = matrix(c(f11, l*f11, k*f22, f22), ncol = 2, nrow = 2)
  Fm1 = solve(F)
  
  A2 = 0:maxCN
  A1 = matrix(A2, ncol = maxCN + 1, nrow = maxCN + 1)
  A2 = t(A1)
  
  E1 = c(alphamax*A1 + (1-alphamax))
  E2 = c(alphamax*A2 + (1-alphamax))
  
  ### Starting value positions of lattice points in observed scale
  e1 = F[1,1] * E1 + F[1,2]*E2
  e2 = F[2,1] * E1 + F[2,2]*E2
  #  points(e2, e1, col = 'blue')
  
  ### Plot starting values on transformed scale
  par(mfrow = c(1,1))
  x = cbind(mat$a1, mat$a2)
  a1 = Fm1[1,1] * x[,1] + Fm1[1,2]*x[,2]
  a2 = Fm1[2,1] * x[,1] + Fm1[2,2]*x[,2]
  xlim = c(0, quantile(a2, probs = c(1)))
  ylim = c(0, quantile(a1, probs = c(1)))
  plot(a2, a1, pch = 16, col = col.transp('grey', .3), 
       xlab = 'Rotated a2', ylab = 'Rotated a1',
       main = 'Minor and major CN intensity',
       xlim = xlim, ylim = ylim)
  lines(range(E2), range(E1), col = 'red')
  for (j in 1:maxCN){
    lines(range(E2), rep(E1[j], 2), lty = 'dotted', col = 'grey20')
    lines(rep(unique(E2)[j], 2), range(E1), lty = 'dotted', col = 'grey20')
  }
  if (printfile != FALSE) dev.print(png, file = paste(printfile, '_3.png', 
                    sep = ''), width = 640, height = 500)
  
  #Adjust F if allelic balance clusters are not on diagonal:
  #  print('Mark the allelic balance point with maximal CN or Esc:')
  #  posb = locator(1)
  #  if (!is.null(posb)){
  #    Fm1[2,] = Fm1[2,]*posb$y/posb$x
  #    F = solve(Fm1)
  #  }
  #Make sure to stay within boundaries:
  #F = pmin(F, UBf)
  #F = pmax(F, LBf)
  #Fm1 = solve(F)
  
  list(alphamax = alphamax, F=F, LBf = LBf, UBf = UBf)
}


rho.alphamax.f = function(pars, x, w, maxCN=24, rhofunc = 'tukey', 
                          rhoargs = NULL){
  alphamax = pars[1]
  F = matrix(pars[-1], ncol = 2, nrow = 2)
  
  #Position of latice points on observed scale
  A2 = 0:maxCN
  A1 = matrix(A2, ncol = maxCN + 1, nrow = maxCN + 1)
  A2 = t(A1)
  E1 = c(alphamax*A1 + (1-alphamax))
  E2 = c(alphamax*A2 + (1-alphamax))
  e1 = F[1,1] * E1 + F[1,2]*E2
  e2 = F[2,1] * E1 + F[2,2]*E2
  
  #Assign observations to closest lattice point
  alldist = apply(x, 1, dist, y = cbind(e1, e2))
  closest = apply(alldist, 2, which.min)
  
  #Calculate residuals
  e = cbind(e1, e2)[closest,]
  r = dist(x, e)*w
  sum(do.call(rhofunc, args = list(x = r, args = rhoargs)))
}

setWeights = function(obj, weight.constant = 500000, rhofunc = 'tukey',
                      nmad = 3){
  obj$w = 1-exp(-obj$mat$length/weight.constant)
  obj$rhofunc = rhofunc
  obj$nmad = nmad
  evaluate.weights(obj$mat, obj$mystart$alphamax, obj$mystart$F, 
                   nmad = nmad, w=obj$w)
  title(paste('w = 1-exp(-segment length/', weight.constant, ')', 
              sep = ''), adj = 1, outer = T, line = -2, font.main = 1)
  dev.print(png, file=file.path(obj$d, paste(obj$info$sample.name[i],'_3.1_evaluate_weights.png', 
                                             sep = '')),  width=900, height=900)
  obj$weight.constant = weight.constant
  print("Change weight constant until 1'000'000 length segs have almost top weight")
  print("Change nmad until Contrib to rho  gets constant at good Distance to lattice point")
  obj
}


evaluate.weights = function(mat, alpha, F, rhofunc = 'tukey', 
                            nmad = 1, maxCN = 24, w=1){
  #Recalculate transformed lattice points from F and alpha
  x = cbind(mat$a1, mat$a2)
  Fm1 = solve(F)  
  A2 = 0:maxCN
  A1 = matrix(A2, ncol = maxCN + 1, nrow = maxCN + 1)
  A2 = t(A1)
  
  E1 = c(alpha*A1 + (1-alpha))
  E2 = c(alpha*A2 + (1-alpha))
  
  ### Positions of lattice points in observed scale and
  ### each segment's distance to it's lattice point -> weight
  e1 = F[1,1] * E1 + F[1,2]*E2
  e2 = F[2,1] * E1 + F[2,2]*E2
  alldist = apply(x, 1, dist, y = cbind(e1, e2))
  closest = apply(alldist, 2, which.min)
  e = cbind(e1, e2)[closest,]
  
  di = dist(x, e)
  r = di*w
  rhoargs = list(k = nmad*mad(di), c = nmad*mad(di))
  eachrho = do.call(rhofunc, args = list(x = r, args = rhoargs))
  
  #Plot
  par(mfrow=c(2,2), mar = c(4, 4, 3, .1))
  plot(r, eachrho, log = 'xy', xlab = 'r = Distance * Weight', 
       ylab = 'Contribution to rho')
  plot(di, eachrho, log = 'xy', xlab = 'Distance to lattice point', 
       ylab = 'Contribution to rho')
  plot(mat$length, w, log = 'xy', xlab = 'Segment length', 
       ylab = 'Weight')
  plot(di, r, log = 'xy', ylab = 'r = Distance * Weight', 
       xlab = 'Distance to lattice point')
  title(main = 'Evaluation of rho function and weights', line = -1,
        outer = T)
}
rho = function(pars, x, w, EE, rhofunc = 'huber', rhoargs = NULL){
  F = matrix(pars, ncol = 2, nrow = 2)
  ee1 = F[1,1] * EE[,1] + F[1,2]*EE[,2]
  ee2 = F[2,1] * EE[,1] + F[2,2]*EE[,2]
  
  #Assign observations to closest lattice point
  alldist = apply(x, 1, dist, y = cbind(ee1, ee2))
  closest = apply(alldist, 2, which.min)
  
  #Calculate residuals
  e = cbind(ee1, ee2)[closest,]
  r = dist(x, e)*w
  sum(do.call(rhofunc, args = list(x = r, args = rhoargs)))
}

huber = function(x, args, ...){
  k = args$k
  out =  x*x/2
  ix = which(abs(x)>k)
  out[ix] = k*(abs(x[ix])-k/2)
  out
}
tukey = function(x, args, ...){
  c = args$c
  out =  c*c/6*(1-(1-(x/c)^2)^3)
  ix = which(abs(x)>c)
  out[ix] = c*c/6
  out
}
leastsq = function(x, ...){
  x^2
}


getdist = function(mat, alpha, F, maxCN = 24){
  x = cbind(mat$a1, mat$a2)
  Fm1 = solve(F)  
  A2 = 0:maxCN
  A1 = matrix(A2, ncol = maxCN + 1, nrow = maxCN + 1)
  A2 = t(A1)
  
  E1 = c(alpha*A1 + (1-alpha))
  E2 = c(alpha*A2 + (1-alpha))
  
  ### Positions of lattice points in observed scale and
  ### each segment's distance to it's lattice point -> weight
  e1 = F[1,1] * E1 + F[1,2]*E2
  e2 = F[2,1] * E1 + F[2,2]*E2
  alldist = apply(x, 1, dist, y = cbind(e1, e2))
  closest = apply(alldist, 2, which.min)
  e = cbind(e1, e2)[closest,]
  
  di = dist(x, e)
  di
}

optimizeRotation = function(obj){
  mat = obj$mat
  mystart = obj$mystart
  di = getdist(mat, mystart$alphamax, mystart$F) 
  obj$rhoargs = list(c = mad(di)*obj$nmad)
  x = cbind(mat$a1, mat$a2)
  
  pars = solnp(pars = c(mystart$alphamax, mystart$F), fun = rho.alphamax.f, x = x,
               w = obj$w, maxCN = 24, rhofunc = obj$rhofunc, 
               rhoargs = obj$rhoargs)$pars
  
  alpha = pars[1]
  F = matrix(pars[-1], ncol = 2, nrow = 2)
  Fm1 = solve(F)
  a1p = obj$mat$a1*Fm1[1,1] + obj$mat$a2 * Fm1[1,2]
  a2p = obj$mat$a1*Fm1[2,1] + obj$mat$a2 * Fm1[2,2]
  xlim = c(0, quantile(a2p, probs = .95, na.rm = T)*1.2)
  ylim = c(0, quantile(a1p, probs = .95, na.rm = T)*1.2)
  plot.transformed(mat, alpha, F, w=obj$w, rhoargs=obj$rhoargs, rhofunc = obj$rhofunc, 
                   xlim = xlim, ylim = ylim)
  dev.print(png, file=file.path(obj$d, paste(info$sample.name[i],'_3.2_rotated_grid.png', 
                                             sep = '')),  width=640, height=500)
  
  #All possible alphas and Fs
  maxCN = 24
  
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
  obj$alphas = alphas
  obj$fs = fs
  obj$rotation.pars = pars
  obj
}

####################
###Old
####################
transform.startval = function(a1, a2, alpha, maxCN=24){
  
  #Standard setup
  A2 = 0:maxCN
  A1 = matrix(A2, ncol = maxCN + 1, nrow = maxCN + 1)
  A2 = t(A1)
  ix.na = A1>A2
  A1[ix.na] = NA
  A2[ix.na] = NA
  
  E1 = alpha*A1 + (1-alpha)
  E2 = alpha*A2 + (1-alpha)
  
  #Plot
  xlim = c(quantile(a2, probs = c(.025,.975)))
  ylim = c(quantile(a1, probs = c(.025,.975)))
  plot(a2, a1, pch = 16, col = col.transp('grey', .3), 
       xlim = xlim, ylim = ylim, xlab = 'a2', ylab = 'a1')
  abline(0,1, col = 'red')
  title(main = i)
  
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
  
  f11 = dy/alpha
  f22 = dx/alpha
  
  F = matrix(c(f11, l*f11, k*f22, f22), ncol = 2, nrow = 2)
  
  ### Starting value positions of lattice points
  e1 = F[1,1] * E1 + F[1,2]*E2
  e2 = F[2,1] * E1 + F[2,2]*E2
  ee1 = c(e1)[!is.na(c(e1))]
  ee2 = c(e2)[!is.na(c(e2))]
  
  ###Assign observations to closest lattice point
  x = cbind(a1, a2)
  alldist = apply(x, 1, dist, y = cbind(ee1, ee2))
  closest = apply(alldist, 2, which.min)
  points(a2, a1, pch = 16, col = col.transp(closest, .3), cex = .7)
  points(ee2, ee1)
  F
}


##################################################################
##################################################################






##################################################################
##################################################################
##
## Read mutations and set rotation/scaling
##
##################################################################
##################################################################


readMutations = function(obj){
  
  ###Read mutations
  maf.fn = paste(getwd(), 
                 '/AROMA/rawData/responsify/Exome/MAF/AllVarintsMutectOncotatorPlusMutsigCategNo5214_mindepth10_withvaf', '_', 
                 obj$info$sample.name[obj$i], '.maf', sep = '')
  
  maf = read.delim(maf.fn, row.names = NULL, stringsAsFactors = FALSE, 
                   check.names = FALSE, na.strings = c("NA", "---"), 
                   blank.lines.skip = TRUE, comment.char = "#")
  maf$depth = maf$t_alt_count + maf$t_ref_count
  
  #Assign mutations to segments
  mat = obj$mat
  maf$segment = NA
  for (seg in 1:nrow(mat)){
    ix = which(maf$Start_position >= mat$Start.bp[seg] & 
                 maf$Start_position <= mat$End.bp[seg]
               & as.character(maf$Chromosome) == as.character(mat$Chromosome[seg]))
    maf$segment[ix] = seg
  } 
  muts = unique(subset(maf, select = c('Chromosome','Start_position','End_position',
                                       't_alt_count','t_ref_count','segment', 'depth',
                                       'effect', 'Hugo_Symbol')))
  print(paste(nrow(muts), 'mutations in sample'))
  muts$af = muts$t_alt_count/(muts$t_alt_count+muts$t_ref_count)
  muts = subset(muts, !is.na(segment) )
  muts = subset(muts, muts$depth>=20)
  print(paste('Keeping', nrow(muts), 'mutations with read depth at least 20') )
  obj$muts = muts
  obj
}

plotAfByCN = function(obj, ns = 1:4){
  mat = obj$mat
  muts = obj$muts
  i = obj$i
  
  par(mfrow = c(2,2))
  for (n in ns){
    F = matrix(obj$rotation.pars[-1], ncol = 2, nrow = 2)
    alpha = obj$alphas[n]
    F = F*obj$fs[n]
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
  dev.print(png, file=file.path(obj$d, paste(obj$info$sample.name[i],'_4.1_minor_major_af.png', 
                                         sep = '')),  width=900, height=900)
}


plotAfByPosition = function(obj, ns = 2:3, restrict = FALSE){
  ###Plot by genome position
  
  par(mfrow = c(2,1))
  mat = obj$mat
  muts = obj$muts
  i = obj$i
  for (n in ns){
    F = matrix(obj$rotation.pars[-1], ncol = 2, nrow = 2)
    alpha = obj$alphas[n]
    F = F*obj$fs[n]
    Fm1 = solve(F)
    c1 = mat$a1*Fm1[1,1] + mat$a2 * Fm1[1,2]
    c2 = mat$a1*Fm1[2,1] + mat$a2 * Fm1[2,2]
    c1 = (c1 - (1 - alpha))/alpha
    c2 = (c2 - (1 - alpha))/alpha 
    mult = muts$af*(alpha*((c1 + c2)[muts$segment])+(1-alpha)*2)/alpha
    plot.CNmut(plotmat = mat, plotmuts = muts, c1 = c1, c2 = c2, mult = mult,
               main = paste('n =', n, ', purity', round(alpha*100), '%'),
               restrict = restrict)
  }
  dev.print(png, file=file.path(obj$d, paste(obj$info$sample.name[i],'_4.2_CNmuts_preliminary.png', 
                                         sep = '')),  width=900, height=900)  
}

setRotation = function(obj, n){
  mat = obj$mat
  alpha = obj$alphas[n]
  F = matrix(obj$rotation.pars[-1], ncol = 2, nrow = 2)
  F = F*obj$fs[n]
  Fm1 = solve(F)
  mat$a1r = obj$mat$a1*Fm1[1,1] + obj$mat$a2 * Fm1[1,2]
  mat$a2r = obj$mat$a1*Fm1[2,1] + obj$mat$a2 * Fm1[2,2]
  xlim = c(0, quantile(mat$a2r, probs = .95, na.rm = T)*1.2)
  ylim = c(0, quantile(mat$a1r, probs = .95, na.rm = T)*1.2)
  plot.transformed(mat, alpha, F = F, w=obj$w, rhoargs=obj$rhoargs, rhofunc = obj$rhofunc, 
                   xlim = xlim, ylim = ylim)
  mtext(paste('Purity', round(alpha*100), '%'),adj = 0.1, outer = T, line = -2)
  dev.print(png, file=file.path(obj$d, paste(info$sample.name[i],'_4.3_rotated_grid_chosen.png', 
                                             sep = '')),  width=640, height=500)
  #Save settings
  obj$alpha = alpha
  obj$F = F
  obj$Fm1 = Fm1
  obj$n = n
  obj$mat = mat
  obj  
} 

setCNs = function(obj, linearize = FALSE){
  alpha = obj$alpha
  mat = obj$mat
  vari1 = mat$a1r
  vari2 = mat$a2r
  if (linearize) {
    vari1 = mat$a1rl
    vari2 = mat$a2rl
  }
  
  c1 = (vari1 - (1 - alpha))/alpha
  c2 = (vari2 - (1 - alpha))/alpha
  mat$c1 = pmin(c1, c2)
  mat$c2 = pmax(c1, c2)
  obj$ploidy = sum((mat$c1 + mat$c2)*mat$W)
  
  par(mfrow = c(1,1), mar = c(3, 3, 2, .1))
  xlim = c(-1, quantile(mat$c2, probs = c(.95))*c(1.2))
  ylim = c(-1, quantile(mat$c1, probs = c(.95))*c(1.2))
  
  plot.cn(x = mat$c2, y = mat$c1, mycols = col.transp('blue', alpha = .5),
          xlim = xlim, ylim = ylim, xlab = 'Major CN', ylab = 'Minor CN')
  mtext(paste('Ploidy =', format(round(obj$ploi,1),nsmall = 1)), adj = 0)
  dev.print(png, file=file.path(obj$d, 
                                paste(info$sample.name[i],'_5.4_CNs.png', 
                                      sep = '')),  width=400, height=400)
  obj$mat = mat
  obj
}
setMultiplicity = function(obj){
  muts = obj$muts
  #Express vaf as number of copies per cell (cellular multiplicity)
  alpha = obj$alpha
  muts$mult = muts$af/alpha*
    (alpha*((mat$c1 + mat$c2)[muts$segment])+(1-alpha)*2)
  obj$muts = muts
  obj
}

##################################################################
##################################################################











##################################################################
##################################################################
##
## Linearization
##
##################################################################
##################################################################

getpeaks = function(val, alpha, maxCN=10){
  EEs = alpha*c(0:maxCN) + (1-alpha) 
  delta = EEs[2] - EEs[1]
  peaks = numeric(length(EEs))
  expected <- EEs[1]
  
  den = density(val,kernel="gaussian", bw = .03)
  ix = abs(den$x - expected) < (delta/4)
  tp = try((den$x[ix])[which.max(den$y[ix])], silent = T)
  peaks[1] = ifelse(is.na(tp[1]) | inherits(tp, 'try-error'), expected, tp)
  
  for (j in 2:length(EEs)){
    expected = peaks[j-1] + delta
    ix = abs(den$x - expected) < (delta/4)
    tp = try((den$x[ix])[which.max(den$y[ix])], silent = T)
    peaks[j] = ifelse(is.na(tp[1]) | inherits(tp, 'try-error'), expected, tp)
    delta = peaks[j] - peaks[j-1]
  }
  peaks
}

linearize = function(obj, vari = 'a1r', 
                     output = FALSE, logis = FALSE, asymporig = FALSE, none = FALSE, plot = FALSE){
  y = obj$mat[,vari]
  if (none) return(y)
  alpha = obj$alpha
  maxCN = 6
  bw = .03
  if (output & logis & asymporig) return(print('Cannot return both logis and asymporig estimates.'))
  #Find peaks
  EEs = alpha*c(0:maxCN) + (1-alpha) 
  ypeaks = getpeaks(y, alpha = alpha, maxCN = maxCN)
  par(mfrow = c(2,1), mar = c(3,3,1,1))
  yden = density(y,kernel="gaussian", bw = bw)
  yden$y = log(yden$y+1)
  plot(1, 1, xlim = c(0, max(EEs)), ylim = range(yden$y), 
       type = 'n', yaxt = 'n')
  title(ylab = 'Smoothed frequency', xlab = 'Array intensity', line = 2)
  legend('topleft', lty = c('solid', 'dotted'), col = c(1,2), 
         legend = c('Suggested peak positons','Linear peak positions'),
         bty = 'n', cex = .8)
  polygon(yden$x, yden$y, col = col.transp(1, .5), border = NA)
  abline(v = ypeaks)
  abline(v = EEs, col = 'red', lty = 'dotted')
  
  plot(EEs, ypeaks, xlim = c(0, max(EEs)), ylim = c(0, max(EEs)))
  title(xlab = 'Linear peak position', ylab = 'Suggested peak position', line = 2)
  abline(0,1, col = 'red')
  df = data.frame(x = EEs, y = ypeaks)
  
  #Logistic regression
  leg = 'No linear adjustment'
  cols = 2
  if (asymporig){
    try({
      out = NULL
      leg = c(leg, 'Asymptotic regression through origin')
      cols = c(cols, 4)
      mod <- nls(ypeaks ~ SSasympOrig(EEs, Asym, lrc), data = df)
      g = function(x) coef(mod)['Asym']*(1-exp(-x*exp(coef(mod)['lrc'])))
      lines(EEs, g(EEs), col = 'blue')
      gm1 = function(y, mod) {
        y = pmin(coef(mod)['Asym']*.9, y)
        - exp(-coef(mod)['lrc'])*log(1-y/coef(mod)['Asym'])
        
        breakat = .9
        ybreak = coef(mod)['Asym']*breakat
        xbreak = - exp(-coef(mod)['lrc'])*log(1-ybreak/coef(mod)['Asym'])
        ix = y>ybreak
        out = numeric(length(y))
        out[!ix] = - exp(-coef(mod)['lrc'])*log(1-y[!ix]/coef(mod)['Asym'])
        yvec = ybreak*.99
        xvec = - exp(-coef(mod)['lrc'])*log(1-yvec/coef(mod)['Asym'])
        yvec = c(yvec, ybreak)
        xvec = c(xvec, xbreak)
        linmod = lm(xvec~yvec)
        out[ix] = predict(linmod, newdata = data.frame(yvec = y[ix]))
        out
      }
      out = gm1(y, mod)
      if (is.null(out)) leg[length(leg)] = paste(leg[length(leg)], 'not converging')
    }, silent = T)
  }
  if (logis){
    try({
      out = NULL
      leg = c(leg, 'Logistic regression')
      cols = c(cols, 3)
      mod <- nls(ypeaks ~ SSlogis(EEs, Asym, xmid, scal), data = df)
      
      g = function(x) coef(mod)['Asym']/(1+exp((coef(mod)['xmid']-x)/coef(mod)['scal']))
      lines(EEs, g(EEs), col = 'green')
      gm1 = function(y, mod) {
        breakat = .9
        ybreak = coef(mod)['Asym']*breakat
        xbreak = coef(mod)['xmid'] - coef(mod)['scal']*
          log((coef(mod)['Asym']-ybreak)/
                ybreak)
        ix = y>ybreak
        out = numeric(length(y))
        out[!ix] = coef(mod)['xmid'] - coef(mod)['scal']*
          log((coef(mod)['Asym']-y[!ix])/y[!ix])
        yvec = ybreak*.99
        xvec = coef(mod)['xmid'] - coef(mod)['scal']*
          log((coef(mod)['Asym']-yvec)/yvec)
        yvec = c(yvec, ybreak)
        xvec = c(xvec, xbreak)
        linmod = lm(xvec~yvec)
        out[ix] = predict(linmod, newdata = data.frame(yvec = y[ix]))
        out
      }
      out = gm1(y, mod)
      if (is.null(out)) leg[length(leg)] = paste(leg[length(leg)], 'not converging')
    }, silent = T)
  }
  legend('topleft', lty = rep('solid', length(leg)), col = c(2,4,3), 
         legend = leg,
         bty = 'n', cex = .8)  
  if (length(grep('1', vari))>0 & plot) dev.print(png, file=file.path(obj$d, 
                                                                      paste(obj$info$sample.name[obj$i],'_5.1_nonlinearlity.png', 
                                                                            sep = '')),  width=900, height=900)
  if (length(grep('2', vari))>0 & plot) dev.print(png, file=file.path(obj$d, 
                                                                      paste(obj$info$sample.name[obj$i],'_5.2_nonlinearlity.png', 
                                                                            sep = '')),  width=900, height=900)
  if (output) return(out)
  NULL
}

plotLin = function(obj){
  mat = obj$mat
  alpha = obj$alpha
  maxCN = 24
  
  #Calculate lattice points from alpha
  A2 = 0:maxCN
  A1 = matrix(A2, ncol = maxCN + 1, nrow = maxCN + 1)
  A2 = t(A1)
  E1 = c(alpha*A1 + (1-alpha))
  E2 = c(alpha*A2 + (1-alpha))
  EEs = unique(E1)
  
  #Plot CN lattice
  par(mfrow = c(2,2))
  xlim = c(0, quantile(mat$a2r, probs = .95)*1.2)
  ylim = c(0, quantile(mat$a1r, probs = .95)*1.2)
  plot.cn(x = mat$a2r, y = mat$a1r, 
          xlim = xlim, ylim = ylim, grid = FALSE,
          xlab = 'Rotated a2',
          ylab = 'Rotated a1')
  for (j in 1:length(EEs)){
    lines(range(EEs), rep(EEs[j], 2), lty = 'dotted', col = 'grey20')
    lines(rep(EEs[j], 2), range(EEs), lty = 'dotted', col = 'grey20')
  }
  plot.cn(x = mat$a2rl, y = mat$a1rl, 
          xlim = xlim, ylim = ylim, grid = FALSE,
          xlab = 'Rotated a2 after non-linearity fix',
          ylab = 'Rotated a1 after non-linearity fix')
  for (j in 1:length(EEs)){
    lines(range(EEs), rep(EEs[j], 2), lty = 'dotted', col = 'grey20')
    lines(rep(EEs[j], 2), range(EEs), lty = 'dotted', col = 'grey20')
  }
  
  ### OVERLAY PLOT Without reattenuation
  ### Each segment's distance to it's lattice point
  x = cbind(mat$a1r, mat$a2r)
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
  x = cbind(mat$a1rl, mat$a2rl)
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
  dev.print(png, file=file.path(obj$d, paste(info$sample.name[i],'_5.3_linearization.png', 
                                             sep = '')),  width=900, height=900)  
}


nonlinAdj.asympOrig = function(y, alpha = 1, maxCN = 10, bw = .03){
  EEs = alpha*c(0:maxCN) + (1-alpha) 
  ypeaks = getpeaks(y, alpha = alpha, maxCN = maxCN)
  par(mfrow = c(2,1), mar = c(3,3,1,1))
  #Single scale density y
  yden = density(y,kernel="gaussian", bw = bw)
  yden$y = log(yden$y+1)
  plot(1, 1, xlim = c(0, max(EEs)), ylim = range(yden$y), 
       type = 'n', yaxt = 'n')
  title(ylab = 'Smoothed frequency', xlab = 'Array intensity', line = 2)
  legend('topleft', lty = c('solid', 'dotted'), col = c(1,2), 
         legend = c('Suggested peak positons','Linear peak positions'),
         bty = 'n', cex = .8)
  polygon(yden$x, yden$y, col = col.transp(1, .5), border = NA)
  abline(v = ypeaks)
  abline(v = EEs, col = 'red', lty = 'dotted')
  
  plot(EEs, ypeaks, xlim = c(0, max(EEs)), ylim = c(0, max(EEs)))
  title(xlab = 'Linear peak position', ylab = 'Suggested peak position', line = 2)
  abline(0,1, col = 'red')
  df = data.frame(x = EEs, y = ypeaks)
  mod <- nls(ypeaks ~ SSasympOrig(EEs, Asym, lrc), data = df)
  g = function(x) coef(mod)['Asym']*(1-exp(-x*exp(coef(mod)['lrc'])))
  lines(EEs, g(EEs), col = 'blue')
  legend('topleft', lty = c('solid', 'solid'), col = c(2,4), 
         legend = c('Diagonal through (0,0)','Fitted trend'),
         bty = 'n', cex = .8)  
  gm1 = function(y, mod) {
    y = pmin(coef(mod)['Asym']*.9, y)
    - exp(-coef(mod)['lrc'])*log(1-y/coef(mod)['Asym'])
    
    breakat = .9
    ybreak = coef(mod)['Asym']*breakat
    xbreak = - exp(-coef(mod)['lrc'])*log(1-ybreak/coef(mod)['Asym'])
    ix = y>ybreak
    out = numeric(length(y))
    out[!ix] = - exp(-coef(mod)['lrc'])*log(1-y[!ix]/coef(mod)['Asym'])
    yvec = ybreak*.99
    xvec = - exp(-coef(mod)['lrc'])*log(1-yvec/coef(mod)['Asym'])
    yvec = c(yvec, ybreak)
    xvec = c(xvec, xbreak)
    linmod = lm(xvec~yvec)
    out[ix] = predict(linmod, newdata = data.frame(yvec = y[ix]))
    out
  }
  par(mfrow = c(1,1))
  gm1(y, mod)
}

nonlinAdj.logis = function(y, alpha = 1, maxCN = 10, bw = .03){
  EEs = alpha*c(0:maxCN) + (1-alpha) 
  ypeaks = getpeaks(y, alpha = alpha, maxCN = maxCN)
  par(mfrow = c(2,1), mar = c(3,3,1,1))
  #Single scale density y
  yden = density(y,kernel="gaussian", bw = bw)
  yden$y = log(yden$y+1)
  plot(1, 1, xlim = c(0, max(EEs)), ylim = range(yden$y), 
       type = 'n', yaxt = 'n')
  title(ylab = 'Smoothed frequency', xlab = 'Array intensity', line = 2)
  legend('topleft', lty = c('solid', 'dotted'), col = c(1,2), 
         legend = c('Suggested peak positons','Linear peak positions'),
         bty = 'n', cex = .8)
  polygon(yden$x, yden$y, col = col.transp(1, .5), border = NA)
  abline(v = ypeaks)
  abline(v = EEs, col = 'red', lty = 'dotted')
  
  plot(EEs, ypeaks, xlim = c(0, max(EEs)), ylim = c(0, max(EEs)))
  title(xlab = 'Linear peak position', ylab = 'Suggested peak position', line = 2)
  abline(0,1, col = 'red')
  df = data.frame(x = EEs, y = ypeaks)
  mod <- nls(ypeaks ~ SSlogis(EEs, Asym, xmid, scal), data = df)
  
  g = function(x) coef(mod)['Asym']/(1+exp((coef(mod)['xmid']-x)/coef(mod)['scal']))
  lines(EEs, g(EEs), col = 'blue')
  legend('topleft', lty = c('solid', 'solid'), col = c(2,4), 
         legend = c('Diagonal through (0,0)','Fitted trend'),
         bty = 'n', cex = .8)  
  gm1 = function(y, mod) {
    breakat = .9
    ybreak = coef(mod)['Asym']*breakat
    xbreak = coef(mod)['xmid'] - coef(mod)['scal']*
      log((coef(mod)['Asym']-ybreak)/
            ybreak)
    ix = y>ybreak
    out = numeric(length(y))
    out[!ix] = coef(mod)['xmid'] - coef(mod)['scal']*
      log((coef(mod)['Asym']-y[!ix])/y[!ix])
    yvec = ybreak*.99
    xvec = coef(mod)['xmid'] - coef(mod)['scal']*
      log((coef(mod)['Asym']-yvec)/yvec)
    yvec = c(yvec, ybreak)
    xvec = c(xvec, xbreak)
    linmod = lm(xvec~yvec)
    out[ix] = predict(linmod, newdata = data.frame(yvec = y[ix]))
    out
  }
  par(mfrow = c(1,1))
  gm1(y, mod)
}
##################################################################
##################################################################








##################################################################
##################################################################
##
## CN lattice plots
##
##################################################################
##################################################################
plot.transformed = function(mat, alpha, F, w=NULL, rhofunc = 'leastsq',
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
  r = dist(x, e)*w
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
  
  par(mfrow = c(1,1))
  split.screen(rbind(c(0.2, 0.98, 0.2, 0.9),
                     c(0.1,0.2,0.2, 0.9), 
                     c(0.2, 0.98, 0.1, 0.2),
                     c(0.1, 0.2, 0.1, 0.2)))
  screen(1)
  par(mar = c(0,0,0,0))
  
  
  #Plot frame
  if (is.null(xlim)) xlim = c(0, max(c(a2p, E2), na.rm = T))
  if (is.null(ylim)) ylim = c(0, max(a1p, na.rm = T)*1.2)
  plot(a2p, a1p, pch = 16, 
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
  
  #Single scale density y
  screen(2)
  par(mar = c(0,0,0,0))
  yden = density(a1p,kernel="gaussian", bw = .03)
  roughy = round(sum(diff(yden$y)**2)/diff(yden$x[1:2]), 3)
  yden$y = log(yden$y+1)
  plot(1, 1, ylim = ylim, xlim = range(yden$y), type = 'n', xaxt = 'n')
  polygon(yden$y, yden$x, col = col.transp(1, .5), border = NA)
  title(ylab = 'Rotated a1', line = -1, outer = T)
  
  #Single scale density x
  screen(3)
  par(mar = c(0,0,0,0))
  xden = density(a2p,kernel="gaussian", bw = .03)
  roughx = round(sum(diff(xden$y)**2)/diff(xden$x[1:2]), 3)
  xden$y = log(xden$y+1)
  plot(1, 1, xlim = xlim, ylim = range(xden$y), type = 'n', yaxt = 'n')
  polygon(xden$x, xden$y, col = col.transp(1, .5), border = NA)
  title(xlab = 'Rotated a2', line = -1, outer = T)
  
  screen(4)
  par(mar = c(0,0,0,0))
  plot(1,1, type = 'n', xaxt = 'n', yaxt = 'n', xlim = c(0,1), ylim = c(0,1), bty = 'n')
  text(0, 1, roughx, adj = c(0,1))
  text(1, 0, roughy, adj = c(1,0))
  
  title(main = main, outer = T, line = -1)
  close.screen(all.screens = TRUE)
  c(roughx = roughx, roughy=roughy)
  
}

plot.cn = function(x, y, main = NULL, xlim = NULL, xlab = NULL, ylab = NULL,
                   ylim = NULL, maxCN = 24, mycols = NULL, grid = T){
  E1 = E2 = 0:maxCN
  
  #Plot frame
  if (is.null(xlim)) xlim = c(0, max(c(x, E2), na.rm = T))
  if (is.null(ylim)) ylim = c(0, max(y, na.rm = T)*1.2)
  type = ifelse(is.null(mycols), 'p', 'n')
  plot(x, y, pch = 16, col = col.transp('blue', .3), type = type,
       xlim = xlim, ylim = ylim, xlab = '', ylab = '')
  if (is.vector(mycols)){
    points(x, y, pch = 16, col = mycols, cex = 1)    
  }
  if (is.list(mycols)) {
    points(x, y, pch = 16, col = mycols$mycols, cex = 1)
    
    #Colour scale
    tmp = c(xlim[2]/2, xlim[2])
    tmp = seq(tmp[1], tmp[2], length=length(mycols$cols))
    points(tmp, rep(ylim[1], length(tmp)), col = mycols$cols, pch = 16)
    
  }
  #Gridlines
  if (grid){
    EEs = unique(E1)
    for (j in 1:length(EEs)){
      lines(range(EEs), rep(EEs[j], 2), lty = 'dotted', col = 'grey20')
      lines(rep(EEs[j], 2), range(EEs), lty = 'dotted', col = 'grey20')
    }
  }
  abline(0,1)
  title(main = main, outer = T, line = -1)
  title(xlab = xlab, ylab = ylab, line = 2)
}


####################
###Devel
####################

plot.transformed.2D = function(mat, alpha, F, w=NULL, 
                               rhoargs = NULL, main = NULL, xlim = NULL, ylim = NULL,
                               maxCN = 24, color.by = 'eachrho'){
  require('RColorBrewer')
  require('graphics')
  
  # #Rotate first (observed points' scale did not work)
  Fm1 = solve(F)
  c1 = mat$a1*Fm1[1,1] + mat$a2 * Fm1[1,2]
  c2 = mat$a1*Fm1[2,1] + mat$a2 * Fm1[2,2]
  
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
  
  c = cov.trob(X, center = c(0,0), nu = rhoargs$nu)$cov 
  eachweight = mahalanobis(X, center = c(0,0), cov = c)
  
  #Colour palette
  #display.brewer.all()
  pal = brewer.pal(9, 'YlOrRd')
  #display.brewer.pal(8, 'Set1')
  cols = pal[-c(1:2)]
  cols = colorRampPalette(cols)(100)
  
  #Transformed positions of data
  a1p = c1
  a2p = c2
  
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


####################
###Old
####################

plot.transformed.old = function(a1, a2, w, Fout, EE, rhofunc='leastsq', 
                                rhoargs=list(k=1, c=1), color.by = 'eachweight', 
                                xlim = NULL, ylim = NULL){
  require('RColorBrewer')
  #Recalculate optimized fit weights
  x = cbind(a1, a2)
  ee1 = Fout[1,1] * EE[,1] + Fout[1,2]*EE[,2]
  ee2 = Fout[2,1] * EE[,1] + Fout[2,2]*EE[,2]
  alldist = apply(x, 1, dist, y = cbind(ee1, ee2))
  closest = apply(alldist, 2, which.min)
  e = cbind(ee1, ee2)[closest,]
  r = dist(x, e)/w
  eachweight = do.call(rhofunc, args = list(x = r, args = rhoargs))
  
  #Transform
  Fm1 = solve(Fout)
  a1p = a1*Fm1[1,1] + a2 * Fm1[1,2]
  a2p = a1*Fm1[2,1] + a2 * Fm1[2,2]
  
  #Colour palette
  #display.brewer.all()
  pal = brewer.pal(9, 'YlOrRd')
  #display.brewer.pal(8, 'Set1')
  cols = pal[-c(1:2)]
  cols = colorRampPalette(cols)(100)
  
  #Plot frame
  if (is.null(xlim)) xlim = c(0, max(c(a2p, EE[,2]), na.rm = T))
  if (is.null(ylim)) ylim = c(0, max(a1p, na.rm = T)*1.2)
  plot(a2p, a1p, pch = 16, col = col.transp('grey', .3), 
       xlim = xlim, ylim = ylim, xlab = 'Transformed a2', 
       ylab = 'Transformed a1', type = 'n')
  title(main = paste('Sample', i, rhofunc, 'minfunction'))
  
  #Single scale densities
  EEs = unique(EE[,1])
  xden = density(a2p,kernel="gaussian", bw = .03)
  rough = round(sum(diff(xden$y)**2)/diff(xden$x[1:2]), 3)
  xden$y = log(xden$y+1)
  xden$y = xden$y/max(xden$y)*EEs[1]*.9
  polygon(xden$x, xden$y, col = col.transp(1, .5), border = NA)
  text(sum(range(axTicks(1))*c(.9,.1)), 0, rough, pos = 1, offset = 0)
  
  yden = density(a1p,kernel="gaussian", bw = .03)
  rough = round(sum(diff(xden$y)**2)/diff(xden$x[1:2]), 3)
  yden$y = log(yden$y+1)
  yden$y = yden$y/max(yden$y)*EEs[1]*.9
  polygon(yden$y, yden$x, col = col.transp(1, .5), border = NA)
  text(0, sum(range(axTicks(2))*c(.9,.1)), rough, srt = 90, pos = 3, offset = 0)
  
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
  
  #Colour scale
  if (color.by != 'nocol'){
    tmp = c(xlim[2]/2, xlim[2])
    tmp = seq(tmp[1], tmp[2], length=length(cols))
    points(tmp, rep(-.05, length(tmp)), col = cols, pch = 16)
  }
  
}


plot.lattice = function(x, y, alpha, main = NULL, xlim = NULL, 
                        ylim = NULL, maxCN = 24, colby = 'r'){
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
  
  
  #  dummy = try(close.screen(all.screens = T), silent = T)
  par(mfrow = c(1,1))
  split.screen(rbind(c(0.2, 0.98, 0.2, 0.9),
                     c(0.1,0.2,0.2, 0.9), 
                     c(0.2, 0.98, 0.1, 0.2),
                     c(0.1, 0.2, 0.1, 0.2)))
  screen(1)
  par(mar = c(0,0,0,0))
  
  
  #Plot frame
  if (is.null(xlim)) xlim = c(0, max(c(x, E2), na.rm = T))
  if (is.null(ylim)) ylim = c(0, max(y, na.rm = T)*1.2)
  plot(x, y, pch = 16, col = col.transp('grey', .3), type = 'n',
       xlim = xlim, ylim = ylim, xaxt = 'n', yaxt = 'n')
  
  #Points
  
  if (colby[1] == 'r') {
    colby = r
    colbreaks = seq(min(colby), 
                    max(colby[x<max(E2) & x>min(E2) & y>min(E1) & y<max(E1)]), length = length(cols))
  } else colbreaks = seq(min(colby), max(colby), length = length(cols))
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
  
  
  #Single scale density y
  screen(2)
  par(mar = c(0,0,0,0))
  yden = density(y,kernel="gaussian", bw = .03)
  roughy = round(sum(diff(yden$y)**2)/diff(yden$x[1:2]), 3)
  yden$y = log(yden$y+1)
  plot(1, 1, ylim = ylim, xlim = range(yden$y), type = 'n', xaxt = 'n')
  polygon(yden$y, yden$x, col = col.transp(1, .5), border = NA)
  title(ylab = 'Transformed a1', line = -1, outer = T)
  
  #Single scale density x
  screen(3)
  par(mar = c(0,0,0,0))
  xden = density(x,kernel="gaussian", bw = .03)
  roughx = round(sum(diff(xden$y)**2)/diff(xden$x[1:2]), 3)
  xden$y = log(xden$y+1)
  plot(1, 1, xlim = xlim, ylim = range(xden$y), type = 'n', yaxt = 'n')
  polygon(xden$x, xden$y, col = col.transp(1, .5), border = NA)
  title(xlab = 'Transformed a2', line = -1, outer = T)
  
  screen(4)
  par(mar = c(0,0,0,0))
  plot(1,1, type = 'n', xaxt = 'n', yaxt = 'n', xlim = c(0,1), 
       ylim = c(0,1), bty = 'n')
  text(.5, .5, 'Purity', adj = c(0.5,-.2))
  text(.5, .5, paste(round(alpha*100), '%'), adj = c(.5,1.2))
  
  title(main = main, outer = T, line = -1)
  close.screen(all.screens = TRUE)
  
}
##################################################################
##################################################################








##################################################################
##################################################################
##
## Subclonal mutations and CNs
##
##################################################################
##################################################################

mahalanobis.w = function(X, nu, w){
  c = cov.trob(X, center = c(0,0), nu = nu)$cov
  D = numeric(nrow(X))
  for (j in 1:nrow(X)){
    D[j] = mahalanobis(X[j,], center = c(0,0), cov = c/w[j])
  } 
  D
}


subclonalMutations = function(obj){
  alpha = obj$alpha
  mat = obj$mat
  muts = obj$muts
  
  par(mfrow = c(1,1))
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
  muts$subclonal = subclonal.muts(muts = muts, mat = mat, alpha = alpha)
  
  points(muts$index[muts$subclonal], muts$af[muts$subclonal], 
         col = 'darkgreen', pch = 16)
  dev.print(png, file=file.path(obj$d, 
                                paste(obj$info$sample.name[obj$i],'_7.1_minor_major_af.png', 
                                      sep = '')),  width=640, height=400)
  obj$muts = muts
  obj
}

subclonal.muts = function(muts, mat = NULL, c1 = NULL, c2 = NULL, 
                          alpha = NULL, sign.level = .01){
  if(!is.null(c1)) {
    if (is.null(mat)) mat = data.frame(c1 = c1) else mat$c1 = c1
  }
  if(!is.null(c2)) mat$c2 = c2  
  mat$c1 = pmax(0, mat$c1)
  mat$c2 = pmax(0, mat$c2)
  minor = alpha*mat$c1/(alpha*(mat$c1+mat$c2) + 
                              (1-alpha)*2)
  major = alpha*mat$c2/(alpha*(mat$c1+mat$c2) + 
                              (1-alpha)*2)
  both = alpha*(mat$c1+mat$c2)/(alpha*(mat$c1+mat$c2) + 
                                      (1-alpha)*2)
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
  muts$subclonal
}

subclonalCNdist = function(obj, cut = 3){
  mat = obj$mat
  w = obj$w
  
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
  
  par(mfrow = c(2,2), mar = c(4, 4, 3, .1))
  #Multivariate t robust covariance qq=plots
  nu = 2
  D = mahalanobis.w(X, nu, w)  
  cols = (colors.by(D, pal, span = c(0,10))$mycols)[order(D)]
  qqplot(qexp(ppoints(length(D)), rate = .5), D, col = cols,
         main = 'Exponential QQ-plot',ylim = c(0,20),
         ylab = 'Mahalanobis squared distance',
         xlab = 'Exponential quantiles')
  abline(0, 1, col = 'gray')
  abline(h = cut, lty = 'dotted')
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
  mat$D = D
  
  #Subclonal CN segments given cutoff with 1-pexp(cut, .5) sign.level
  cols = rep('blue', nrow(mat))
  cols[mat$D > cut] = 'red'
  cols = col.transp(cols)
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
  
  #Choose how extreme CNs can be judged to be subclonal: CN subclonality defined here
  xlim = c(-1, quantile(mat$c2, probs = c(.95))*c(1.2))
  ylim = c(-1, quantile(mat$c1, probs = c(.95))*c(1.2))
  plot.cn(x = mat$c2, y = mat$c1, mycols = cols, xlim = xlim, ylim = ylim,
          xlab = 'Major CN', ylab = 'Minor CN')
  dev.print(png, file=file.path(obj$d, 
                                paste(obj$info$sample.name[obj$i],'_6.1_subclonal_CNs.png', 
                                      sep = '')),  width=900, height=900)
  obj$mat = mat
  obj$cut = cut
  obj
}




subclonal.segmuts = function(muts, mat, cn.cut = NULL, cn.cut0 = NULL){
  chroms = as.character(unique(mat$Chromosome))
  out = rep(FALSE, nrow(mat))
  for (ch in 1:length(chroms)){
    chrom = chroms[ch]
    matix = as.character(mat$Chromosome) == chrom
    IRsub = NULL
    endpoint = max(mat$End.bp[matix])
    binsize = 10000
    starts = seq(1, endpoint, binsize)
    ends = c(starts[-1]-1, endpoint)
    submuts = muts$Start_position[as.character(muts$Chromosome) == chrom & muts$subclonal]
    subtab = table(ceiling(submuts/binsize))
    y = rep(0, length(starts))
    y[as.numeric(names(subtab))] = subtab
    y[length(y)] = y[length(y)]*(ends[length(ends)]-starts[length(starts)]+1)/binsize
    x = 1:length(starts)
    smooth = ksmooth(x,y,  bandwidth=2000,
            kernel='normal',x.points = x)
#    plot(smooth$x, smooth$y)
    aboveind = which(smooth$y>0.001)
    if (length(aboveind)>0)  {
      IRsub = reduce(IRanges(start=aboveind, width=1))
      IRsub = IRanges(start = starts[start(IRsub)], end = ends[end(IRsub)])
    }
    IRseg = IRanges(start = mat$Start.bp[matix], end = mat$End.bp[matix])
    cv = coverage(c(IRseg, IRsub))
    out[matix] = viewMeans(Views(cv, IRseg))>1.5
    if (!is.null(cn.cut)){
      for (j in which(matix)){
        minor = mat$alpha[1]*mat$c1[j]/(mat$alpha[1]*(
          mat$c1[j]+mat$c2[j]) + (1-mat$alpha[1])*2)
        mutix = which(as.character(muts$Chromosome) == chrom &
                        muts$Start_position >= mat$Start.bp[j] &
                        muts$Start_position <= mat$End.bp[j])
        af = sum(muts$t_alt_count[mutix])/(sum(muts$t_alt_count[mutix]) +
                                             sum(muts$t_ref_count[mutix]))
        if (is.na(af)) out[j] = FALSE else {
          if (mat$c2[j]>cn.cut & mat$c1[j]<=cn.cut & af >= minor) out[j] = FALSE
        }
      }
    }
  }
  if (!is.null(cn.cut0)) out[mat$c2<cn.cut0 | mat$c1<cn.cut0] = FALSE
  if (!is.null(cn.cut)) out[mat$c2>cn.cut & mat$c1>cn.cut] = FALSE
  out
}  
ITH = function(obj, cn.cut, cn.cut0 = -1){
  mat = obj$mat
  mat$alpha = obj$alpha
  muts = obj$muts
  info = obj$info
  i = obj$i
  mat$subclonal = FALSE
  mat$subclonal[mat$D>obj$cut] = TRUE
  mat$subclonal[mat$c2>cn.cut | mat$c1>cn.cut] = FALSE
  mat$subclonal[mat$c2<cn.cut0 | mat$c1<cn.cut0] = FALSE
  
  #Subclonal segments wrt mutations
  mat$subclonal.muts = subclonal.segmuts(muts = muts,
                                         mat = mat, cn.cut = cn.cut, cn.cut0 = cn.cut0)
  
  #% ITH
  cols = rep('blue', nrow(mat))
  cols[mat$subclonal.muts] = 'darkorange1'
  cols[mat$subclonal] = 'red'
  cols = col.transp(cols)
  par(mfrow = c(1,1), mar = c(3, 3, 2, .1))
  xlim = c(-1, quantile(mat$c2, probs = c(.95))*c(1.2))
  ylim = c(-1, quantile(mat$c1, probs = c(.95))*c(1.2))
  plot.cn(x = mat$c2, y = mat$c1, mycols = cols, xlim = c(0,8), ylim = c(-.5,4),
          xlab = 'Major CN', ylab = 'Minor CN')
  ith = sum(mat$length[mat$subclonal | mat$subclonal.muts])/sum(mat$length)
  mtext(paste(round(ith*100), '% heterogeneity'), adj = 1)
  mtext(paste('Ploidy =', format(round(obj$ploidy[1],1), nsmall = 1)), adj = 0)
  legend('topleft', bty = 'n', pch = 16, col = c('blue','darkorange1', 'red'), 
         legend = c('Segment with integer CN',
                    'Segment with subclonal mutations',
                    'Segment with subclonal CN'))
  dev.print(png, file=file.path(obj$d, 
                                paste(info$sample.name[i],'_6.2_subclonal_CNs.png', 
                                      sep = '')),  width=400, height=400)
  
  obj$mat = mat
  obj$ith = ith
  obj$cn.cut = cn.cut
  obj$cn.cut0 = cn.cut0 
  obj
}

plotCNmutsGrid = function(obj){
  #Plot CN grid with mutations
  mat = obj$mat
  muts = obj$muts
  
  cols = rep('blue', nrow(mat))
  cols[mat$subclonal.muts] = 'darkorange1'
  cols[mat$subclonal] = 'red'
  cols = col.transp(cols, .7)
  xlim = c(-1, quantile(mat$c2, probs = c(.95))*c(1.2))
  ylim = c(-1, quantile(mat$c1, probs = c(.95))*c(1.2))
  plot.cn(x = mat$c2, y = mat$c1, mycols = cols, xlim = xlim, ylim = ylim,
          xlab = 'Major CN', ylab = 'Minor CN')
  mtext(paste(round(obj$ith*100), '% heterogeneity'), adj = 1)
  mtext(paste('Ploidy =', format(round(obj$ploidy,1), nsmall = 1)), adj = 0)
  points(mat$c2[muts$segment], mat$c1[muts$segment], pch = '.', cex = 4)
  points(mat$c2[(muts$segment)[muts$subclonal]], mat$c1[(muts$segment)[muts$subclonal]], 
         pch = 4, cex = 1.5)
  legend('topleft', bty = 'n', pch = c(16, 16, 16, 4, 20), col = c('blue','darkorange1', 'red', 'black', 'black'), 
         legend = c('Segment with integer CN', 'Segment with subclonal mutations', 'Segment with subclonal CN',
                    'Clonal mutation', 'Subclonal mutation'))
  dev.print(png, file=file.path(obj$d, 
                                paste(obj$info$sample.name[obj$i],'_7.2_subclonal_CNs_withMuts.png', 
                                      sep = '')),  width=400, height=400)  
}
plotCNmutsGenome = function(obj, restrictCN = FALSE){
  mat = obj$mat
  muts = obj$muts
  
  #Colour settings
  mutcols = c('green3','darkorange1')
  mycols = rep('blue', nrow(mat))
  mycols[mat$subclonal.muts] = 'darkorange1'
  mycols[mat$subclonal] = 'red'
  
  #Plot Copy numbers and multiplicities whole genome
  
  plot.CNmut(plotmat = mat, plotmuts = muts, mycols = mycols, mutcols = mutcols,
             restrict = restrictCN)
  mtext(paste(nrow(muts), 'mutations'))
  mtext(paste(round(obj$ith*100), '% heterogeneity'), adj = 1)
  mtext(paste('Ploidy =', format(round(obj$ploidy,1), nsmall = 1)), adj = 0)
  dev.print(png, file=file.path(obj$d, 
                                paste(obj$info$sample.name[obj$i],'_7.3_Genome_CNs_allMutations.png', 
                                      sep = '')),  width=940, height=600)
  muts = subset(muts, effect == 'nonsilent')
  plot.CNmut(plotmat = mat, plotmuts = muts, mycols = mycols,
             mutcols = mutcols, restrict = restrictCN)
  mtext(paste(nrow(muts), 'nonsilent mutations'))
  mtext(paste(round(obj$ith*100), '% heterogeneity'), adj = 1)
  mtext(paste('Ploidy =', format(round(obj$ploidy,1), nsmall = 1)), adj = 0)
  
  dev.print(png, file=file.path(obj$d, 
                                paste(obj$info$sample.name[obj$i],'_7.4_Genome_CNs_nonsilentMutations.png', 
                                      sep = '')),  width=940, height=600)  
}


plotCNmutsChrom = function(obj){
  plotmat = obj$mat
  muts = obj$muts
  
  #Colour settings
  mutcols = c('green3','darkorange1')
  mycols = rep('blue', nrow(plotmat))
  mycols[plotmat$subclonal.muts] = 'darkorange1'
  mycols[plotmat$subclonal] = 'red'
  
  #Plot Copy numbers and multiplicities chromosome by chromosome
  pdf(paste(obj$d, '/', obj$info$sample.name[obj$i], '_7.5_CNs_and_Mutations_byChromosome.pdf', 
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
}

###############
# Old
###############
subclonal.CNmuts = function(muts, mat, sign.level = .01, alpha = NULL){
  c1 = pmax(0, mat$c1)
  c2 = pmax(0, mat$c2)
  mat$out = FALSE
  for (k in 1:nrow(mat)){
    minor = alpha*c1[k]/(alpha*(c1[k]+c2[k]) + (1-alpha)*2)
    major = alpha*c2[k]/(alpha*(c1[k]+c2[k]) + (1-alpha)*2)
    both = alpha*(c1[k] + c2[k])/(alpha*(c1[k]+c2[k]) + (1-alpha)*2)
    ix = which(as.character(muts$Chromosome) == as.character(mat$chrom[k]) & 
                 muts$Start_position>=mat$Start.bp[k] &
                 muts$Start_position<=mat$End.bp[k])
    min = maj = both = FALSE
    if (length(ix)>0){
      min = binom.test(x = sum(muts$t_alt_count[ix]), n = sum(muts$t_alt_count[ix]) + sum(muts$t_ref_count[ix]),
                       p = minor)$p.value < sign.level
      maj = binom.test(x = sum(muts$t_alt_count[ix]), n = sum(muts$t_alt_count[ix]) + sum(muts$t_ref_count[ix]),
                       p = major)$p.value < sign.level
      both = binom.test(x = sum(muts$t_alt_count[ix]), n = sum(muts$t_alt_count[ix]) + sum(muts$t_ref_count[ix]),
                        p = both)$p.value < sign.level
    }
    mat$out[k] = min & maj & both
  }
  mat$out
}

##################################################################
##################################################################










##################################################################
##################################################################
##
## Plots of CNs and mutations by genome position
##
##################################################################
##################################################################

plot.CNmut = function(plotmat, plotmuts, ylim = c(-1, 8), c1 = NULL, c2 = NULL, mult = NULL,
                      subclonal = NULL,
                      main = NULL, xlab = 'Position along genome', ylab = 'Copies per cell',
                      mycols = NULL, mutcols = c('green3','red'),
                      restrict = FALSE){
  if (!is.null(c1)) plotmat$c1 = c1
  if (!is.null(c2)) plotmat$c2 = c2
  if (!is.null(mult)) plotmuts$mult = mult
  if (!is.null(subclonal)) plotmuts$subclonal = subclonal
  plotmuts$cols = mutcols[1]
  plotmuts$cols[plotmuts$subclonal] = mutcols[2]
  
  #Restrict CNs to plotting area
  plotmat$c1 = pmax(plotmat$c1, 0)
  plotmat$c2 = pmax(plotmat$c2, 0)
  plotmat$tot = plotmat$c1 + plotmat$c2
  plotmat$c1 = pmin(plotmat$c1, ylim[2])
  plotmat$c2 = pmin(plotmat$c2, ylim[2])
  plotmat$tot = pmin(plotmat$tot, ylim[2])
  
  #Find genome coordinates for x-axes
  chroms = as.character(sort(unique(plotmat$Chromosome)))
  par(mar = c(2, 4, 2, .1))
  tab = tapply(c(plotmat$End.bp,plotmuts$End_position), 
               c(plotmat$Chromosome, plotmuts$Chromosome), max, na.omit = T)
  tab = tab[order(as.numeric(names(tab)))]
  stopifnot(all(names(tab) == chroms)) #Check that chromosome names agree, should be TRUE
  chrlens = c(tab)
  chrends = diffinv(chrlens)[-1]
  chrstarts = c(0,chrends[-length(chrends)])+1
  chrmids = (chrstarts + chrends)/2
  xlim = c(0, chrends[length(chrends)]+1)
  
  #Plot frame and gridlines
  plot(1, 1, type = 'n', ylim = ylim, 
       xlim = xlim, ylab = ylab, 
       xlab = '', xaxt = 'n')
  title(xlab = xlab, line = 1)
  for (k in 0:round(ylim[2])){
    abline(h = k, lty = 'dotted', col = 'grey20')
  }
  vs = unique(c(chrstarts-.5,chrends+.5))
  for (k in 1:length(vs)){
    lines(rep(vs[k],2), c(-.5, ylim[2]), lty = 'dotted', col = 'grey20')
  }
  text(chrmids, rep(-1, length(chrmids)), labels = chroms, adj = c(.5,0),
       cex = .8)
  title(main = main)
  
  if (is.list(mycols)){
    #Colour scale
    colscale = mycols$cols
    tmp = c(xlim[2]/2, xlim[2])
    tmp = seq(tmp[1], tmp[2], length=length(colscale))
    points(tmp, rep(-1.2, length(tmp)), col = colscale, pch = 16)
    
    mycols = mycols$mycols
  }
  
  if (is.null(mycols)) mycols = rep('blue', nrow(plotmat))
  for (j in 1:length(chroms)){
    chrom = chroms[j]
    
    ix.plotmat = which(as.character(plotmat$Chromosome) == chrom)
    mystart = chrstarts[j]-1
    for (k in ix.plotmat){
      same = abs(plotmat$tot[k] - plotmat$c2[k]<.5) | abs(plotmat$c2[k] - plotmat$c1[k])<.5
      if (!restrict | same){
        lines(mystart + c(plotmat$Start.bp[k], plotmat$End.bp[k]), rep(plotmat$tot[k], 2), 
              lty  = 'solid', lwd = 7, col = 'grey', lend = 'butt')
        lines(mystart + c(plotmat$Start.bp[k], plotmat$End.bp[k]), rep(plotmat$c1[k], 2), 
              lty  = 'solid', lwd = 7, col = mycols[k], lend = 'butt')
        lines(mystart + c(plotmat$Start.bp[k], plotmat$End.bp[k]), rep(plotmat$c2[k], 2), 
              lty  = 'solid', lwd = 7, col = mycols[k], lend = 'butt')
        ix = which(as.character(plotmuts$Chromosome) == chrom & 
                     plotmuts$Start_position>=plotmat$Start.bp[k] &
                     plotmuts$Start_position<=plotmat$End.bp[k])
        if (length(ix)>0){
          points(mystart+plotmuts$Start_position[ix], pmin(plotmuts$mult[ix], ylim[2]), 
                 pch = 2,
                 col = plotmuts$cols[ix], cex = .4)
        }
      }
    }    
  }
}

##################################################################
##################################################################









##################################################################
##################################################################
##
## HAPSEG output format
##
##################################################################
##################################################################


hapseg.format = function(seg.dat, mat, xname = 'c1', yname = 'c2'){
  
  #Indata
  mat$a1p = mat[,yname]
  mat$a2p = mat[,xname]
  totalsum = sum(c(mat$a1p, mat$a2p)* rep(mat$W, 2)/2)
  mat$a1p = mat$a1p/totalsum
  mat$a2p = mat$a2p/totalsum
  
  #Reset attenuation constant
  seg.dat$error.model$AT = 0
  
  #New version or hscr values in as.seg.dat
  as.seg.dat = seg.dat$as.seg.dat
  as.seg.dat$which = 1
  nr = as.integer(nrow(as.seg.dat))
  as.seg.dat$which[(nr/2+1):nr] = 2
  upper = as.seg.dat[as.seg.dat$which == 1,]
  tmp = mat[,c('Chromosome','Start.bp','End.bp', 'a1p')]
  upper = merge(upper, tmp, all = T)
  upper$copy_num = upper$a1p
  upper = upper[order(upper$seg.ix),1:9]
  lower = as.seg.dat[as.seg.dat$which == 2,]
  tmp = mat[,c('Chromosome','Start.bp','End.bp', 'a2p')]
  lower = merge(lower, tmp, all = T)
  lower$copy_num = lower$a2p
  lower = lower[order(lower$seg.ix),1:9]
  as.seg.dat = rbind(upper, lower)
  seg.dat$as.seg.dat = as.seg.dat
  
  
  #New version of allele.segs
  allele.segs = seg.dat$allele.segs
  allele.segs = as.data.frame(allele.segs)
  tmp = mat[,c('Chromosome','Start.bp','End.bp', 'a1p', 'a2p')]
  allele.segs = merge(allele.segs, tmp, sort = FALSE)
  allele.segs$A1.Seg.CN = allele.segs$a1p
  allele.segs$A2.Seg.CN = allele.segs$a2p
  allele.segs = as.matrix(allele.segs[,1:9])
  seg.dat$allele.segs = allele.segs
  
  #New version of seg.info
  seg.info = seg.dat$seg.info
  seg.info = as.data.frame(seg.info)
  tmp = mat[,c('Chromosome','Start.bp','End.bp', 'a1p', 'a2p')]
  seg.info = merge(seg.info, tmp, sort = FALSE)
  seg.info$cn = (seg.info$a1p + seg.info$a2p)/2
  seg.info$copy_num = seg.info$cn
  seg.info = as.matrix(seg.info[,1:8])
  seg.dat$seg.info = seg.info
  
  #Now switch over so that A2 is always the largest (keep index of switched ones)
  ix = which(seg.dat$allele.segs[,'A1.Seg.CN'] > seg.dat$allele.segs[,'A2.Seg.CN'])
  a1 = seg.dat$allele.segs[,'A1.Seg.CN']
  a2 = seg.dat$allele.segs[,'A2.Seg.CN']
  seg.dat$allele.segs[,'A1.Seg.CN'] = pmin(a1, a2)
  seg.dat$allele.segs[,'A2.Seg.CN'] = pmax(a1, a2)
  seg.dat$switched.allele.segs = ix
  
  nr = nrow(seg.dat$as.seg.dat)
  a1 = seg.dat$as.seg.dat[1:(nr/2),'copy_num']
  a2 = seg.dat$as.seg.dat[(nr/2+1):(nr),'copy_num']
  seg.dat$as.seg.dat[1:(nr/2),'copy_num'] = pmin(a1, a2)
  seg.dat$as.seg.dat[(nr/2+1):(nr),'copy_num'] = pmax(a1, a2)
  seg.dat
}
##################################################################
##################################################################










##################################################################
##################################################################
##
## Old rotation
##
##################################################################
##################################################################
#A1 vs A2
plot1v2 = function(a1, a2, mat, cl=NULL, clust=NULL, w=NULL, xlab = 'Homologue 2 copy value',
                   ylab = 'Homologue 1 copy value', ...){
  cols = c(2:7, 'gold' , colors()[c(96)])
  for (c in 7:length(cols)) cols[c] = col.transp(cols[c])
  x = a2
  y = a1
  ix = !(y %in% c(Inf, -Inf, NA, NaN) | x %in% c(Inf, -Inf, NA, NaN))
  size = mat$length/sum(mat$length)*100
  Wmod = lm(size~length, data = mat)
  ranx = quantile(x[ix], prob = c(.05, .95))
  ranx = ranx + c(-1, 1)*diff(ranx)/10
  rany = quantile(y[ix], prob = c(.05, .95))
  rany = rany + c(-1, 1)*diff(rany)/10
  plot(x, y, xlim = ranx, ylim = rany, cex = size, xlab = xlab,
       ylab = ylab, ...)
  if (!is.null(clust)){
    ix = which(clust == 'cluster het')
    clust.W = predict(Wmod, data.frame(length=w[ix]))
    points(x[ix], y[ix], col =cols[8], pch = 16, cex = clust.W)
    try({ix = which(clust == 'cluster hom')
         clust.W = predict(Wmod, data.frame(length=w[ix]))
         points(x[ix], y[ix], col =cols[7], pch = 16, cex = clust.W)
    }, silent = T)
    for (c in seq_along(unique(cl$cluster))){
      myclust = unique(cl$cluster)[c]
      if (!(myclust %in% c('het','hom'))){
        ix = which(clust == paste('cluster', myclust))
        clust.W = predict(Wmod, data.frame(length=w[ix]))
        points(x[ix], y[ix], col =cols[c], pch = 16, cex = clust.W)
      }
    }
    grid()
    clusters = unique(cl$cluster)
    nc = length(clusters)
    clusters = clusters[c(nc-1, nc, 1:(nc-2))]
    lcols = cols[c(7, 8, 1:(nc-2))]
    legend('bottomright', col = lcols, pch = 16, legend = c('CN2 Hom','CN2 Het',
                                                            paste('Cluster', 1:(nc-2))), bty = 'n')
  }
}

histCol = function(toplot = c('A1','A2'), mat = NULL, cl, clust, myab = NULL, xlim = c(0,3), xlab = 'Haplytype copy values'){
  cols = c(2:7, 'gold' , colors()[c(96)])
  for (c in 7:length(cols)) cols[c] = col.transp(cols[c])
  nclust = length(unique(cl$cluster))
  colvec = c(as.numeric(cols[c(1:(nclust-2))]), 6, 7)
  if (is.null(myab)){
    myab = mat
    size = mat$length/sum(mat$length)*100
    Wmod = lm(size~length, data = mat)
    myab$W = size
    myab$color.by = 8
    for (c in seq_along(unique(cl$cluster))){
      myclust = unique(cl$cluster)[c]
      mycl = rep(0, nrow(myab))
      tmpclust = c(clust, rep(0, nrow(myab)-length(clust)))
      mycl[tmpclust == paste('cluster', myclust)] = w[tmpclust == paste('cluster', myclust)]
      if (myclust == 'het') myclust = nclust - 1
      if (myclust == 'hom') myclust = nclust
      
      clust.W = predict(Wmod, data.frame(length=mycl))
      clust.W[mycl == 0] = 0
      myab$W = pmax(0, myab$W - clust.W)
      ix = mycl>0
      toadd = myab[ix,]
      toadd$W = clust.W[ix]
      if (length(ix[ix])>0){
        toadd$color.by= colvec[c]
        if (myclust == nclust-1) toadd$color.by = colvec[nclust-1]
        if (myclust == nclust) toadd$color.by = colvec[nclust]
      }
      myab = rbind(myab, toadd)
    }
  }
  x = unlist(myab[,toplot])
  l = rep(myab$W, length = length(x))
  colseg = rep(myab$color.by, length = length(x))
  
  ABSOLUTE:::PlotSeglenHist(x, l, color.by = colseg, use.pal = 2:8, 
                            ylab = '', 
                            data.whiskers = F,
                            xlab = xlab, xlim = xlim)
  grid()
  clusters = unique(cl$cluster)
  nc = length(clusters)
  clusters = clusters[c(nc-1, nc, 1:(nc-2))]
  lcols = colvec[c(nc, nc-1, 1:(nc-2))]
  legend('topright', col = lcols, pch = 16, legend = c('CN2 Hom','CN2 Het',
                                                       paste('Cluster', 1:(nc-2))), bty = 'n')
  myab
}

minx = function(pars, X, Y, xclust, xw, myclust){
  require(limma)
  a = pars[1]
  b = pars[2]
  c = 0
  d = 1
  xp = a*X + b*Y
  yp = c*X + d*Y
  groups = as.character(myclust)
  gdat = data.frame(group = unique(groups), y. = NA, cn = NA, weight = NA)
  for (i in seq_along(gdat$group)){
    ix = which(groups == gdat$group[i])
    gdat$y.[i] =weighted.mean(x = xp[ix], w = xw[ix])
    gdat$cn[i] = unique(xclust[ix])
    gdat$weight[i] = sum(xw[ix])
  }
  cndat = data.frame(cn = unique(xclust), y.. = NA)
  for (i in seq_along(cndat$cn)){
    ix = which(xclust == cndat$cn[i])
    cndat$y..[i] = weighted.mean(x = xp[ix], w = xw[ix], na.rm = T)
  }
  gdat = merge(gdat, cndat, all = T)
  gdat$ss = (gdat$y. - gdat$y..)^2
  sum(gdat$ss)
}
eqfunx = function(pars, X, Y, xclust, xw, myclust){
  A = pars[1]
  B = pars[2]
  c = 0
  d = 1
  xp = A*X + B*Y
  mean(xp)
}
miny = function(pars, X, Y, yclust, yw, myclust){
  require(limma)
  a = 1
  b = 0
  c = pars[1]
  d = pars[2]
  xp = a*X + b*Y
  yp = c*X + d*Y
  groups = as.character(myclust)
  gdat = data.frame(group = unique(groups), y. = NA, cn = NA, weight = NA)
  for (i in seq_along(gdat$group)){
    ix = which(groups == gdat$group[i])
    gdat$y.[i] =weighted.mean(x = yp[ix], w = yw[ix])
    gdat$cn[i] = unique(yclust[ix])
    gdat$weight[i] = sum(yw[ix])
  }
  cndat = data.frame(cn = unique(yclust), y.. = NA)
  for (i in seq_along(cndat$cn)){
    ix = which(yclust == cndat$cn[i])
    cndat$y..[i] = weighted.mean(x = yp[ix], w = yw[ix], na.rm = T)
  }
  gdat = merge(gdat, cndat, all = T)
  gdat$ss = (gdat$y. - gdat$y..)^2
  sum(gdat$ss)
}
eqfuny = function(pars, X, Y, yclust, yw, myclust){
  b = 1
  a = 0
  c = pars[1]
  d = pars[2]
  yp = c*X + d*Y
  mean(yp)
}
minxy = function(pars, X, Y, xclust, yclust, xw, yw, mat, myclust){
  require(limma)
  A = 0
  B = pars[1]
  C = 0
  D = pars[2]
  xp = A + B*X
  yp = C + D*Y
  groups = c(as.character(interaction('x', xclust)), as.character(interaction('y', yclust)))
  gdat = data.frame(group = unique(groups), y. = NA, cn = NA, weight = NA)
  for (i in seq_along(gdat$group)){
    ix = which(groups == gdat$group[i])
    gdat$y.[i] =weighted.mean(x = c(xp, yp)[ix], w = c(xw, yw)[ix])
    if (!is.na(gdat$group[i])) gdat$cn[i] = unique(c(xclust, yclust)[ix])
    gdat$weight[i] = sum(c(xw, yw)[ix])
  }
  cndat = data.frame(cn = unique(c(xclust, yclust)), y.. = NA)
  for (i in seq_along(cndat$cn)){
    ix = which(c(xclust, yclust) == cndat$cn[i])
    cndat$y..[i] = weighted.mean(x = c(xp, yp)[ix], w = c(xw, yw)[ix], na.rm = T)
  }
  gdat = merge(gdat, cndat, all = T)
  gdat = subset(gdat, !is.na(cn))
  gdat$ss = (gdat$y. - gdat$y..)^2
  sum(gdat$ss*gdat$weight)
}
eqfuncxy = function(pars, X, Y, xclust, yclust, xw, yw, mat, myclust){
  A = 0
  B = pars[1]
  C = 0
  D = pars[2]
  a2 = A + B*mat$A2new
  a1 = C + D*mat$A1new
  sum(c(a1, a2)* rep(mat$W, 2)/2)
}

getKnownClustersAndOverlaps = function(mat, d){
  subsegs = read.delim(file.path(d,'Subsegs.txt'))
  homsegs = read.delim(file.path(d, 'Homsegs.txt'))
  homsegs$cluster = 'hom'
  hetsegs = read.delim(file.path(d, 'Hetsegs.txt'))
  hetsegs$cluster = 'het'
  cl = rbind(subsegs, homsegs, hetsegs)
  cl = subset(cl, !(Chromosome %in% c('chrX', 'chrY')))
  cl$Chromosome = as.character(cl$Chromosome)
  nc = nchar(cl$Chromosome, type = 'c')
  cl$Chromosome = as.numeric(substr(cl$Chromosome, 4, nc))
  cl$Chromosome = factor(cl$Chromosome, levels = as.character(1:22))
  
  ###Find cluster overlapped segments 
  clust = rep('', nrow(mat))
  w = rep(0, nrow(mat))
  for (c in seq_along(unique(cl$cluster))){
    cc = subset(cl, cluster == unique(cl$cluster)[c])
    for (ch in seq_along(levels(cc$Chromosome))){
      mat.ix = which(mat$Chromosome == levels(cc$Chromosome)[ch])
      cc.ix = which(cc$Chromosome == levels(cc$Chromosome)[ch])
      IRmat = IRanges(start = mat$Start.bp[mat.ix], end = mat$End.bp[mat.ix])
      IRcc = IRanges(start = cc$Start[cc.ix], end = cc$End[cc.ix])
      if (length(IRcc)>0 & length(IRmat)>0){
        cov = coverage(c(IRmat, IRcc)) - 1
        ls = viewSums(Views(cov, IRmat))
        w[mat.ix[ls>0]] = ls[ls>0]
        clust[mat.ix[ls>0]]= paste('cluster', unique(cl$cluster)[c])
      } 
    }  
  }
  list(clust = clust, w = w, cl = cl)
}
##################################################################
##################################################################







##################################################################
##################################################################
##
## ABSOLUTE and HAPSEG functions
##
##################################################################
##################################################################


myRunAbsolute = function (seg.dat.fn, sigma.p, max.sigma.h, min.ploidy, max.ploidy, 
                        primary.disease, platform, sample.name, results.dir, max.as.seg.count, 
                        max.non.clonal, max.neg.genome, copy_num_type, maf.fn = NULL, 
                        min.mut.af = NULL, output.fn.base = NULL, verbose = FALSE) 
{
  kQ = 15
  kTauDom = c(min.ploidy, max.ploidy)
  kSigmaHDom = c(0, max.sigma.h)
  kPiSomThetaQ = c(100, 50, rep(2, (kQ - 2)))
  mut_class_w = list(SM = 0.5, GL = 0, SC = 0.5, OL = 0.001, 
                     Pi_SM = 15, Pi_SC = 15)
  platform = match.arg(platform, c("SNP_6.0", "Illumina_WES", 
                                   "SNP_250K_STY"))
  if (platform %in% c("SNP_6.0", "SNP_250K_STY")) {
    filter_segs = FALSE
  } else if (platform == "Illumina_WES") {
    filter_segs = TRUE
  } else {
    stop("Unsupported platform: ", platform)
  }
  min_probes = 10
  max_sd = 100
  copy_num_type = match.arg(copy_num_type, c("allelic", "total"))
  if (copy_num_type == "total") {
    pi_theta_qz = list(W = c(1, 25, 100, 25, 10, 5, rep(1, 
                                                        kQ - 6), 10), PC = 0.05)
    ABSOLUTE:::set_total_funcs()
  } else if (copy_num_type == "allelic") {
    pi_theta_qz = list(W = c(25, 100, 25, 10, 5, rep(1, kQ - 
                                                       5), 10), PC = 0.05)
    ABSOLUTE:::set_allelic_funcs()
  }  else {
    stop("Unsupported copy number type: ", copy_num_type)
  }
  SubclonalPost <<- ABSOLUTE:::UnifSubclonalPost
  sigma.h = sigma.p
  tmp.dir = file.path(results.dir, "tmp")
  dir.create(tmp.dir, recursive = TRUE, showWarnings = FALSE)
  dir.create(results.dir, recursive = TRUE, showWarnings = FALSE)
  seg.dat = MakeSegObj(seg.dat.fn, min_probes = min_probes, 
                       max_sd = max_sd, filter_segs = filter_segs, verbose = verbose)
  seg.dat[["primary.disease"]] = primary.disease
  seg.dat[["group"]] = ABSOLUTE:::DetermineGroup(primary.disease)
  seg.dat[["platform"]] = platform
  seg.dat[["sample.name"]] = as.character(sample.name)
  if (is.null(seg.dat$array.name)) {
    seg.dat$array.name = seg.dat$sample.name
  }
  seg.dat[["maf.fn"]] = maf.fn
  seg.dat[["obs.scna"]] = ExtractSampleObs(seg.dat)
  seg = seg.dat[["obs.scna"]][["segtab"]]
  e.cr = sum(seg[, "W"] * seg[, "copy_num"])
  if (verbose) {
    print(paste("Expected copy-ratio = ", round(e.cr, 5), 
                sep = ""))
  }
  mode.res = list(mode.flag = NA)
  if (nrow(seg) > max.as.seg.count) {
    mode.res[["mode.flag"]] = "OVERSEG"
  }
  if ((e.cr < 0.75) || (e.cr > 1.25)) {
    mode.res[["mode.flag"]] = "E_CR_SCALE"
  }
  if (is.na(mode.res[["mode.flag"]])) {
    maf = NULL
    if ((!is.null(maf.fn)) && (file.exists(maf.fn))) {
      maf = read.delim(maf.fn, row.names = NULL, stringsAsFactors = FALSE, 
                       check.names = FALSE, na.strings = c("NA", "---"), 
                       blank.lines.skip = TRUE, comment.char = "#")
    } else {
      if (verbose) {
        print(paste("MAF file: ", maf.fn, " not found.", 
                    sep = ""))
      }
    }
    data(ChrArmsDat, package = "ABSOLUTE")
    if ((!is.null(maf)) && (nrow(maf) > 0)) {
      if (is.null(min.mut.af)) {
        stop("min.mut.af must be defined")
      }
      mut.cn.dat = ABSOLUTE:::CreateMutCnDat(maf, seg.dat, min.mut.af, 
                                             verbose = verbose)
      mut = list(mut_CN_dat = mut.cn.dat, Pi_som_theta_Q = kPiSomThetaQ, 
                 mut_class_W = mut_class_w)
    } else {
      mut = NA
    }
    mode.res = ABSOLUTE:::FitSample(seg.dat, mut, kQ, pi_theta_qz, sigma.h, 
                                    kTauDom, kSigmaHDom, chr.arms.dat, verbose = verbose)
    mode.res = myapply_subclonal_scna_model(segobj = seg.dat, mode_res = mode.res, 
                                            verbose = verbose)
    if (inherits(mode.res, "try-error")) {
      mode.res = list(mode.flag = "FAIL")
    }
  }
  if (is.na(mode.res[["mode.flag"]])) {
    bad.ix = ABSOLUTE:::GenomeHetFilter(seg.dat[["obs.scna"]], mode.res, 
                                        max.non.clonal, max.neg.genome, kQ, verbose = verbose)
    if (sum(bad.ix) == nrow(mode.res[["mode.tab"]])) {
      mode.res = list(mode.flag = "ALPHA_TAU_DOM")
    } else {
      mode.res = ABSOLUTE:::ReorderModeRes(mode.res, !bad.ix)
    }
  }
  if (is.na(mode.res[["mode.flag"]])) {
    data(ChrArmPriorDb, package = "ABSOLUTE")
    model.id = ifelse(seg.dat[["group"]] %in% names(train.obj), 
                      seg.dat[["group"]], "Primary")
    mode.res = ABSOLUTE:::ApplyChrArmPrior(mode.res, model.id, train.obj, 
                                           verbose = verbose)
    if ((!is.null(maf)) && (nrow(maf) > 0)) {
      if (is.null(min.mut.af)) {
        stop("maf was not NULL, but min.mut.af was")
      }
      seg.dat[["mut.cn.dat"]] = mut.cn.dat
      mode.res = ABSOLUTE:::ApplySomaticMutsModel(mode.res, seg.dat$obs_SCNA, 
                                                  mut.cn.dat, kPiSomThetaQ, mut_class_w, kQ, verbose = verbose)
    }
    mode.res = ABSOLUTE:::WeighSampleModes(mode.res)
    mode.res[["call.status"]] = ABSOLUTE:::GetCallStatus(mode.res, seg.dat[["obs.scna"]][["W"]])
  }
  seg.dat[["mode.res"]] = mode.res
  if (is.null(output.fn.base)) {
    output.fn.base = ifelse(is.null(seg.dat$array.name), 
                            sample.name, seg.dat$array.name)
  }
  file.base = paste(output.fn.base, ".ABSOLUTE", sep = "")
  if (is.na(mode.res[["mode.flag"]])) {
    sample.pdf.fn = file.path(results.dir, paste(file.base, 
                                                 "plot.pdf", sep = "_"))
    ABSOLUTE:::AbsoluteResultPlot(sample.pdf.fn, seg.dat)
  } else {
    if (verbose) {
      print("Mode flag is NA, not generating plots. Sample has failed ABSOLUTE")
    }
  }
  seg.dat$version = 1.1
  save(seg.dat, file = file.path(results.dir, paste(file.base, 
                                                    "RData", sep = ".")))
  return(TRUE)
}


myapply_subclonal_scna_model = function (segobj, mode_res, verbose = FALSE) 
{
  Q = dim(mode_res$theta.q.tab)[2]
  M = nrow(mode_res$mode.tab)
  n_seg = nrow(segobj$obs.scna$segtab)
  if (verbose) {
    print(paste("Evaluating subclonal SCNAs in ", M, " purity/ploidy modes: ", 
                sep = ""))
  }
  for (i in seq_len(M)) {
    res = myget_scna_cancer_cell_fractions(segobj, mode_res, 
                                           mode_ix = i)
    if (i == 1) {
      n_col = ncol(res$subclonal_SCNA_tab)
      subclonal_scna_tab = array(NA, dim = c(M, n_seg, 
                                             n_col))
      dimnames(subclonal_scna_tab)[[3]] = colnames(res$subclonal_SCNA_tab)
      n_col = ncol(res$CCF_dens)
      ccf_dens = array(NA, dim = c(M, n_seg, n_col))
      dimnames(ccf_dens)[[3]] = colnames(res$CCF_dens)
    }
    subclonal_scna_tab[i, , ] = res$subclonal_SCNA_tab
    ccf_dens[i, , ] = res$CCF_dens
  }
  mode_res$subclonal_SCNA_res = list(subclonal_SCNA_tab = subclonal_scna_tab, 
                                     CCF_dens = ccf_dens)
  return(mode_res)
}



myget_scna_cancer_cell_fractions = function (segobj, mode_res, mode_ix, pr_subclonal_threshold = 0.2) 
{
  mode_tab = mode_res$mode.tab
  pr_subclonal = mode_res$seg.z.tab[mode_ix, ]
  seg_q_tab = mode_res$seg.q.tab[mode_ix, , ]
  seg_qz_tab = mode_res$seg.qz.tab[mode_ix, , ]
  obs = ExtractSampleObs(segobj)
  obs$error.model$fit.at = mode_tab[mode_ix, "AT"]
  n_seg = nrow(seg_q_tab)
  sigma_h = mode_tab[mode_ix, "sigma.h.hat"]
  tau = mode_tab[mode_ix, "tau"]
  alpha = mode_tab[mode_ix, "alpha"]
  delta = alpha/(2 * (1 - alpha) + alpha * tau)
  b = 2 * (1 - alpha)/(2 * (1 - alpha) + alpha * tau)
  comb = InvTxData(obs$error.model, GetCopyRatioComb(15, delta, 
                                                     b, obs$error.model))
  seg_qz_hat = apply(seg_qz_tab, 1, which.max) - 1
  modal_cn = which.max(colSums(seg_q_tab)) - 1
  scna_ix = seg_qz_hat != modal_cn
  subclonal_ix = pr_subclonal >= pr_subclonal_threshold
  res = myget_subclonal_states(obs, subclonal_ix, modal_cn, comb)
  qc = res$qc
  qs = res$qs
  ccf_hat = rep(NA, n_seg)
  ccf_hat[!subclonal_ix] = 1
  ccf_ci95 = matrix(NA, nrow = n_seg, ncol = 2)
  ccf_ci95[!subclonal_ix, ] = 1
  ccf_grid = seq(0.01, 1, length = 100)
  ccf_dens = matrix(NA, nrow = n_seg, ncol = length(ccf_grid))
  colnames(ccf_dens) = ccf_grid
  cr_dens = matrix(NA, nrow = n_seg, ncol = length(ccf_grid))
  cr_grid = matrix(NA, nrow = n_seg, ncol = length(ccf_grid))
  for (i in 1:nrow(ccf_dens)) {
    if (!subclonal_ix[i]) {
      next
    }
    if (is.na(qs[i])) {
      next
    }
    d = qs[i] - qc[i]
    dd = (comb[qs[i] + 1] - comb[qc[i] + 1])
    cr_grid[i, ] = (dd * ccf_grid) + comb[qc[i] + 1]
    cr_dens[i, ] = GetScnaStderrGridDensity(obs, grid = cr_grid[i, 
                                                                ], sigma_h, i)[i,]
    ccf_dens[i, ] = cr_dens[i, ]/sum(cr_dens[i, ])
    ccf_hat[i] = ccf_grid[which.max(ccf_dens[i, ])]
    ecdf = cumsum(ccf_dens[i, ])
    ccf_ci95[i, ] = approx(x = ecdf, y = ccf_grid, xout = c(0.025, 
                                                            0.975))$y
  }
  nix1 = is.na(ccf_ci95[, 1])
  ccf_ci95[nix1, 1] = min(ccf_grid)
  nix2 = is.na(ccf_ci95[, 2])
  ccf_ci95[nix2, 2] = max(ccf_grid)
  ix = ccf_ci95[, 2] > ccf_grid[length(ccf_grid) - 1]
  ccf_ci95[ix, 2] = 1
  colnames(ccf_ci95) = c("CI95_low", "CI95_high")
  subclonal_scna_tab = cbind(CCF_hat = ccf_hat, subclonal_ix = subclonal_ix, 
                             SCNA_ix = scna_ix, ccf_ci95, Pr_subclonal = pr_subclonal, 
                             qs = qs, qc = qc)
  res = list(subclonal_SCNA_tab = subclonal_scna_tab, CCF_dens = ccf_dens)
  return(res)
}


myget_subclonal_states = function (obs, subclonal_ix, modal_cn, comb) 
{
  cr_set = InvTxData(obs$error.model, obs$d.tx)
  qc = rep(modal_cn, length(cr_set))
  qs = rep(NA, length(cr_set))
  del_ix = subclonal_ix & (cr_set < comb[modal_cn + 1])
  gain_ix = subclonal_ix & (cr_set > comb[modal_cn + 1])
  Q = length(comb)
  del_vals = (c(0:(modal_cn - 1)))
  for (i in seq_along(del_vals)) {
    qs[del_ix & (cr_set > comb[del_vals[i] + 1])] = del_vals[i]
  }
  qs[del_ix & cr_set < comb[1]] = NA
  qc[del_ix] = qs[del_ix] + 1
  gain_vals = rev(c((modal_cn + 1):Q))
  for (i in seq_along(gain_vals)) {
    qs[gain_ix & (cr_set < comb[gain_vals[i] + 1])] = gain_vals[i]
  }
  qc[gain_ix] = qs[gain_ix] - 1
  res = list(qs = qs, qc = qc)
  return(res)
}

myMakeSegObj = function(seg.dat.fn, min_probes = NA, 
                        max_sd = NA, filter_segs = FALSE, verbose = FALSE){
  load(seg.dat.fn)
  seg.dat
}
myRunAbsolute2 = function (seg.dat.fn, sigma.p, max.sigma.h, min.ploidy, max.ploidy, 
                          primary.disease, platform, sample.name, results.dir, max.as.seg.count, 
                          max.non.clonal, max.neg.genome, copy_num_type, maf.fn = NULL, 
                          min.mut.af = NULL, output.fn.base = NULL, verbose = FALSE, b.res = .125, d.res = .125) 
{
  kQ = 15
  kTauDom = c(min.ploidy, max.ploidy)
  kSigmaHDom = c(0, max.sigma.h)
  kPiSomThetaQ = c(100, 50, rep(2, (kQ - 2)))
  mut_class_w = list(SM = 0.5, GL = 0, SC = 0.5, OL = 0.001, 
                     Pi_SM = 15, Pi_SC = 15)
  platform = match.arg(platform, c("SNP_6.0", "Illumina_WES", 
                                   "SNP_250K_STY"))
  if (platform %in% c("SNP_6.0", "SNP_250K_STY")) {
    filter_segs = FALSE
  } else if (platform == "Illumina_WES") {
    filter_segs = TRUE
  } else {
    stop("Unsupported platform: ", platform)
  }
  min_probes = 10
  max_sd = 100
  copy_num_type = match.arg(copy_num_type, c("allelic", "total"))
  if (copy_num_type == "total") {
    pi_theta_qz = list(W = c(1, 25, 100, 25, 10, 5, rep(1, 
                                                        kQ - 6), 10), PC = 0.05)
    ABSOLUTE:::set_total_funcs()
  } else if (copy_num_type == "allelic") {
    pi_theta_qz = list(W = c(25, 100, 25, 10, 5, rep(1, kQ - 
                                                       5), 10), PC = 0.05)
    ABSOLUTE:::set_allelic_funcs()
  }  else {
    stop("Unsupported copy number type: ", copy_num_type)
  }
  SubclonalPost <<- ABSOLUTE:::UnifSubclonalPost
  sigma.h = sigma.p
  tmp.dir = file.path(results.dir, "tmp")
  dir.create(tmp.dir, recursive = TRUE, showWarnings = FALSE)
  dir.create(results.dir, recursive = TRUE, showWarnings = FALSE)
  seg.dat = myMakeSegObj(seg.dat.fn, min_probes = min_probes, 
                       max_sd = max_sd, filter_segs = filter_segs, verbose = verbose)
  seg.dat[["primary.disease"]] = primary.disease
  seg.dat[["group"]] = ABSOLUTE:::DetermineGroup(primary.disease)
  seg.dat[["platform"]] = platform
  seg.dat[["sample.name"]] = as.character(sample.name)
  if (is.null(seg.dat$array.name)) {
    seg.dat$array.name = seg.dat$sample.name
  }
  seg.dat[["maf.fn"]] = maf.fn
  seg.dat[["obs.scna"]] = ExtractSampleObs(seg.dat)
  seg = seg.dat[["obs.scna"]][["segtab"]]
  e.cr = sum(seg[, "W"] * seg[, "copy_num"])
  if (verbose) {
    print(paste("Expected copy-ratio = ", round(e.cr, 5), 
                sep = ""))
  }
  mode.res = list(mode.flag = NA)
  if (nrow(seg) > max.as.seg.count) {
    mode.res[["mode.flag"]] = "OVERSEG"
  }
  if ((e.cr < 0.75) || (e.cr > 1.25)) {
    mode.res[["mode.flag"]] = "E_CR_SCALE"
  }
  if (is.na(mode.res[["mode.flag"]])) {
    maf = NULL
    if ((!is.null(maf.fn)) && (file.exists(maf.fn))) {
      maf = read.delim(maf.fn, row.names = NULL, stringsAsFactors = FALSE, 
                       check.names = FALSE, na.strings = c("NA", "---"), 
                       blank.lines.skip = TRUE, comment.char = "#")
    } else {
      if (verbose) {
        print(paste("MAF file: ", maf.fn, " not found.", 
                    sep = ""))
      }
    }
    data(ChrArmsDat, package = "ABSOLUTE")
    if ((!is.null(maf)) && (nrow(maf) > 0)) {
      if (is.null(min.mut.af)) {
        stop("min.mut.af must be defined")
      }
      mut.cn.dat = ABSOLUTE:::CreateMutCnDat(maf, seg.dat, min.mut.af, 
                                             verbose = verbose)
      mut = list(mut_CN_dat = mut.cn.dat, Pi_som_theta_Q = kPiSomThetaQ, 
                 mut_class_W = mut_class_w)
    } else {
      mut = NA
    }
    mode.res = myFitSample(seg.obj = seg.dat, mut, Q = kQ, pi_theta_qz, sigma.h, 
                                    tau.dom = kTauDom, sigma.h.dom = kSigmaHDom, 
                           chr.arms.dat, verbose = verbose, d.res = d.res, b.res = b.res)
    mode.res = myapply_subclonal_scna_model(segobj = seg.dat, mode_res = mode.res, 
                                            verbose = verbose)
    if (inherits(mode.res, "try-error")) {
      mode.res = list(mode.flag = "FAIL")
    }
  }
  if (is.na(mode.res[["mode.flag"]])) {
    bad.ix = ABSOLUTE:::GenomeHetFilter(seg.dat[["obs.scna"]], mode.res, 
                                        max.non.clonal, max.neg.genome, kQ, verbose = verbose)
    if (sum(bad.ix) == nrow(mode.res[["mode.tab"]])) {
      mode.res = list(mode.flag = "ALPHA_TAU_DOM")
    } else {
      mode.res = ABSOLUTE:::ReorderModeRes(mode.res, !bad.ix)
    }
  }
  if (is.na(mode.res[["mode.flag"]])) {
    data(ChrArmPriorDb, package = "ABSOLUTE")
    model.id = ifelse(seg.dat[["group"]] %in% names(train.obj), 
                      seg.dat[["group"]], "Primary")
    mode.res = ABSOLUTE:::ApplyChrArmPrior(mode.res, model.id, train.obj, 
                                           verbose = verbose)
    if ((!is.null(maf)) && (nrow(maf) > 0)) {
      if (is.null(min.mut.af)) {
        stop("maf was not NULL, but min.mut.af was")
      }
      seg.dat[["mut.cn.dat"]] = mut.cn.dat
      mode.res = ABSOLUTE:::ApplySomaticMutsModel(mode.res, seg.dat$obs_SCNA, 
                                                  mut.cn.dat, kPiSomThetaQ, mut_class_w, kQ, verbose = verbose)
    }
    mode.res = ABSOLUTE:::WeighSampleModes(mode.res)
    mode.res[["call.status"]] = ABSOLUTE:::GetCallStatus(mode.res, seg.dat[["obs.scna"]][["W"]])
  }
  seg.dat[["mode.res"]] = mode.res
  if (is.null(output.fn.base)) {
    output.fn.base = ifelse(is.null(seg.dat$array.name), 
                            sample.name, seg.dat$array.name)
  }
  file.base = paste(output.fn.base, ".ABSOLUTE", sep = "")
  if (is.na(mode.res[["mode.flag"]])) {
    sample.pdf.fn = file.path(results.dir, paste(file.base, 
                                                 "plot.pdf", sep = "_"))
    ABSOLUTE:::AbsoluteResultPlot(sample.pdf.fn, seg.dat)
  } else {
    if (verbose) {
      print("Mode flag is NA, not generating plots. Sample has failed ABSOLUTE")
    }
  }
  seg.dat$version = 1.1
  save(seg.dat, file = file.path(results.dir, paste(file.base, 
                                                    "RData", sep = ".")))
  return(TRUE)
}

myFindLocationModes = function (obs, mut, Q, theta.qz, sigma.h, tau.dom, verbose = FALSE,
                                b.res = 0.125, d.res = 0.125) 
{
  kAlphaDom <- c(0, 1)
  kDom1 <- c(0, 1)
  kDom2 <- log(c(0.08, 1.05))
  lambda.qz.res <- ABSOLUTE:::GetLambda(theta.qz, obs[["W"]], Q + 1, verbose = verbose)
  pz <- theta.qz[Q + 1]
  mode.tab <- ABSOLUTE:::MargModeFinder(obs, mut, kDom1, kDom2, Q, lambda.qz.res, 
                                        pz, sigma.h, tau.dom, verbose = verbose, b.res = b.res,
                                        d.res = d.res)
  if (nrow(mode.tab) == 0) {
    return(list(mode.flag = "DELTA_B_DOM"))
  }
  mode.tab <- cbind(mode.tab, NA, NA, NA)
  for (i in 1:nrow(mode.tab)) {
    par <- mode.tab[i, c(1, 2)]
    if (obs[["platform"]] == "ARRAY") {
      cur.par = par
      if (is.null(obs$error.model$at)) {
        stop("at/AT issue, contact Jeff Gentry jgentry@broadinstitute.org")
      }
      mode.params = list(b = par[1], delta = par[2], AT = obs$error.model$at)
      obs$error.model$fit.at = obs$error.model$at
    }
    else {
      cur.par <- par
      mode.params <- list(b = par[1], delta = par[2], AT = NA)
    }
    LL <- ABSOLUTE:::CombLL(cur.par, Q = Q, obs = obs, dom1 = kDom1, 
                 dom2 = kDom2, lambda.qz.res, pz, sigma.h)
    mode.hess <- hessian(ABSOLUTE:::CombLL, x = cur.par, method = "Richardson", 
                         Q = Q, obs = obs, dom1 = kDom1, dom2 = kDom2, lambda.qz.res = lambda.qz.res, 
                         pz = pz, sigma.h = sigma.h)
    if (!is.na(mode.hess[1])) {
      mode.curv <- ABSOLUTE:::CalcModeLogCurv(par, mode.hess, verbose = verbose)
    }
    else {
      LL <- NA
      mode.curv <- NA
    }
    mode.tab[i, ] <- c(mode.params[["b"]], mode.params[["delta"]], 
                       mode.params[["AT"]], LL, mode.curv)
  }
  b <- mode.tab[, 1]
  delta <- exp(mode.tab[, 2])
  at <- mode.tab[, 3]
  res <- ABSOLUTE:::GetAlphaAndTau(b, delta)
  alpha <- res[["alpha"]]
  tau <- res[["tau"]]
  LL <- mode.tab[, 4]
  mode.curv <- mode.tab[, 5]
  mode.tab <- cbind(alpha, tau, at, b, delta, LL, mode.curv)
  colnames(mode.tab) <- c("alpha", "tau", "AT", "b", "delta", 
                          "LL", "mode_curv")
  if (verbose) {
    print(mode.tab)
  }
  mode.ix <- (mode.tab[, "alpha"] >= kAlphaDom[1] & mode.tab[, 
                                                             "alpha"] <= kAlphaDom[2] & mode.tab[, "tau"] >= tau.dom[1] & 
                mode.tab[, "tau"] <= tau.dom[2])
  if (verbose) {
    print(paste("removing ", sum(!mode.ix), " / ", length(mode.ix), 
                " modes outside of alpha/tau range.", sep = ""))
  }
  mode.tab <- mode.tab[mode.ix, , drop = FALSE]
  if (nrow(mode.tab) == 0) {
    return(list(mode.flag = "ALPHA_TAU_DOM"))
  }
  return(list(mode.tab = mode.tab))
}

myFitSample = function (seg.obj, mut, Q, pi_theta_qz, sigma.h, tau.dom, sigma.h.dom, 
                        chr.arms.dat, verbose = FALSE, b.res = 0.125, d.res = 0.125) 
{
  kThetaQz = rep(1/(Q + 1), Q + 1)
  obs = seg.obj[["obs.scna"]]
  res = myFindLocationModes(obs, mut, Q, kThetaQz, sigma.h, tau.dom, 
                            verbose = verbose, b.res = b.res, d.res = d.res)
    if (!is.null(res[["mode.flag"]])) {
    return(list(mode.flag = res[["mode.flag"]]))
  } else {
    mode.tab = res[["mode.tab"]]
  }
  if (is.na(mode.tab)) {
    return(list(mode.flag = "ERROR"))
  }
  n.modes = nrow(mode.tab)
  theta.qz.hat = matrix(NA, nrow = n.modes, ncol = Q + 1)
  theta.q.tab = array(NA, dim = c(n.modes, Q))
  if (verbose) {
    print(paste("Optimizing LL(data, theta.qz, sigma.h | comb) for ", 
                n.modes, " modes: ", sep = ""))
  }
  seg.z.tab = array(NA, dim = c(n.modes, length(obs[["W"]])))
  seg.qz.tab = array(NA, dim = c(n.modes, length(obs[["W"]]), 
                                 Q + 1))
  seg.q.tab = array(NA, dim = c(n.modes, length(obs[["W"]]), 
                                Q))
  if (obs[["data.type"]] == "ALLELIC") {
    ab.tab = array(NA, dim = c(n.modes, Q + 1))
    chr.arm.tab = array(NA, dim = c(n.modes, 2, nrow(chr.arms.dat), 
                                    Q))
  }
  if (obs[["data.type"]] == "TOTAL") {
    ab.tab = NULL
    chr.arm.tab = array(NA, dim = c(n.modes, 1, nrow(chr.arms.dat), 
                                    Q))
  }
  dimnames(chr.arm.tab)[[3]] = rownames(chr.arms.dat)
  mode.tab = cbind(mode.tab, NA, NA, NA, NA, NA, NA, NA, NA)
  colnames(mode.tab)[c((ncol(mode.tab) - 7):ncol(mode.tab))] = c("genome mass", 
                                                                 "sigma.h.hat", "theta.z.hat", "frac.het", "log.partition", 
                                                                 "entropy", "clust_LL", "post_LL")
  mode.tab = cbind(mode.tab, somatic.mut.ll = NA)
  for (i in 1:n.modes) {
    delta = mode.tab[i, "delta"]
    b = mode.tab[i, "b"]
    obs[["error.model"]][["fit.at"]] = mode.tab[i, "AT"]
    comb = GetCopyRatioComb(Q, delta, b, obs[["error.model"]])
    res = ABSOLUTE:::OptThetaQzSigmaH(obs, comb, sigma.h.dom, pi_theta_qz, 
                           verbose = verbose)
    mode.tab[i, "sigma.h.hat"] = res[["sigma.h.hat"]]
    mode.tab[i, "log.partition"] = res[["LL"]]
    mode.tab[i, "post_LL"] = res[["LL"]]
    lambda.qz.res = res[["lambda.qz.res"]]
    mode.tab[i, "theta.z.hat"] = res$theta_qz_hat[(Q + 1)]
    res = ABSOLUTE:::GetSegQzPost(obs[["d.tx"]], obs[["d.stderr"]], 
                       obs[["W"]], comb, NA, mode.tab[i, "sigma.h.hat"], 
                       theta.qz = res[["theta_qz_hat"]])
    seg.z.tab[i, ] = res[["QZ"]][, (Q + 1)]
    seg.qz.tab[i, , ] = res[["QZ"]]
    seg.q.tab[i, , ] = res[["Q"]]
    theta.q.tab[i, ] = colSums(res[["Q"]] * obs[["W"]])
    mode.tab[i, "theta.z.hat"] = sum(seg.qz.tab[i, , (Q + 
                                                        1)] * obs$W)
    mode.tab[i, "frac.het"] = sum(obs[["W"]] * seg.z.tab[i, 
                                                         ])
    if (obs[["data.type"]] == "ALLELIC") {
      ab.tab[i, ] = ABSOLUTE:::CalcAbDistr(obs, res[["QZ"]])
      mode.tab[i, "genome mass"] = 2 * sum(c((1:Q) - 1) * 
                                             colSums(res[["Q"]] * obs[["W"]]))
    }
    if (obs[["data.type"]] == "TOTAL") {
      mode.tab[i, "genome mass"] = 1 * sum(c((1:Q) - 1) * 
                                             colSums(res[["Q"]] * obs[["W"]]))
    }
    mode.tab[i, "entropy"] = ABSOLUTE:::CalcFitEntropy(obs, res[["QZ"]])
    chr.arm.tab[i, , , ] = CalcChrArmDistr(seg.obj, res[["Q"]], 
                                           chr.arms.dat)
  }
  return(list(mode.tab = mode.tab, mode.posts = NA, theta.q.tab = theta.q.tab, 
              b = theta.qz.hat, seg.z.tab = seg.z.tab, seg.qz.tab = seg.qz.tab, 
              seg.q.tab = seg.q.tab, ab.tab = ab.tab, chr.arm.tab = chr.arm.tab, 
              mode.flag = NA))
}
myCreateReviewObject = function (obj.name, absolute.files, indv.results.dir, copy_num_type, 
                                 plot.modes = TRUE, verbose = FALSE, myModes = NULL) 
{
  if (copy_num_type == "total") {
    ABSOLUTE:::set_total_funcs()
  } else if (copy_num_type == "allelic") {
    ABSOLUTE:::set_allelic_funcs()
  } else {
    stop("Unsupported copy number type: ", copy_num_type)
  }
  dir.create(indv.results.dir, recursive = TRUE)
  nm = file.path(indv.results.dir, paste(obj.name, ".PP-modes", 
                                         sep = ""))
  modesegs.fn = paste(nm, ".data.RData", sep = "")
  pdf.fn = paste(nm, ".plots.pdf", sep = "")
  failed.pdf.fn = paste(nm, "FAILED_plots.pdf", sep = ".")
  failed.tab.fn = paste(nm, "FAILED_tab.txt", sep = ".")
  call.tab.fn = file.path(indv.results.dir, paste(obj.name, 
                                                  "PP-calls_tab.txt", sep = "."))
  segobj.list = vector(mode = "list", length = length(absolute.files))
  failed.list = vector(mode = "list", length = length(absolute.files))
  so.ix = 1
  fa.ix = 1
  sids = character()
  for (i in seq_along(absolute.files)) {
    seg.out.fn = absolute.files[i]
    if (!file.exists(seg.out.fn)) {
      if (verbose) {
        cat("\n")
        print(paste("sample #", i, " result not found", 
                    sep = ""))
      }
      next
    }
    load(absolute.files[i])
    SID = seg.dat$sample.name
    if (SID %in% sids) {
      stop("Sample names must be unique")
    }
    sids = c(sids, SID)
    if (is.null(seg.dat$array.name)) {
      seg.dat$array.name = SID
    }
    if (is.na(seg.dat[["mode.res"]][["mode.flag"]])) {
      segobj.list[[so.ix]] = seg.dat
      names(segobj.list)[so.ix] = SID
      so.ix = so.ix + 1
      if (verbose) {
        cat(".")
      }
    }    else {
      if (verbose) {
        print(paste("Failed list got updated", i))
      }
      failed.list[[fa.ix]] = seg.dat
      names(failed.list)[fa.ix] = SID
      failed.list[[fa.ix]][["sample.name"]] = seg.dat[["sample.name"]]
      fa.ix = fa.ix + 1
      if (verbose) {
        cat("-")
      }
    }
  }
  segobj.list = segobj.list[c(1:(so.ix - 1))]
  failed.list = failed.list[c(1:(fa.ix - 1))]
  mode.ent = rep(NA, length(segobj.list))
  names(mode.ent) = names(segobj.list)
  for (i in seq_along(segobj.list)) {
    mtab = segobj.list[[i]][["mode.res"]][["mode.tab"]]
    ix = which.max(mtab[, "post_LL"])
    mode.ent[i] = mtab[ix, "entropy"]
  }
  samples = names(sort(mode.ent))
  segobj.list = segobj.list[samples]
  save(segobj.list, file = modesegs.fn)
  ABSOLUTE:::PrintPpCallTable(segobj.list, call.tab.fn)
  if (plot.modes) {
    n.print = length(myModes)
    myPlotModes(segobj.list, pdf.fn, n.print = n.print, myModes = myModes)
  }
  if (!is.null(failed.list[[1]])) {
    ABSOLUTE:::PlotFailedSamples(failed.list, failed.pdf.fn)
    ABSOLUTE:::PrintFailedTable(failed.list, failed.tab.fn)
  }
  return(TRUE)
}



myPlotModes = function (segobj.list, pdf.fn, n.print = NA, write.tab = TRUE, 
                        debug.info = TRUE, add = FALSE, fig.mode = FALSE, sideways = FALSE, 
                        low.w = FALSE, called.mode.ix = NA, verbose = FALSE, myModes = myModes) 
{
  pal <- colorRampPalette(c("chocolate4", "cyan4"))
  if (fig.mode) {
    add <- TRUE
    debug.info <- FALSE
  }
  Q <- 15
  alpha.dom <- c(0, 1)
  tau.dom <- c(1, 10)
  mode.colors <- ABSOLUTE:::GetModeColors()
  x.max <- 2.25
  binW <- 0.025
  if (!add) {
    pdf(pdf.fn, 3.5 * (1 + 2), 15)
    par(mfrow = c(5, 4))
  }
  for (s in seq_len(length(segobj.list))) {
    mode.tab <- segobj.list[[s]][["mode.res"]][["mode.tab"]]
    seg.z.tab <- segobj.list[[s]][["mode.res"]][["seg.z.tab"]]
    if (nrow(mode.tab) == 0) {
      for (i in seq_len(n.print + 2)) {
        frame()
      }
      next
    }
    obs <- ExtractSampleObs(segobj.list[[s]])
    SN <- segobj.list[[s]][["sample.name"]]
    if (is.na(n.print)) {
      s.n.print <- nrow(mode.tab)
    }
    else {
      s.n.print <- min(nrow(mode.tab), n.print)
    }
    s.mix <- s.n.print
    if (!fig.mode) {
      model.id <- segobj.list[[s]][["group"]]
      ABSOLUTE:::PostPlot(mode.tab[myModes,], mode.colors, alpha.dom, tau.dom, 
                          sample.name = SN, obs, called.mode.ix, debug.info, call.status = 
                            segobj.list[[s]][["mode.res"]][["call.status"]], 
                          model.id = model.id)
      ABSOLUTE:::PpModeScorerBarplot(mode.tab[myModes,], mode.colors, obs, n.print)
      if (length(segobj.list) == 1) {
        frame()
        frame()
      }
    }
    n.plot <- min(n.print, NROW(mode.tab))
    for (i in seq_len(n.plot)) {
      tau <- mode.tab[myModes[i], "tau"]
      alpha <- mode.tab[myModes[i], "alpha"]
      delta <- alpha/(2 * (1 - alpha) + alpha * tau)
      b <- 2 * (1 - alpha)/(2 * (1 - alpha) + alpha * tau)
      mode.info <- mode.tab[myModes[i], ]
      obs[["error.model"]][["fit.at"]] <- mode.tab[myModes[i], "AT"]
      comb <- InvTxData(obs[["error.model"]], GetCopyRatioComb(Q, 
                                                               delta, b, obs[["error.model"]]))
      seg.dat <- InvTxData(obs[["error.model"]], obs[["d.tx"]])
      last.plot <- ifelse(i == n.plot, TRUE, FALSE)
      mid.plot <- ifelse(i == ceiling(n.plot/2), TRUE, 
                         FALSE)
      first.plot <- ifelse(i == 1, TRUE, FALSE)
      if ((length(segobj.list) > 1) && (!is.null(segobj.list[[s]][["mut.cn.dat"]])) && 
            (i > 1)) {
        frame()
        frame()
      }
      if (!is.null(segobj.list[[s]][["mode.res"]][["ab.tab"]])) {
        comb.ab <- segobj.list[[s]][["mode.res"]][["ab.tab"]][myModes[i], ]
      }
      else {
        comb.ab <- NA
      }
      ABSOLUTE:::ModeCombPlot(seg.dat, obs[["W"]], obs[["d.stderr"]], 
                              mode.info, seg.z.tab[myModes[i], ], comb.ab, mode.colors[i], 
                              comb, last.plot, mid.plot, first.plot, x.max, 
                              debug.info, fig.mode, sideways, pal)
      if (!is.na(called.mode.ix)) {
        if (i == called.mode.ix) {
          mtext(text = "*", side = 3, at = 0, col = "red", 
                line = -1, cex = 3 * par("cex") * par("cex.axis"))
        }
      }
      if (!is.null(segobj.list[[s]][["mut.cn.dat"]])) {
        mut.cn.dat <- segobj.list[[s]][["mut.cn.dat"]]
        modeled <- segobj.list[[s]][["mode.res"]][["muts.post.prs"]][, 
                                                                     , myModes[i], drop = TRUE]
        if (is.null(dim(modeled))) {
          modeled <- matrix(modeled, nrow = 1)
          colnames(modeled) <- dimnames(segobj.list[[s]][["mode.res"]][["muts.post_prs"]])[[2]]
        }
        modeled.mut.dat <- cbind(mut.cn.dat, modeled)
        ABSOLUTE:::PlotSomaticMutDensities(modeled.mut.dat, segobj.list[[s]], 
                                           1, mode.colors[i], min.cov = 3, verbose = verbose)
      }
    }
    if ((length(segobj.list) > 1) && (is.null(segobj.list[[s]][["mut.cn.dat"]]))) {
      if (i < n.print) {
        for (j in seq_len(n.print - i)) {
          frame()
        }
      }
      frame()
    }
  }
  if (!add) {
    dev.off()
  }
}



myRunHapSeg = function (plate.name, array.name, seg.fn, snp.fn, genome.build, 
                        results.dir, platform, use.pop, impute.gt, plot.segfit, merge.small, 
                        merge.close, min.seg.size, normal, out.p, seg.merge.thresh, 
                        use.normal, adj.atten, phased.bgl.dir, force.diploid = normal, 
                        drop.x = FALSE, drop.y = TRUE, calls.fn = NULL, mn.sample = NULL, 
                        out.file = NULL, calibrate.data = NULL, clusters.fn = NULL, 
                        prev.theta.fn = NULL, snp.file.parser = DefaultSnpFileParser, 
                        clusters.file.parser = DefaultClustersFileParser, verbose = FALSE) 
{
  if (!HAPSEG:::CheckGenomeBuild(genome.build)) {
    stop("Unsupported genome build: ", genome.build)
  }
  if ((!impute.gt) && (verbose)) {
    print("Note: impute.gt is set to FALSE. It is recommended that this be set to TRUE")
  }
  if (!identical(adj.atten, FALSE)) {
    print("Note: The adj.atten argument has been disabled and will be removed in a future release")
  }
  tmp.dir = HAPSEG:::CreateTmpDir(results.dir)
  on.exit(try(unlink(tmp.dir, recursive = TRUE), silent = TRUE), 
          add = TRUE)
  iams.res = HAPSEG:::InitAndMergeSmall(array.name, genome.build, use.pop, 
                                        use.normal, normal, impute.gt, platform, seg.fn, snp.fn, 
                                        drop.x, drop.y, calls.fn, mn.sample, min.seg.size, merge.close, 
                                        out.p, merge.small, force.diploid, calibrate.data, clusters.fn, 
                                        prev.theta.fn, snp.file.parser = snp.file.parser, clusters.file.parser = clusters.file.parser, 
                                        verbose = verbose)
  if (merge.close == TRUE) {
    mcaf.res = HAPSEG:::MergeCloseAndFit(iams.res, out.p, seg.merge.thresh, 
                                         impute.gt, platform, force.diploid = force.diploid, 
                                         verbose = verbose)
  }
  else {
    if (verbose) {
      print("Merged Loci: nulled ")
    }
    seg.dat[["merged.loci"]] = NA
    seg.dat[["final.merge.prob"]] = matrix(NA, ncol = 2, 
                                           nrow = length(iams.res[["h.d"]]) - 1)
    seg.dat[["seg.expected.phase"]] = iams.res[["em.fit"]][["seg.expected.phase"]]
    mcaf.res = iams.res
    mcaf.res[["seg.dat"]] = seg.dat
  }
  if (impute.gt == TRUE) {
    gt.res = HAPSEG:::ImputeGt(mcaf.res, platform, tmp.dir, plate.name, 
                               out.p, phased.bgl.dir, verbose = verbose)
    gt.res = HAPSEG:::PostImputeSegFit(gt.res, out.p, platform, force.diploid = force.diploid, 
                                       verbose = verbose)
  }
  else {
    gt.res = mcaf.res
  }
  seg.dat = HAPSEG:::FinishHapsegAndSave(gt.res, plate.name, normal, 
                                         array.name, results.dir, platform, out.file)
  if (plot.segfit) {
    if (nrow(seg.dat$merged.loci) == 0) seg.dat$merged.loci = NA
    HAPSEG:::DoPlots(seg.dat, results.dir, platform, verbose = verbose)
  }
  return(TRUE)
}

ABSOLUTE:::OptThetaQzSigmaH
function (obs, comb, sigma.h.dom, pi_theta_qz, verbose = FALSE) 
{
  objective <- function(par, obs, comb, lambda) {
    sigma.h <- par
    max.q <- length(comb)
    theta.qz.hat <- ABSOLUTE:::GetThetaQzPost(obs[["d.tx"]], obs[["d.stderr"]], 
                                   obs[["W"]], comb, sigma.h, pi_theta_qz)
    theta.qz.map <- theta.qz.hat
    theta.qz.map <- theta.qz.map/sum(theta.qz.map)
    theta.qz.hat <- theta.qz.map
    if (verbose) {
      cat("S")
    }
    LL <- ABSOLUTE:::CalcNormLoglik(obs[["d.tx"]], obs[["d.stderr"]], 
                         obs[["W"]], comb, NA, theta.qz.hat, sigma.h)
    if (is.nan(LL)) {
      LL <- -Inf
    }
    return(LL)
  }
  max_q = length(comb)
  use.sigma.h <- mean(sigma.h.dom)
  theta.qz.hat <- ABSOLUTE:::GetThetaQzPost(obs[["d.tx"]], obs[["d.stderr"]], 
                                 obs[["W"]], comb, use.sigma.h, pi_theta_qz)
  lambda.qz.res <- ABSOLUTE:::GetLambda(theta.qz.hat, obs[["W"]], max.q + 
                               1, verbose = verbose)
  init.lambda <- lambda.qz.res[["Lambda"]]
  res <- optimize(f = objective, interval = sigma.h.dom, tol = 0.001, 
                  maximum = TRUE, obs = obs, comb = comb, lambda = init.lambda)
  sigma.h.hat <- res[["maximum"]]
  theta.qz.hat <- ABSOLUTE:::GetThetaQzPost(obs[["d.tx"]], obs[["d.stderr"]], 
                                 obs[["W"]], comb, sigma.h.hat, pi_theta_qz)
  lambda.qz.res <- ABSOLUTE:::GetLambda(theta.qz.hat, obs[["W"]], max.q + 
                               1, verbose = verbose)
  LL <- ABSOLUTE:::CalcNormLoglik(obs[["d.tx"]], obs[["d.stderr"]], obs[["W"]], 
                       comb, lambda.qz.res, NA, sigma.h.hat)
  if (is.nan(LL)) {
    LL <- -Inf
    if (verbose) {
      print("Warning: NaN loglik in opt_theta_Z_sigma_H")
    }
  }
  if (verbose) {
    cat(paste("\tsigma.h.hat=", round(sigma.h.hat, 5), sep = ""))
    cat("\n")
  }
  return(list(LL = LL, sigma.h.hat = sigma.h.hat, theta_qz_hat = theta.qz.hat, 
              lambda.qz.res = lambda.qz.res))
}


