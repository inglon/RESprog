################################################################
################################################################
# Script: ProgR_5.4_Heterogeneity.R
# Author: Ingrid Lonnstedt
# Date:  23/05/2013
# R version: R 3.0
# Details: Derive basic heterogeneity estimates
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

prefix.raw = paste(getwd(), "/AROMA/rawData/responsify/GenomeWideSNP_6/", sep='')
prefix.ann = paste(getwd(), "/AROMA/annotationData/chipTypes/GenomeWideSNP_6/", sep='')
prefix.out = paste(getwd(), "/ABSOLUTE/HAPSEG/", sep='')


date = format(Sys.Date())
source(file.path(getwd(), 'programs', 'ProgR_0.1_Functions.R'))



################################################################
################################################################
#
# Format intensities for HAPSEG
#
################################################################
################################################################
info = read.delim(paste(prefix.ann, '../../Infokey.txt', sep = ''), stringsAsFactors = F)
info$sample = substr(info$Sample.ID, 1, 4)
info$sample.name = paste(rownames(info), info$sample, sep = '.')
info$myarray = info$arrname

######Find segments with 2m1 and 2m0

pdf(paste(prefix.out, 'Find_segments_with_2m1_2m0_', date, '.pdf', sep=''), paper = 'a4r',
    width = 11, height = 8) 
par(mfrow=c(2,2), mar = c(2,2,3,0)+.1)
for (i in 1:56){
ta = read.delim(file.path(getwd(), 'TAPS', info$sample.name[i], 
                          paste(info$sample.name[i], 'segmentCN.txt', sep = '_')))
c2 = subset(ta, Cn == 2 & !is.na(imba))
plot(imba~log2, data = ta, xlim = c(-1,2), main = info$sample.name[i])
points(imba~log2, data = c2, col = 'red')
tmp = sort(kmeans(c2$imba, 2)$centers)
het = tmp[1]
loh = tmp[2]
mid = (het + loh)/2  

loh.segs = subset(c2, imba>mid)
if (nrow(loh.segs)>1){
den = density(loh.segs$imba)
loh.r = range(den$x[which(den$y>max(den$y)*2/3)])
loh.segs = subset(c2, imba>loh.r[1] & imba<loh.r[2], select = c('Chromosome','Start','End'))
points(imba~log2, data=c2, subset = imba>loh.r[1] & imba<loh.r[2], col = 'green')
} else {
  loh.segs = subset(c2, select = c('Chromosome','Start','End'))
  points(imba~log2, data=c2, subset = imba>mid, col = 'green')
}

het.segs = subset(c2, imba<mid)
den = density(het.segs$imba)
het.r = range(den$x[which(den$y>max(den$y)/2)])
het.segs = subset(c2, imba>het.r[1] & imba<het.r[2], select = c('Chromosome','Start','End'))
points(imba~log2, data=c2, subset = imba>het.r[1] & imba<het.r[2], col = 'blue')

d = file.path(getwd(), 'ABSOLUTE', 'HAPSEG', info$sample.name[i])
dir.create(d, recursive= T, showWarnings = F)
write.table(het.segs, file.path(d, 'Hetsegs.txt'), row.names = F, quote = F, sep = '\t')
write.table(loh.segs, file.path(d, 'Homsegs.txt'), row.names = F, quote = F, sep = '\t')
}
dev.off()
     
######Find allele intensities

for (i in 1:56){   
library(aroma.affymetrix)
dataSet <- "responsify";
tags <- "ACC,ra,-XY,BPN,-XY,RMA,FLN,-XY,CMTN,v2";
chipType <- "GenomeWideSNP_6";
ds <- AromaUnitTotalCnBinarySet$byName(dataSet, tags=tags, chipType=chipType);
df <- getFile(ds, which(getNames(ds) == info$arrname2[i]));
bs <- AromaUnitFracBCnBinarySet$byName(dataSet, tags=tags, chipType=chipType);
bf <- getFile(bs, which(getNames(ds) == info$arrname2[i]));
ugp <- getAromaUgpFile(df)
unf <- getUnitNamesFile(ugp)
unitNames <- getUnitNames(unf)
chromosome <- ugp[, 1, drop = TRUE]
position <- ugp[, 2, drop = TRUE]
data = data.frame(unitNames = unitNames, chromosome=chromosome, position = position, 
                  total = df[,1], fracB = bf[,1], stringsAsFactors = F)
names(data)[4:5] = c('total','fracB')
d = file.path(getwd(), 'ABSOLUTE', 'HAPSEG', info$sample.name[i])
dir.create(d, recursive= T, showWarnings = F)
write.table(data, file.path(d, 'Total,fracB.txt'), row.names = F, quote = F, sep = '\t')
}


######Find calibration constants
i = 24
i = 21
i = 26

library(IRanges)
pdf(paste(prefix.out, 'Find_calibration_constants_', date, '.pdf', sep=''), paper = 'a4') 
for (i in 1:56){
  d = file.path(getwd(), 'ABSOLUTE', 'HAPSEG', info$sample.name[i])
  hets = read.delim(file.path(d, 'Hetsegs.txt'), stringsAsFactors = F)
  data = read.delim(file.path(d, 'Total,fracB.txt'), stringsAsFactors = F)
  
  nc = nchar(hets$Chromosome, type = 'c')
  hets$Chromosome = substr(hets$Chromosome, 4, nc)
  hets$Chromosome[hets$Chromosome == 'X'] = '23'
  hets$Chromosome[hets$Chromosome == 'Y'] = '24'
  hets$Chromosome[hets$Chromosome == 'M'] = '25'
  hets$Chromosome = as.integer(hets$Chromosome)
  
  subd = subset(data, !is.na(chromosome) & !is.na(position) & !is.na(fracB) & !is.na(total) & !is.na(unitNames))
  subd$cn2 = 'no'
  subd$chromosome = factor(subd$chromosome)
  hets$Chromosome = factor(hets$Chromosome, levels = levels(subd$chromosome))
  h = RangedData(IRanges(start = hets$Start, end = hets$End),  
                 space = hets$Chromosome, universe = 'hg19') 
  d = RangedData(IRanges(start = subd$position, width = 1),  
                 space = subd$chromosome, universe = 'hg19')
  hits = findOverlaps(query = h, subject = d)
  for (ch in levels(subd$chromosome)){
    myhits = subjectHits(hits[[ch]])
    ix = which(subd$chromosome == ch)
    subd$cn2[ix[myhits]] = 'het' 
  }
  
  subd = subset(subd, cn2 %in% c('het'))
#  subd$baf = pmax(subd$fracB, 0)
#  subd$baf = pmin(subd$baf, 1)
  subd$baf = subd$fracB
  
  subd$a = subd$total*(1-subd$baf)
  subd$b = subd$total*subd$baf
  
  par(mfrow=c(1,2))

  X = subd$a
  q = quantile(X, c(.05, .95))
  X = X[X> q[1] & X<q[2]]
  den = density(X)
#  plot(den, main = 'A intensities')
  hst = hist(X, breaks = 500, freq = F, xlab = 'Allele A intensities',
             main = '')
  lines(den, col = 'red', lwd = 2)
  A0hist = hst$mids[which.max(hst$density)] #[1] 0.0025
  A0med = median(X[X<.5]) #[1] 0.04358323
  A0kernel = den$x[which.max(den$y)]  #[1] 0.03768932
  A0max = quantile(X[X<.5], probs = .975) #0.3266564 
  
  X = subd$b
  q = quantile(X, c(.05, .95))
  X = X[X> q[1] & X<q[2]]
  den = density(X)
#  plot(den, main = 'B intensities')
  hst = hist(X, breaks = 500, freq = F, xlab = 'Allele B intensities',
             main = '')
  lines(den, col = 'red', lwd = 2)
  B0hist = hst$mids[which.max(hst$density)] #[1] 0.0025
  B0med = median(X[X<.5]) #[1] 0.06436947
  B0kernel = den$x[which.max(den$y)]  #[1] 0.05391259
  B0max = quantile(X[X<.5], probs = .975) #0.3625948 
  d = file.path(getwd(), 'ABSOLUTE', 'HAPSEG', info$sample.name[i])
  dev.print(png, file=file.path(d, paste(info$sample.name[i], '_Calibration.png', 
                                         sep = '')), 
            width=640, height=540)
  
  #A1
  x = subd$a[subd$baf>.3 & subd$baf<0.7 & subd$cn2 == 'het']
  q = quantile(x, c(.05, .95))
  x = x[x> q[1] & x<q[2]]
  if (length(x)>1) {
    den = density(x, bw = 0.05)
    plot(den, main = 'A1: A intensities at AB')
    A1 = den$x[which.max(den$y)]
  } else {
    A1 = ifelse(length(x)>0, x, 1)
    
  }
  
  #B1
  x = subd$b[subd$baf>.3 & subd$baf<0.7 & subd$cn2 == 'het']
  q = quantile(x, c(.05, .95))
  x = x[x> q[1] & x<q[2]]
  if (length(x)>1) {
    den = density(x, bw = 0.05)
    plot(den, main = 'B1: B intensities at AB')
    B1 = den$x[which.max(den$y)]
  } else {
    B1 = ifelse(length(x)>0, x, 1)
    
  }
  
  title(main = info$sample.name[i], outer = T, line = -1)
  
  ###Standardize the A and B intensities and write to file for HAPSEG
  A0 = 0
  B0 = 0
  data$baf = pmax(data$fracB, 0)
  data$baf = pmin(data$fracB, 1)
  data$a = data$total*(1-data$baf)
  data$b = data$total*(data$baf)
  data$A = (data$a-A0)/(A1-A0)
  data$B = (data$b-B0)/(B1-B0)
  
  dat = subset(data, select = c('A','B'))
  rownames(dat) = data$unitNames
  d = file.path(getwd(), 'ABSOLUTE', 'HAPSEG', info$sample.name[i])
  #write.table(data, file.path(d, 'snpfile.txt'), quote = F, sep = '\t')
  save(dat, file = file.path(d, 'snpfile.RData'))
}
dev.off()

#No calibration
data$A = data$a
data$B = data$b
dat = subset(data, select = c('A','B'))
rownames(dat) = data$unitNames
d = file.path(getwd(), 'ABSOLUTE', 'HAPSEG', info$sample.name[i])
save(dat, file = file.path(d, 'snpfileNoCalibration.RData'))

#histmax calibration
data$baf = pmax(data$fracB, 0)
data$baf = pmin(data$fracB, 1)
data$a = data$total*(1-data$baf)
data$b = data$total*(data$baf)
data$A = (data$a-A0hist)/(A1-A0hist)
data$B = (data$b-B0hist)/(B1-B0hist)
dat = subset(data, select = c('A','B'))
rownames(dat) = data$unitNames
d = file.path(getwd(), 'ABSOLUTE', 'HAPSEG', info$sample.name[i])
save(dat, file = file.path(d, 'snpfileCalibration5.RData'))

#Kernel calibration
#data$baf = pmax(data$fracB, 0)
#data$baf = pmin(data$fracB, 1)
#data$a = data$total*(1-data$baf)
#data$b = data$total*(data$baf)
data$A = (data$a-A0kernel)/(A1-A0kernel)
data$B = (data$b-B0kernel)/(B1-B0kernel)
dat = subset(data, select = c('A','B'))
rownames(dat) = data$unitNames
d = file.path(getwd(), 'ABSOLUTE', 'HAPSEG', info$sample.name[i])
save(dat, file = file.path(d, 'snpfileCalibration6.RData'))


#Median calibration
#data$baf = pmax(data$fracB, 0)
#data$baf = pmin(data$fracB, 1)
#data$a = data$total*(1-data$baf)
#data$b = data$total*(data$baf)
data$A = (data$a-A0med)/(A1-A0med)
data$B = (data$b-B0med)/(B1-B0med)
dat = subset(data, select = c('A','B'))
rownames(dat) = data$unitNames
d = file.path(getwd(), 'ABSOLUTE', 'HAPSEG', info$sample.name[i])
save(dat, file = file.path(d, 'snpfileCalibration7.RData'))


#Max calibration
#data$baf = pmax(data$fracB, 0)
#data$baf = pmin(data$fracB, 1)
#data$a = data$total*(1-data$baf)
#data$b = data$total*(data$baf)
data$A = (data$a-A0max)/(A1-A0max)
data$B = (data$b-B0max)/(B1-B0max)
dat = subset(data, select = c('A','B'))
rownames(dat) = data$unitNames
d = file.path(getwd(), 'ABSOLUTE', 'HAPSEG', info$sample.name[i])
save(dat, file = file.path(d, 'snpfileCalibration8.RData'))


############################
#Example figure for sample 15
#A1
x = subd$a[subd$cn2 == 'het']
A = x
q = quantile(x, c(.05, .95))
x = x[x> q[1] & x<q[2]]
den = density(x)
plot(den, main = 'Allele A normalized intensities in CN=2 normal segments',
     xlab = 'SNP array intensity', cex.main = .8)
#mtext('What is A1: A intensity when there is exactly one copy of A?', cex = .8)
den$x[which.max(den$y)]
[1] 0.05578445

#B1
x = subd$b[ subd$cn2 == 'het']
B=x
q = quantile(x, c(.05, .95))
x = x[x> q[1] & x<q[2]]
den = density(x)
plot(den, main = 'Allele B normalized intensities in CN=2 normal segments',cex.main = .8,
     xlab = 'SNP array intensity')
#mtext('What is B1: B intensity when there is exactly one copy of B?', cex = .8)
den$x[which.max(den$y)]
[1] 0.05214338

par(mfrow = c(1,1))
smoothScatter(A, B, xlim = c(0,2), ylim = c(0,2), xlab = 'Allele A intensities',
              ylab = 'Allele B intensities', 
              main = 'Normalized SNP array intensities for normal CN=2 segments')
d = file.path(getwd(), 'ABSOLUTE', 'HAPSEG', info$sample.name[i])
dev.print(png, file=file.path(d, paste(info$sample.name[i], '_Chosen_subclones.png', sep = '')), 
          width=640, height=640)




###### Segment Logratios for calibrated data

library(aroma.affymetrix)
tot = NULL
for (i in 1:56){
  d = file.path(getwd(), 'ABSOLUTE', 'HAPSEG', info$sample.name[i])
  load(file.path(d, 'snpfile.RData'))
  tot = cbind(tot, apply(dat, 1, sum))
  rm(dat)
}
den = apply(tot, 1, median, na.rm = T)

dataSet <- "responsify";
tags <- "ACC,ra,-XY,BPN,-XY,RMA,FLN,-XY,CMTN,v2";
chipType <- "GenomeWideSNP_6";
ds <- AromaUnitTotalCnBinarySet$byName(dataSet, tags=tags, chipType=chipType);
df <- getFile(ds, which(getNames(ds) == info$arrname2[1]));
ugp <- getAromaUgpFile(df)
chromosome <- ugp[, 1, drop = TRUE]
position <- ugp[, 2, drop = TRUE]

data = data.frame(chromosome=chromosome, position = position, 
                  stringsAsFactors = F)

library(PSCBS)
library(IRanges)
for (i in 1:56){
  print(i)
  #Prepare data for this sample
  fitname = paste('fit', i, sep = '')
  filename <- paste(getwd(), '/output/SNP/CBSfits/', fitname, '.RData', sep = '');
  fit <- loadObject(filename);
  seg = getSegments(fit)
  seg$Sample = info$myarray[i]
  seg = data.frame(Sample = seg$Sample, Chromosome = seg$chromosome, "Start Position" = seg$tcnStart,
                   "End Position" = seg$tcnEnd, "Num markers" = seg$tcnNbrOfSNPs, 
                   check.names = F)
  seg = subset(seg, !is.na(Chromosome) & Chromosome != 25)
  seg$Chromosome = factor(seg$Chromosome, levels = as.character(1:24))
  seg$cn = NA
  d = file.path(getwd(), 'ABSOLUTE', 'HAPSEG', info$sample.name[i])
  load(file.path(d, 'snpfile.RData'))
  data$cn = (apply(dat, 1, sum)/den)
  ix = !is.na(data$chromosome) & !is.na(data$position) & !is.na(data$cn) & data$chromosome != 25
  h = RangedData(IRanges(start = seg[,'Start Position'], end = seg[,'End Position']),  
                 space = seg$Chromosome, universe = 'hg19') 
  d = RangedData(IRanges(start = data$position[ix], width = 1),  
                 space = data$chromosome[ix], universe = 'hg19', cn = data$cn[ix])
  hits = findOverlaps(query = h, subject = d)
  
  for (ch in levels(seg$Chromosome)){
    cn = d$cn[space(d) == as.numeric(ch)]
    seg$cn[seg$Chromosome == ch & seg[,'Num markers']>0] = tapply(cn[subjectHits(hits[[ch]])], queryHits(hits[[ch]]), 
                                          median, na.rm = T)
  }
  seg$cn = log2(seg$cn)
  d = file.path(getwd(), 'ABSOLUTE', 'HAPSEG', info$sample.name[i], 'Segdat.txt')
  write.table(seg, d, row.names = F, quote = F, sep = '\t')
}



################################################################
################################################################
#
# HAPSEG
#
################################################################
################################################################
library(HAPSEG)
#source(paste(getwd(), '/programs/Hapseg.r', sep = ''))
library(foreach)
library(doMC)
registerDoMC(20)
info = read.delim(paste(prefix.ann, '../../Infokey.txt', sep = ''))
info$sample = substr(info$Sample.ID, 1, 4)
info$myarray = info$arrname
info$sample.name = paste(info$ind, info$sample, sep = '.')

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

DoHapseg <- function(myarray, mysample) {
  registerDoSEQ()
  plate.name <- mysample
  phased.bgl.dir <- "/usr/local/bioinfsoftware/R/current/lib64/R/library/ABSOLUTE/etc/phasedBGL/hg19"
  d = file.path(getwd(), 'ABSOLUTE', 'HAPSEG', mysample)
  snp.fn <- file.path(d, 'snpfile.RData')
  seg.fn = file.path(d, 'Segdat.txt')
  results.dir <- paste(prefix.out, mysample, sep = '')
  if (!file.exists(results.dir)) {
    dir.create(results.dir, recursive=TRUE)
  }
  print(paste("Starting array", mysample, "at", results.dir))
  log.dir <- paste(prefix.out, 'log', sep = '')
  if (!file.exists(log.dir)) {
    dir.create(log.dir, recursive=TRUE)
  }
  sink(file=file.path(log.dir, paste(mysample, ".hapseg.out.txt", sep="")))
  myRunHapSeg(plate.name, myarray, seg.fn, snp.fn, genome.build = 'hg19', results.dir, platform = "SNP_6.0", 
            use.pop='CEPH', 
            impute.gt=TRUE, plot.segfit=TRUE, merge.small=FALSE, merge.close=TRUE, min.seg.size=2, 
            normal=FALSE, 
            out.p=0.001, seg.merge.thresh=1e-10, use.normal = FALSE, adj.atten = FALSE, 
            phased.bgl.dir, calibrate.data=FALSE, 
            drop.x=FALSE, verbose=TRUE)
  sink()
}
foreach (i=1:56, .combine=c) %dopar% {
  DoHapseg(myarray = info$myarray[i], mysample = info$sample.name[i])
}

DoHapseg2 <- function(myarray, mysample) {
  registerDoSEQ()
  plate.name <- mysample
  phased.bgl.dir <- "/usr/local/bioinfsoftware/R/current/lib64/R/library/ABSOLUTE/etc/phasedBGL/hg19"
  d = file.path(getwd(), 'ABSOLUTE', 'HAPSEG', mysample)
  snp.fn <- file.path(d, 'snpfile.RData')
  seg.fn = file.path(d, 'Segdat.txt')
  results.dir <- file.path(getwd(), 'ABSOLUTE', 'HAPSEG2', mysample)
  if (!file.exists(results.dir)) {
    dir.create(results.dir, recursive=TRUE)
  }
  print(paste("Starting array", mysample, "at", results.dir))
  log.dir <- file.path(getwd(), 'ABSOLUTE', 'HAPSEG2', 'log')
  if (!file.exists(log.dir)) {
    dir.create(log.dir, recursive=TRUE)
  }
  sink(file=file.path(log.dir, paste(mysample, ".hapseg.out.txt", sep="")))
  myRunHapSeg(plate.name, myarray, seg.fn, snp.fn, genome.build = 'hg19', results.dir, platform = "SNP_6.0", 
              use.pop='CEPH', 
              impute.gt=TRUE, plot.segfit=TRUE, merge.small=FALSE, merge.close=TRUE, min.seg.size=2, 
              normal=FALSE, 
              out.p=0.001, seg.merge.thresh=1e-80, use.normal = FALSE, adj.atten = FALSE, 
              phased.bgl.dir, calibrate.data=FALSE, 
              drop.x=FALSE, verbose=TRUE)
  sink()
}

foreach (i=c(15), .combine=c) %dopar% {
  DoHapseg2(myarray = info$myarray[i], mysample = paste(i,info$sample[i], sep = '.'))
}
foreach (i = c(2, 3, 12:14,16, 20:27, 29, 30, 34, 36, 51, 53, 56), .combine=c) %dopar% {
  DoHapseg2(myarray = info$myarray[i], mysample = paste(i,info$sample[i], sep = '.'))
}
foreach (i = setdiff(1:56, c(2, 3, 12:16, 20:27, 29, 30, 34, 36, 51, 53, 56)), .combine=c) %dopar% {
  DoHapseg2(myarray = info$myarray[i], mysample = paste(i,info$sample[i], sep = '.'))
}

DoHapseg3 <- function(myarray, mysample) {
  registerDoSEQ()
  plate.name <- mysample
  phased.bgl.dir <- "/usr/local/bioinfsoftware/R/current/lib64/R/library/ABSOLUTE/etc/phasedBGL/hg19"
  d = file.path(getwd(), 'ABSOLUTE', 'HAPSEG', mysample)
  snp.fn <- file.path(d, 'snpfileNoCalibration.RData')
  seg.fn = file.path(d, 'Segdat.txt')
  results.dir <- file.path(getwd(), 'ABSOLUTE', 'HAPSEG3', mysample)
  if (!file.exists(results.dir)) {
    dir.create(results.dir, recursive=TRUE)
  }
  print(paste("Starting array", mysample, "at", results.dir))
  log.dir <- file.path(getwd(), 'ABSOLUTE', 'HAPSEG3', 'log')
  if (!file.exists(log.dir)) {
    dir.create(log.dir, recursive=TRUE)
  }
  sink(file=file.path(log.dir, paste(mysample, ".hapseg.out.txt", sep="")))
  myRunHapSeg(plate.name, myarray, seg.fn, snp.fn, genome.build = 'hg19', results.dir, platform = "SNP_6.0", 
              use.pop='CEPH', 
              impute.gt=TRUE, plot.segfit=TRUE, merge.small=FALSE, merge.close=TRUE, min.seg.size=2, 
              normal=FALSE, 
              out.p=0.001, seg.merge.thresh=1e-80, use.normal = FALSE, adj.atten = FALSE, 
              phased.bgl.dir, calibrate.data=FALSE, 
              drop.x=FALSE, verbose=TRUE)
  sink()
}
foreach (i=c(15), .combine=c) %dopar% {
  DoHapseg3(myarray = info$myarray[i], mysample = paste(i,info$sample[i], sep = '.'))
}

DoHapseg4 <- function(myarray, mysample) {
  registerDoSEQ()
  plate.name <- mysample
  phased.bgl.dir <- "/usr/local/bioinfsoftware/R/current/lib64/R/library/ABSOLUTE/etc/phasedBGL/hg19"
  d = file.path(getwd(), 'ABSOLUTE', 'HAPSEG', mysample)
  snp.fn <- file.path(d, 'snpfile1Calibration.RData')
  seg.fn = file.path(d, 'Segdat.txt')
  results.dir <- file.path(getwd(), 'ABSOLUTE', 'HAPSEG4', mysample)
  if (!file.exists(results.dir)) {
    dir.create(results.dir, recursive=TRUE)
  }
  print(paste("Starting array", mysample, "at", results.dir))
  log.dir <- file.path(getwd(), 'ABSOLUTE', 'HAPSEG4', 'log')
  if (!file.exists(log.dir)) {
    dir.create(log.dir, recursive=TRUE)
  }
  sink(file=file.path(log.dir, paste(mysample, ".hapseg.out.txt", sep="")))
  myRunHapSeg(plate.name, myarray, seg.fn, snp.fn, genome.build = 'hg19', results.dir, platform = "SNP_6.0", 
              use.pop='CEPH', 
              impute.gt=TRUE, plot.segfit=TRUE, merge.small=FALSE, merge.close=TRUE, min.seg.size=2, 
              normal=FALSE, 
              out.p=0.001, seg.merge.thresh=1e-80, use.normal = FALSE, adj.atten = FALSE, 
              phased.bgl.dir, calibrate.data=FALSE, 
              drop.x=FALSE, verbose=TRUE)
  sink()
}
foreach (i=c(15), .combine=c) %dopar% {
  DoHapseg4(myarray = info$myarray[i], mysample = paste(i,info$sample[i], sep = '.'))
}

DoHapsegj <- function(myarray, mysample, j) {
  registerDoSEQ()
  plate.name <- mysample
  phased.bgl.dir <- "/usr/local/bioinfsoftware/R/current/lib64/R/library/ABSOLUTE/etc/phasedBGL/hg19"
  d = file.path(getwd(), 'ABSOLUTE', 'HAPSEG', mysample)
  snp.fn <- file.path(d, paste('snpfileCalibration', j,'.RData', sep = ''))
  seg.fn = file.path(d, 'Segdat.txt')
  results.dir <- file.path(getwd(), 'ABSOLUTE', paste('HAPSEG', j, sep = ''), mysample)
  if (!file.exists(results.dir)) {
    dir.create(results.dir, recursive=TRUE)
  }
  print(paste("Starting array", mysample, "at", results.dir))
  log.dir <- file.path(getwd(), 'ABSOLUTE', paste('HAPSEG', j, sep = ''), 'log')
  if (!file.exists(log.dir)) {
    dir.create(log.dir, recursive=TRUE)
  }
  sink(file=file.path(log.dir, paste(mysample, ".hapseg.out.txt", sep="")))
  myRunHapSeg(plate.name, myarray, seg.fn, snp.fn, genome.build = 'hg19', results.dir, platform = "SNP_6.0", 
              use.pop='CEPH', 
              impute.gt=TRUE, plot.segfit=TRUE, merge.small=FALSE, merge.close=TRUE, min.seg.size=2, 
              normal=FALSE, 
              out.p=0.001, seg.merge.thresh=1e-80, use.normal = FALSE, adj.atten = FALSE, 
              phased.bgl.dir, calibrate.data=FALSE, 
              drop.x=FALSE, verbose=TRUE)
  sink()
}
foreach (j=5:8, .combine=c) %dopar% {
  i = 15
  DoHapsegj(myarray = info$myarray[i], mysample = paste(i,info$sample[i], sep = '.'),
            j = j)
}

################################################################
################################################################
#
# Format exome mutation file for ABSOLUTE
#
################################################################
################################################################
# Create sample_list for Pearl script TSV2MAF_pearl_script.rtf
# and reformat mutations file to MAF (Broad's format)
setwd(paste(getwd(), '/RESPONSIFY', sep=''))
#setwd('/Users/lonnstedt/Documents/RESPONSIFY')
prefix.in = paste(getwd(), "/rawData/Exome/", sep='')
prefix.raw = paste(getwd(), "/datasets/", sep='')

###Sample_list
mut = read.delim(paste(prefix.in, 'Loi_130516_24samples/ALL_130516_Filter_CBS_MinDepth_20.txt', 
                       sep = ''))
samplesT = unique(mut$SAMPLE)
samples = sort(unique(substr(as.character(samplesT), 1, 4)))
sampT = paste(samples, 'T', sep = '')
sampN = paste(samples, 'N', sep = '')
out = rbind(data.frame(sample = sampT, type = 'Tumour', pairsamp = sampN, bam = 'dummy'), 
            data.frame(sample = sampN, type = 'Normal', pairsamp = sampT, bam = 'dummy'))
write.table(out, file = paste(prefix.in, 'sample_list.txt', sep = ''), row.names = F,col.names = F,
            quote = F, sep = '\t')

###Pearl calls
###Run from Exome library
###This never worked with our system, but Jason Li run it for us.
perl parse_tsv_to_maf_format.pl sample_list.txt HUGO_EntrezIDs.txt Loi_130516_24samples/ALL_130516_Filter_CBS_MinDepth_20.txt ALL_130516_Filter_CBS_MinDepth_20.maf


### Find additional variables to merge with the maf file
### AND split maf file into one file for each sample
mut = read.delim(paste(prefix.in, 'Loi_130516_24samples/ALL_130516_Filter_CBS_MinDepth_20.txt', 
                       sep = ''))
dat = data.frame(Tumor_Sample_Barcode = mut$SAMPLE,
                 Chromosome = mut$X.CHROM, Start_position = mut$POS, 
                 t_alt_count = mut$X4688T_Variant_reads,
                 t_ref_count = mut$X4688T_Read_depth - mut$X4688T_Variant_reads)
dat = unique(dat)
mafbase = 'ALL_130516_Filter_CBS_MinDepth_20'
maf = read.delim(paste(prefix.in, mafbase,'.maf', sep = ''), 
                 row.names = NULL, stringsAsFactors = FALSE, 
                 check.names = FALSE, na.strings = c("NA", "---"), 
                 blank.lines.skip = TRUE, comment.char = "#")
names(maf)[names(maf) == 'Start_Position'] = 'Start_position'
names(maf)[names(maf) == 'End_Position'] = 'End_position'
nrow(maf)
maf = merge(maf, dat, all = T, sort = F)
nrow(maf)

info = read.delim(paste(prefix.ann, '../../Infokey.txt', sep = ''))
info$sample = substr(as.character(info$Sample.ID), 1, 4)
for (i in 1:nrow(info)){
  samp = info$sample[i]
  out = subset(maf, Tumor_Sample_Barcode == paste(samp, 'T', sep = ''))
  if (nrow(out)>0) write.table(out, file = paste(prefix.in,'MAF/', 
              mafbase, '_', i,'.', samp, '.maf', sep = ''), row.names = F, quote = F,
                               sep = '\t')
}


#maf = read.delim(maf.fn, 
#                 row.names = NULL, stringsAsFactors = FALSE, 
#                 check.names = FALSE, na.strings = c("NA", "---"), 
#                 blank.lines.skip = TRUE, comment.char = "#")





################################################################
################################################################
#
# Run ABSOLUTE without scaling, just for comparison
#
################################################################
################################################################
library(ABSOLUTE)
library(foreach)
library(doMC)
info = read.delim(paste(prefix.ann, '../../Infokey.txt', sep = ''))
info$sample = substr(info$Sample.ID, 1, 4)
info$myarray = info$arrname
info$sample.name = paste(rownames(info),info$sample, sep = '.')

DoAbsolute <- function(myarray, mysample, i) {
  registerDoSEQ()
  sample.name <- mysample
  maf.fn = paste(getwd(), 
     '/AROMA/rawData/responsify/Exome/MAF/ALL_130516_Filter_CBS_MinDepth_20', '_', 
                 mysample, '.maf', sep = '')
  seg.dat.fn <- paste(getwd(), '/ABSOLUTE/HAPSEG2/', mysample,'/',mysample,'_',myarray,
                      "_segdat.RData", sep="")
  results.dir <- paste(getwd(), '/ABSOLUTE/RunAbsolute2', sep = '')
  print(paste("Starting myarray", mysample, "at", results.dir))
  log.dir <- paste(getwd(), '/ABSOLUTE/RunAbsolute2/abs_logs', sep = '') 
  if (!file.exists(log.dir)) {
    dir.create(log.dir, recursive=TRUE)
  }
  if (!file.exists(results.dir)) {
    dir.create(results.dir, recursive=TRUE)
  }
  sink(file=file.path(log.dir, paste(mysample, ".abs.out.txt", sep="")))
  myRunAbsolute(seg.dat.fn, sigma.p = 0, max.sigma.h = 0.02, min.ploidy = 0.95, max.ploidy = 6, 
              primary.disease = "Breast Cancer", 
              platform = 'SNP_6.0', sample.name, results.dir, max.as.seg.count = 1500, 
              max.non.clonal = 0, 
              max.neg.genome = 0, copy_num_type = 'allelic', verbose=TRUE, maf.fn = maf.fn, 
              min.mut.af = 0.03, output.fn.base = mysample)
  sink()
}


foreach (i = c(2, 3, 12:14,16, 20:27, 29, 30, 34, 36, 51, 53, 56), .combine=c) %dopar% {
  DoAbsolute(myarray = info$myarray[i], mysample = paste(i,info$sample[i], sep = '.'))
}
foreach (i = c(24:27, 29, 30, 34, 36, 51, 53, 56), .combine=c) %dopar% {
  DoAbsolute(myarray = info$myarray[i], mysample = paste(i,info$sample[i], sep = '.'))
}

d = paste(getwd(), '/ABSOLUTE/RunAbsolute2', sep = '')
obj.name = 'ALL'
CreateReviewObject(obj.name=obj.name, absolute.files=paste(d, '/', 
                    info$sample.name[c(2, 3, 12:16, 20:27, 29, 30, 34, 36, 51, 53, 56)], 
                                                           '.ABSOLUTE.RData', sep = ''), 
                   indv.results.dir = file.path(d, obj.name), copy_num_type = "allelic", 
                   plot.modes = TRUE, verbose=TRUE)



DoAbsolute3 <- function(myarray, mysample, i) {
  registerDoSEQ()
  sample.name <- mysample
  maf.fn = paste(getwd(), 
                 '/AROMA/rawData/responsify/Exome/MAF/ALL_130516_Filter_CBS_MinDepth_20', '_', 
                 mysample, '.maf', sep = '')
  seg.dat.fn <- paste(getwd(), '/ABSOLUTE/HAPSEG3/', mysample,'/',mysample,'_',myarray,
                      "_segdat.RData", sep="")
  results.dir <- paste(getwd(), '/ABSOLUTE/RunAbsolute3', sep = '')
  print(paste("Starting myarray", mysample, "at", results.dir))
  log.dir <- paste(getwd(), '/ABSOLUTE/RunAbsolute3/abs_logs', sep = '') 
  if (!file.exists(log.dir)) {
    dir.create(log.dir, recursive=TRUE)
  }
  if (!file.exists(results.dir)) {
    dir.create(results.dir, recursive=TRUE)
  }
  sink(file=file.path(log.dir, paste(mysample, ".abs.out.txt", sep="")))
  myRunAbsolute(seg.dat.fn, sigma.p = 0, max.sigma.h = 0.02, min.ploidy = 0.95, max.ploidy = 6, 
                primary.disease = "Breast Cancer", 
                platform = 'SNP_6.0', sample.name, results.dir, max.as.seg.count = 1500, 
                max.non.clonal = 0, 
                max.neg.genome = 0, copy_num_type = 'allelic', verbose=TRUE, maf.fn = maf.fn, 
                min.mut.af = 0.03, output.fn.base = mysample)
  sink()
}
foreach (i = c(15), .combine=c) %dopar% {
  DoAbsolute3(myarray = info$myarray[i], mysample = paste(i,info$sample[i], sep = '.'))
}

DoAbsolute4 <- function(myarray, mysample, i) {
  registerDoSEQ()
  sample.name <- mysample
  maf.fn = paste(getwd(), 
                 '/AROMA/rawData/responsify/Exome/MAF/ALL_130516_Filter_CBS_MinDepth_20', '_', 
                 mysample, '.maf', sep = '')
  seg.dat.fn <- paste(getwd(), '/ABSOLUTE/HAPSEG4/', mysample,'/',mysample,'_',myarray,
                      "_segdat.RData", sep="")
  results.dir <- paste(getwd(), '/ABSOLUTE/RunAbsolute4', sep = '')
  print(paste("Starting myarray", mysample, "at", results.dir))
  log.dir <- paste(getwd(), '/ABSOLUTE/RunAbsolute4/abs_logs', sep = '') 
  if (!file.exists(log.dir)) {
    dir.create(log.dir, recursive=TRUE)
  }
  if (!file.exists(results.dir)) {
    dir.create(results.dir, recursive=TRUE)
  }
  sink(file=file.path(log.dir, paste(mysample, ".abs.out.txt", sep="")))
  myRunAbsolute(seg.dat.fn, sigma.p = 0, max.sigma.h = 0.02, min.ploidy = 0.95, max.ploidy = 6, 
                primary.disease = "Breast Cancer", 
                platform = 'SNP_6.0', sample.name, results.dir, max.as.seg.count = 1500, 
                max.non.clonal = 0, 
                max.neg.genome = 0, copy_num_type = 'allelic', verbose=TRUE, maf.fn = maf.fn, 
                min.mut.af = 0.03, output.fn.base = mysample)
  sink()
}
foreach (i = c(15), .combine=c) %dopar% {
  DoAbsolute4(myarray = info$myarray[i], mysample = paste(i,info$sample[i], sep = '.'))
}

DoAbsolutej <- function(myarray, mysample, j) {
  registerDoSEQ()
  sample.name <- mysample
  maf.fn = paste(getwd(), 
                 '/AROMA/rawData/responsify/Exome/MAF/ALL_130516_Filter_CBS_MinDepth_20', '_', 
                 mysample, '.maf', sep = '')
  seg.dat.fn <- paste(getwd(), '/ABSOLUTE/HAPSEG',j,'/', mysample,'/',mysample,'_',
                      myarray,
                      "_segdat.RData", sep="")
  results.dir <- paste(getwd(), '/ABSOLUTE/RunAbsolute', j, sep = '')
  print(paste("Starting myarray", mysample, "at", results.dir))
  log.dir <- paste(getwd(), '/ABSOLUTE/RunAbsolute', j, '/abs_logs', sep = '') 
  if (!file.exists(log.dir)) {
    dir.create(log.dir, recursive=TRUE)
  }
  if (!file.exists(results.dir)) {
    dir.create(results.dir, recursive=TRUE)
  }
  sink(file=file.path(log.dir, paste(mysample, ".abs.out.txt", sep="")))
  myRunAbsolute(seg.dat.fn, sigma.p = 0, max.sigma.h = 0.02, min.ploidy = 0.95, max.ploidy = 6, 
                primary.disease = "Breast Cancer", 
                platform = 'SNP_6.0', sample.name, results.dir, max.as.seg.count = 1500, 
                max.non.clonal = 0, 
                max.neg.genome = 0, copy_num_type = 'allelic', verbose=TRUE, maf.fn = maf.fn, 
                min.mut.af = 0.03, output.fn.base = mysample)
  sink()
}
foreach (j=5:8, .combine=c) %dopar% {
  i = 15
  DoAbsolutej(myarray = info$myarray[i], mysample = paste(i,info$sample[i], sep = '.'),
            j = j)
}


################################################################
################################################################
#
# Check HAPSEG solution
#
################################################################
################################################################
info = read.delim(paste(prefix.ann, '../../Infokey.txt', sep = ''), stringsAsFactors = F)
info$sample.name = paste(info$ind, info$sample, sep = '.')
################################################################
######Find clusters of segments which we like to check

for (i in 1:56){
  
  i = 56
  ta = read.delim(file.path(getwd(), 'TAPS', info$sample.name[i], 
                            paste(info$sample.name[i], 'segmentCN.txt', sep = '_')))
  plot(imba~log2, data = ta, xlim = c(-1,2), main = info$sample.name[i])
  
  d0.imba = diff(range(ta$imba, na.rm = T))/200
  d0.log = diff(range(ta$log2, na.rm = T))/1000
  out = NULL
  
  n = 1 ###MANUALLY INSERT THE NUMBER OF CLUSTERS YOU LIKE HERE:
  pos = locator(n)
  
  ###FOR EACH CLUSTER (CHANGE c BETWEEN 1, 2, 3, ..., n)
  c = 1
  plot(imba~log2, data = ta, xlim = c(-1,2), main = info$sample.name[i])
  di = 0
  dl = 0
  
  #########RUN THESE LINES UNTIL JUST TOO MANY SEGMENTS HIGHLIGHTED
  if (di>0) last.tmp = subset(ta, abs(ta$imba - pos$y[c])<di & abs(ta$log2 - pos$x[c])<dl)
  di = di + d0.imba
  dl = dl + d0.log
  tmp = subset(ta, abs(ta$imba - pos$y[c])<di & abs(ta$log2 - pos$x[c])<dl)
  points(imba~log2, data = tmp, col = c+1, pch = 16)
  
  #########THEN SAVE THE PREVIOUS SEGMENTS (RUN THIS ONCE)
  out = rbind(out, cbind(last.tmp, cluster = c))    
  plot(imba~log2, data = ta, xlim = c(-1,2), main = info$sample.name[i], bg = 'white')
  points(imba~log2, data = last.tmp, col = c+1, pch = 16)
  
  
  
  plot(imba~log2, data = ta, xlim = c(-1,2), main = info$sample.name[i])
  points(imba~log2, data = out, col = out$cluster+1, pch = 16)
  d = file.path(getwd(), 'ABSOLUTE', 'Choose_solution')
  dev.print(png, file=file.path(d, paste(info$sample.name[i], '_Chosen_subclones.png', sep = '')), 
            width=640, height=640)
  d = file.path(getwd(), 'ABSOLUTE', 'HAPSEG', info$sample.name[i])
  dir.create(d, recursive = T, showWarnings = F)
  write.table(out[,c('Chromosome','Start','End','cluster')], file.path(d, 'Subsegs.txt'), 
              row.names = F, quote = F, sep = '\t')
}


################################################################
######Plot Allele1 vs Allele2 from TAPS coloured by the chosen clusters
library(PSCBS)

i = 15
d = file.path(getwd(), 'AROMA', 'reports', 'SNP_segmentation', 'CBSfits')
da = loadObject(paste(d, '/fit', i, '.RData', sep = ''))
seg = da$output
plot(dhMean~log(tcnMean), data = seg,  main = info$sample.name[i], xlim = c(-.2, 2))
plot(c1Mean~c2Mean, data = seg, main = info$sample.name[i])
grid()
mat = seg
names(mat)[names(mat) == 'c1Mean'] = 'A1'
names(mat)[names(mat) == 'c2Mean'] = 'A2'
names(mat)[names(mat) == 'tcnStart'] = 'Start.bp'
names(mat)[names(mat) == 'tcnEnd'] = 'End.bp'
names(mat)[names(mat) == 'chromosome'] = 'Chromosome'

mat$length = mat$End.bp - mat$Start.bp
mat = subset(mat, !is.na(length) & Chromosome<23)

###Read selected segments
d = file.path(getwd(), 'ABSOLUTE', 'HAPSEG', info$sample.name[i])
tmp = getKnownClustersAndOverlaps(mat = mat, d = d)
cl = tmp$cl
clust = tmp$clust
w = tmp$w  


#A1 vs A2
plot1v2(a1 = mat$A1, a2=mat$A2, mat = mat, clust = clust, cl = cl, w = w)

d = file.path(getwd(), 'ABSOLUTE','ScaledAbsolute')
dev.print(png, file=file.path(d, paste(info$sample.name[i], '.............png', sep = '')), 
          width=640, height=640)
#This plot, too has the problem with CN 2 clusters not on a vertical line.

################################################################
######Plot HAPSEG HSCRs coloured by the chosen clusters

d = file.path(getwd(), 'ABSOLUTE', 'Choose_clusters', 'Cluster_overview')
pdf(paste(d, 'Cluster_overview2', '.pdf', sep=''), 
    paper = 'a4', height = 12, width = 8)  

for (i in 1:56){
  
#  d = file.path(getwd(), 'ABSOLUTE', 'Choose_clusters', 'Cluster_overview')
#  pdf(paste(d, '/', info$sample.name[i],'_Cluster_overview2', '.pdf', sep=''), 
#      paper = 'a4', height = 12, width = 8)  
  
  
  ### Load data for this sample
  d = file.path(getwd(), 'ABSOLUTE', 'HAPSEG2', info$sample.name[i])
  f = load(file.path(d, paste(info$sample.name[i], '_', info$arrname[i], '_segdat.RData', sep = '')))
  ab = as.data.frame(seg.dat$allele.segs)
  ab$Chromosome = factor(ab$Chromosome, levels = as.character(1:22))
  names(ab)[6:7] = c('A1','A2')
  
  ###Read selected segments
  d = file.path(getwd(), 'ABSOLUTE', 'HAPSEG', info$sample.name[i])
  tmp = getKnownClustersAndOverlaps(mat = ab, d = d)
  cl = tmp$cl
  clust = tmp$clust
  w = tmp$w  
  
  ### Plots of HAPSEG output
  ###
  par(mfrow = c(3,2))
  #imba by logratio from hscrs
  plot1v2(a1 = 2*abs(ab$A1/(ab$A1+ab$A2) - 0.5), a2=log(ab$A1+ab$A2), mat = ab, clust = clust, cl = cl, w = w,
          xlab = 'log(Total Copy Number)', ylab = 'Alleleic imbalance 2*|BAF - 0.5|')
  
  #A1 vs A2
  plot1v2(a1 = ab$A1, a2=ab$A2, mat = ab, clust = clust, cl = cl, w = w)
  
  ### Plots of histograms
  myab = histCol(mat = ab, cl = cl, clust = clust, xlab = 'Pooled haplytype copy values')
  myab = histCol(toplot = 'A1', cl = cl, clust=clust, myab = myab, xlab = 'Haplytype 1 copy values')
  frame()
  myab = histCol(toplot = 'A2', cl = cl, clust=clust, myab = myab, xlab = 'Haplytype 2 copy values')
  
   
  title(main = info$sample.name[i], outer = T, line = -1)
#  dev.off()
}
dev.off()



################################################################
################################################################
#
# Sample 15
#
################################################################
################################################################
info = read.delim(paste(prefix.ann, '../../Infokey.txt', sep = ''), stringsAsFactors = F)
info$sample.name = paste(info$ind, info$sample, sep = '.')
source(file.path(getwd(), 'programs', 'make_HAPSEG_seg_obj.R'))
library(Rsolnp)

i = 15
d = file.path(getwd(), 'ABSOLUTE', 'HAPSEG3')
d = paste(d,'/', info$sample.name[i], sep = '')
#load(file.path(d, paste(info$sample.name[i], '_', info$arrname[i], '_segdat.RData', sep = '')))
seg.dat.fn = file.path(d, paste(info$sample.name[i], '_', info$arrname[i], '_segdat.RData', sep = ''))
seg.dat = HAPSEGMakeSegObj(seg.dat.fn, gender = 'Female', filter_segs=FALSE, min_probes=10, max_sd=100, 
                           verbose=TRUE)  
as.seg.dat = seg.dat$as.seg.dat
nr = as.integer(nrow(as.seg.dat))
mat = as.seg.dat[1:(nr/2),]
names(mat)[6] = 'A1'
mat$A2 = as.seg.dat$copy_num[(nr/2+1):(nr)]
mat$W = mat$W*2
  
###Read selected segments
tmp = getKnownClustersAndOverlaps(mat = mat, d = file.path(getwd(),'ABSOLUTE', 'HAPSEG', info$sample.name[i]))
cl = tmp$cl
clust = tmp$clust
w = tmp$w

#Define data for optimization
y = mat$A1
x = mat$A2
xclust = rep(NA, nrow(mat))
yclust = rep(NA, nrow(mat))
cv = read.delim(file.path(getwd(), 'ABSOLUTE', 'Choose_clusters', 'ClusterCNs.txt'))
cv = cv[cv$sample == i,]

ix = which(clust == 'cluster hom')
if (length(ix)>0) xclust[ix] = cv$xhom
if (length(ix)>0) yclust[ix] = cv$yhom

ix = which(clust == 'cluster het')
if (length(ix)>0) xclust[ix] = cv$xhet
if (length(ix)>0) yclust[ix] = cv$yhet

ix = which(clust == 'cluster 1')
if (length(ix)>0) xclust[ix] = cv$x1
if (length(ix)>0) yclust[ix] = cv$y1

ix = which(clust == 'cluster 2')
if (length(ix)>0) xclust[ix] = cv$x2
if (length(ix)>0) yclust[ix] = cv$y2

ix = which(clust == 'cluster 3')
if (length(ix)>0) xclust[ix] = cv$x3
if (length(ix)>0) yclust[ix] = cv$y3

ix = which(clust == 'cluster 4')
if (length(ix)>0) xclust[ix] = cv$x4
if (length(ix)>0) yclust[ix] = cv$y4


#Save the cluster data to use
dd = data.frame(x=x, y=y, xclust = xclust, yclust=yclust, xw=w/sum(w), yw = w/sum(w),
                myclust = clust)
dd = subset(dd, !(is.na(xclust) & is.na(yclust)))


##########################################################
#Optimize one issue at the time


### Scale x
###
dx = subset(dd, !is.na(xclust))
pars = solnp(pars = c(1,0), fun = minx, myclust = dx$myclust,
             X=dx$x, Y=dx$y, xclust = dx$xclust, xw = dx$xw, 
             eqfun = eqfunx, eqB = mean(dx$x), LB = c(0, -1000))$pars
a = pars[1]
b = pars[2]
mat$A2new = a*mat$A2 + b*mat$A1
dd$x = a*dd$x + b*dd$y
param = c(a, b)

### Scale y
###
dy = subset(dd, !is.na(yclust))
pars = solnp(pars = c(0,1), fun = miny, myclust = dy$myclust,
             X=dy$x, Y=dy$y, yclust = dy$yclust, yw = dy$yw, 
             eqfun = eqfuny, eqB = mean(dy$y), LB = c(-1000, 0))$pars
c = pars[1]
d = pars[2]
mat$A1new = c*mat$A2new + d*mat$A1
dd$y = c*dd$x + d*dd$y
param = c(param, c, d)

### Scale  so that x and y fit together
###
e.cr1 = sum(c(mat$A1, mat$A2)* rep(mat$W, 2)/2)
pars = solnp(pars = c(1, 1), fun = minxy, myclust = dd$myclust, mat = mat,
             eqfun = eqfuncxy, eqB = 1,
             X=dd$x, Y=dd$y, xclust = dd$xclust, xw = dd$xw, yclust = dd$yclust, yw = dd$yw)$pars
B = pars[1]
D = pars[2]
mat$A2last =  B*mat$A2new
mat$A1last =  D*mat$A1new
param = c(param, B, D)
names(param) = c('a','b','c','d', 'B','D')
dd$x = B*dd$x
dd$y = D*dd$y



###########################################################
# Figures of outcome
plot1v2(a1 = mat$A1, a2 = mat$A2, mat = mat, cl = cl, clust = clust, w = w, xlab = 'Relative CN homologue 2',
        ylab = 'Relative CN homologue 1')
d = file.path(getwd(), 'ABSOLUTE','ScaledAbsolute3')
dev.print(png, file=file.path(d, paste(info$sample.name[i], '_A1vA2orig.png', sep = '')), 
          width=640, height=640)
plot1v2(a1 = mat$A1, a2 = mat$A2new, mat = mat, cl = cl, clust = clust, w = w, xlab = 'Relative CN homologue 2',
        ylab = 'Relative CN homologue 1')
dev.print(png, file=file.path(d, paste(info$sample.name[i], '_A1vA2scaledX.png', sep = '')), 
          width=640, height=640)
plot1v2(a1 = mat$A1new, a2 = mat$A2new, mat = mat, cl = cl, clust = clust, w = w, xlab = 'Relative CN homologue 2',
        ylab = 'Relative CN homologue 1')
dev.print(png, file=file.path(d, paste(info$sample.name[i], '_A1vA2scaledXY.png', sep = '')), 
          width=640, height=640)
plot1v2(a1 = mat$A1last, a2 = mat$A2last, mat = mat, cl = cl, clust = clust, w = w, xlab = 'Relative CN homologue 2',
        ylab = 'Relative CN homologue 1')
dev.print(png, file=file.path(d, paste(info$sample.name[i], '_A1vA2scaledAll.png', sep = '')), 
          width=640, height=640)



#Histograms
myab = histCol(mat = mat, cl = cl, xlab = 'Relative CN (homologues pooled)')
#d = file.path(getwd(), 'ABSOLUTE','ScaledAbsolute8')
dev.print(png, file=file.path(d, paste(info$sample.name[i], '_Historig.png', sep = '')), 
          width=640, height=640)

myab = histCol(toplot = c('A1last','A2last'), myab = myab, cl = cl, xlab = 'Relative CN (homologues pooled)')
dev.print(png, file=file.path(d, paste(info$sample.name[i], '_HistscaledAll.png', sep = '')), 
          width=640, height=640)

###########################################################


  

#Reset attenuation constant
seg.dat$error.model$AT = 0

#New version or hscr values in as.seg.dat
as.seg.dat$which = 1
as.seg.dat$which[(nr/2+1):nr] = 2
upper = as.seg.dat[as.seg.dat$which == 1,]
tmp = mat[,c('Chromosome','Start.bp','End.bp', 'A1last')]
upper = merge(upper, tmp, all = T)
upper$copy_num = upper$A1last
upper = upper[order(upper$seg.ix),1:9]
lower = as.seg.dat[as.seg.dat$which == 2,]
tmp = mat[,c('Chromosome','Start.bp','End.bp', 'A2last')]
lower = merge(lower, tmp, all = T)
lower$copy_num = lower$A2last
lower = lower[order(lower$seg.ix),1:9]
as.seg.dat = rbind(upper, lower)
seg.dat$as.seg.dat = as.seg.dat


#New version of allele.segs
allele.segs = seg.dat$allele.segs
allele.segs = as.data.frame(allele.segs)
tmp = mat[,c('Chromosome','Start.bp','End.bp', 'A1last', 'A2last')]
allele.segs = merge(allele.segs, tmp, sort = F)
allele.segs$A1.Seg.CN = allele.segs$A1last
allele.segs$A2.Seg.CN = allele.segs$A2last
allele.segs = as.matrix(allele.segs[,1:9])
seg.dat$allele.segs = allele.segs

#New version of seg.info
seg.info = seg.dat$seg.info
seg.info = as.data.frame(seg.info)
tmp = mat[,c('Chromosome','Start.bp','End.bp', 'A1last', 'A2last')]
seg.info = merge(seg.info, tmp, sort = F)
seg.info$cn = (seg.info$A1last + seg.info$A2last)/2
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

#Print the new object
d = file.path(getwd(), 'ABSOLUTE', 'HAPSEG8', info$sample.name[i])
save(seg.dat, file = file.path(d, paste(info$sample.name[i], '_', info$arrname[i], 
                                        '_scaled_segdat.RData', sep = '')))



################################################################
################################################################
#
# Run ABSOLUTE with scaled haplotype copy values
#
################################################################
################################################################
library(ABSOLUTE)
library(foreach)
library(doMC)
info = read.delim(paste(prefix.ann, '../../Infokey.txt', sep = ''))
info$sample = substr(info$Sample.ID, 1, 4)
info$myarray = info$arrname
info$sample.name = paste(rownames(info),info$sample, sep = '.')



DoAbsolute <- function(myarray, mysample, i) {
  registerDoSEQ()
  sample.name <- mysample
  maf.fn = paste(getwd(), 
                 '/AROMA/rawData/responsify/Exome/MAF/ALL_130516_Filter_CBS_MinDepth_20', '_', 
                 mysample, '.maf', sep = '')
  seg.dat.fn <- paste(getwd(), '/ABSOLUTE/HAPSEG2/', mysample,'/',mysample,'_',myarray,
                      "_scaled_segdat.RData", sep="")
  results.dir <- paste(getwd(), '/ABSOLUTE/ScaledAbsolute2', sep = '')
  print(paste("Starting myarray", mysample, "at", results.dir))
  log.dir <- paste(getwd(), '/ABSOLUTE/ScaledAbsolute2/abs_logs', sep = '') 
  if (!file.exists(log.dir)) {
    dir.create(log.dir, recursive=TRUE)
  }
  if (!file.exists(results.dir)) {
    dir.create(results.dir, recursive=TRUE)
  }
  sink(file=file.path(log.dir, paste(mysample, ".abs.out.txt", sep="")))
  myRunAbsolute(seg.dat.fn, sigma.p = 0, max.sigma.h = 0.01, min.ploidy = 2, max.ploidy = 4, 
              primary.disease = "Breast Cancer", 
              platform = 'SNP_6.0', sample.name, results.dir, max.as.seg.count = 1500, 
              max.non.clonal = 0, 
              max.neg.genome = 0, copy_num_type = 'allelic', verbose=TRUE, maf.fn = maf.fn, 
              min.mut.af = 0.03, output.fn.base = mysample)
  sink()
}


i = 15
DoAbsolute(myarray = info$myarray[i], mysample = info$sample.name[i], i)

d = paste(getwd(), '/ABSOLUTE/ScaledAbsolute2', sep = '')
obj.name = '15'
CreateReviewObject(obj.name=obj.name, absolute.files=paste(d, '/', 
                     info$sample.name[15], '.ABSOLUTE.RData', sep = ''), 
                   indv.results.dir = file.path(d, obj.name), copy_num_type = "allelic", 
                   plot.modes = TRUE, verbose=TRUE)



DoAbsolute3 <- function(myarray, mysample, i) {
  registerDoSEQ()
  sample.name <- mysample
  maf.fn = paste(getwd(), 
                 '/AROMA/rawData/responsify/Exome/MAF/ALL_130516_Filter_CBS_MinDepth_20', '_', 
                 mysample, '.maf', sep = '')
  seg.dat.fn <- paste(getwd(), '/ABSOLUTE/HAPSEG3/', mysample,'/',mysample,'_',myarray,
                      "_scaled_segdat.RData", sep="")
  results.dir <- paste(getwd(), '/ABSOLUTE/ScaledAbsolute3', sep = '')
  print(paste("Starting myarray", mysample, "at", results.dir))
  log.dir <- paste(getwd(), '/ABSOLUTE/ScaledAbsolute3/abs_logs', sep = '') 
  if (!file.exists(log.dir)) {
    dir.create(log.dir, recursive=TRUE)
  }
  if (!file.exists(results.dir)) {
    dir.create(results.dir, recursive=TRUE)
  }
  sink(file=file.path(log.dir, paste(mysample, ".abs.out.txt", sep="")))
  myRunAbsolute2(seg.dat.fn, sigma.p = 0, max.sigma.h = 0.01, min.ploidy = 2, max.ploidy = 12, 
                primary.disease = "Breast Cancer", 
                platform = 'SNP_6.0', sample.name, results.dir, max.as.seg.count = 1500, 
                max.non.clonal = 0, 
                max.neg.genome = 0, copy_num_type = 'allelic', verbose=TRUE, maf.fn = maf.fn, 
                min.mut.af = 0.03, output.fn.base = mysample)
  sink()
}


i = 15
DoAbsolute3(myarray = info$myarray[i], mysample = info$sample.name[i], i)



DoAbsolutej <- function(myarray, mysample, j) {
  registerDoSEQ()
  sample.name <- mysample
  maf.fn = paste(getwd(), 
                 '/AROMA/rawData/responsify/Exome/MAF/ALL_130516_Filter_CBS_MinDepth_20', '_', 
                 mysample, '.maf', sep = '')
  seg.dat.fn <- paste(getwd(), '/ABSOLUTE/HAPSEG',j,'/', mysample,'/',mysample,'_',
                      myarray,
                      "_scaled_segdat.RData", sep="")
  results.dir <- paste(getwd(), '/ABSOLUTE/ScaledAbsolute', j, sep = '')
  print(paste("Starting myarray", mysample, "at", results.dir))
  log.dir <- paste(getwd(), '/ABSOLUTE/ScaledAbsolute', j, '/abs_logs', sep = '') 
  if (!file.exists(log.dir)) {
    dir.create(log.dir, recursive=TRUE)
  }
  if (!file.exists(results.dir)) {
    dir.create(results.dir, recursive=TRUE)
  }
  sink(file=file.path(log.dir, paste(mysample, ".abs.out.txt", sep="")))
  myRunAbsolute2(seg.dat.fn, sigma.p = 0, max.sigma.h = 0.01, min.ploidy = 2.5, max.ploidy = 12, 
                primary.disease = "Breast Cancer", 
                platform = 'SNP_6.0', sample.name, results.dir, max.as.seg.count = 1500, 
                max.non.clonal = 0, 
                max.neg.genome = 0, copy_num_type = 'allelic', verbose=TRUE, maf.fn = maf.fn, 
                min.mut.af = 0.03, output.fn.base = mysample)
  sink()
}
foreach (j=8, .combine=c) %dopar% {
  i = 15
  DoAbsolutej(myarray = info$myarray[i], mysample = paste(i,info$sample[i], sep = '.'),
              j = j)
}

d = paste(getwd(), '/ABSOLUTE/ScaledAbsolute3', sep = '')
obj.name = '15'
CreateReviewObject(obj.name=obj.name, absolute.files=paste(d, '/', 
                                                           info$sample.name[15], '.ABSOLUTE.RData', sep = ''), 
                   indv.results.dir = file.path(d, obj.name), copy_num_type = "allelic", 
                   plot.modes = TRUE, verbose=TRUE)

obj.name = '15'
d = paste(getwd(), '/ABSOLUTE/ScaledAbsolute3', sep = '')
calls.path = file.path(d, obj.name, paste(obj.name, ".PP-calls_tab.txt", sep = ''))
modes.path = file.path(d, obj.name, paste(obj.name, ".PP-modes.data.RData", sep = ''))
output.path = file.path(d, obj.name)
ExtractReviewedResults(reviewed.pp.calls.fn = calls.path, analyst.id = "IL", 
                       modes.fn = modes.path, out.dir.base = output.path, 
                       obj.name = "AfterScaling", copy_num_type="allelic")

obj.name = '15_subset'
myCreateReviewObject(obj.name=obj.name, absolute.files=paste(d, '/', info$sample.name[c(15)], 
                                                             '.ABSOLUTE.RData', sep = ''), 
                     indv.results.dir = file.path(d, obj.name), copy_num_type = "allelic", 
                     plot.modes = TRUE, verbose=TRUE, myModes = c(1:8,12))
obj.name = '15_subset'
d = paste(getwd(), '/ABSOLUTE/ScaledAbsolute3', sep = '')
calls.path = file.path(d, obj.name, paste(obj.name, ".PP-calls_tab.txt", sep = ''))
modes.path = file.path(d, obj.name, paste(obj.name, ".PP-modes.data.RData", sep = ''))
output.path = file.path(d, obj.name)
ExtractReviewedResults(reviewed.pp.calls.fn = calls.path, analyst.id = "IL", 
                       modes.fn = modes.path, out.dir.base = output.path, 
                       obj.name = "AfterScaling", copy_num_type="allelic")


#######Junk
foreach (i = c(2, 3, 12:16, 20:27, 29, 30, 34, 36, 51, 53, 56), .combine=c) %dopar% {
  DoAbsolute(myarray = info$myarray[i], mysample = info$sample.name[i], i)
}

registerDoSEQ()
sample.name <- mysample
maf.fn = paste(getwd(), 
               '/AROMA/rawData/responsify/Exome/MAF/ALL_130516_Filter_CBS_MinDepth_20', '_', 
               mysample, '.maf', sep = '')
seg.dat.fn <- paste(prefix.out, mysample,'/',mysample,'_',myarray,
                    "_segdat.RData", sep="")
results.dir <- paste(getwd(), '/ABSOLUTE/RunAbsolute2', sep = '')
print(paste("Starting myarray", mysample, "at", results.dir))
log.dir <- paste(getwd(), '/ABSOLUTE/RunAbsolute2/abs_logs', sep = '') 
if (!file.exists(log.dir)) {
  dir.create(log.dir, recursive=TRUE)
}
if (!file.exists(results.dir)) {
  dir.create(results.dir, recursive=TRUE)
}
              sigma.p = 0
              max.sigma.h = 0.01
              min.ploidy = 0.95
              max.ploidy = 6
              primary.disease = "Breast Cancer"
               platform = 'SNP_6.0'
             max.as.seg.count = 1500
              max.non.clonal = 0
              max.neg.genome = 0
               copy_num_type = 'allelic'
              verbose=TRUE
               min.mut.af = 0.03
output.fn.base = mysample
copy_num_type = 'allelic'



