################################################################
################################################################
# Script: ProgR_5.3_Heterogeneity.R
# Author: Ingrid Lonnstedt
# Date:  23/05/2013
# R version: R 3.0
# Details: SNP analysis
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
  homs = read.delim(file.path(d, 'Homsegs.txt'), stringsAsFactors = F)
  hets = read.delim(file.path(d, 'Hetsegs.txt'), stringsAsFactors = F)
  data = read.delim(file.path(d, 'Total,fracB.txt'), stringsAsFactors = F)
  
  nc = nchar(homs$Chromosome, type = 'c')
  homs$Chromosome = substr(homs$Chromosome, 4, nc)
  homs$Chromosome[homs$Chromosome == 'X'] = '23'
  homs$Chromosome[homs$Chromosome == 'Y'] = '24'
  homs$Chromosome[homs$Chromosome == 'M'] = '25'
  homs$Chromosome = as.integer(homs$Chromosome)
  
  nc = nchar(hets$Chromosome, type = 'c')
  hets$Chromosome = substr(hets$Chromosome, 4, nc)
  hets$Chromosome[hets$Chromosome == 'X'] = '23'
  hets$Chromosome[hets$Chromosome == 'Y'] = '24'
  hets$Chromosome[hets$Chromosome == 'M'] = '25'
  hets$Chromosome = as.integer(hets$Chromosome)
  
  subd = subset(data, !is.na(chromosome) & !is.na(position) & !is.na(fracB) & !is.na(total) & !is.na(unitNames))
  subd$cn2 = 'no'
  subd$chromosome = factor(subd$chromosome)
  homs$Chromosome = factor(homs$Chromosome, levels = levels(subd$chromosome))
  h = RangedData(IRanges(start = homs$Start, end = homs$End),  
                 space = homs$Chromosome, universe = 'hg19') 
  d = RangedData(IRanges(start = subd$position, width = 1),  
                 space = subd$chromosome, universe = 'hg19')
  hits = findOverlaps(query = h, subject = d)
  for (ch in levels(subd$chromosome)){
    myhits = subjectHits(hits[[ch]])
    ix = which(subd$chromosome == ch)
    subd$cn2[ix[myhits]] = 'hom' 
  }
  
  hets$Chromosome = factor(hets$Chromosome, levels = levels(subd$chromosome))
  h = RangedData(IRanges(start = hets$Start, end = hets$End),  
                 space = hets$Chromosome, universe = 'hg19') 
  hits = findOverlaps(query = h, subject = d)
  for (ch in levels(subd$chromosome)){
    myhits = subjectHits(hits[[ch]])
    ix = which(subd$chromosome == ch)
    subd$cn2[ix[myhits]] = 'het' 
  }
  
  subd = subset(subd, cn2 %in% c('hom','het'))
  subd$baf = pmax(subd$fracB, 0)
  subd$baf = pmin(subd$fracB, 1)
  
  subd$a = subd$total*(1-subd$baf)
  subd$b = subd$total*subd$baf
  
  par(mfrow=c(2,2))
  #A0
  x = subd$a[subd$baf>.5 & subd$cn2 == 'hom']
  q = quantile(x, c(.05, .9))
  x = x[x> q[1] & x<q[2]]
  if (length(x)>1) {
    den = density(x)
    plot(den, main = 'A0: A intensities at BB')
    A0 = den$x[which.max(den$y)]
  } else {
    A0 = ifelse(length(x)>0, x, 0)
  }                    
  
  #B0
  x = subd$b[subd$baf<.5 & subd$cn2 == 'hom']
  q = quantile(x, c(.05, .9))
  x = x[x> q[1] & x<q[2]]
  if (length(x)>1) {
    den = density(x)
    plot(den, main = 'B0: B intensities at AA')
    B0 = den$x[which.max(den$y)]
  } else {
    B0 = ifelse(length(x)>0, x, 0)
  }
  
  #A1
  x = subd$a[subd$baf>.3 & subd$baf<0.7 & subd$cn2 == 'het']
  q = quantile(x, c(.05, .95))
  x = x[x> q[1] & x<q[2]]
  if (length(x)>1) {
    den = density(x)
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
    den = density(x)
    plot(den, main = 'B1: B intensities at AB')
    B1 = den$x[which.max(den$y)]
  } else {
    B1 = ifelse(length(x)>0, x, 1)
    
  }
  
  title(main = info$sample.name[i], outer = T, line = -1)
  
  ###Standardize the A and B intensities and write to file for HAPSEG
  
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


#Discard:
for (arr in 28:56){
  sample.name <- info$sample.name[arr]
  mydir = paste(prefix.out, sample.name, sep ='')
  for (j in 1:22){
    subdir = file.path(mydir,paste('chr', j, sep = ''))
    allfiles = list.files(subdir)
    tmp = unlist(lapply(strsplit(allfiles, split = '_'), function(x){x[3]}))
    nc = nchar(tmp, type = 'c')
    tmp = substr(tmp, 1, nc-4)
    ix = sort(as.numeric(tmp))
    conseq = (min(ix)):(max(ix))
    gap = setdiff(conseq, ix)
    if (length(gap)>0){
      delix = min(ix):(min(gap)-1)
      delfiles = paste('HAPSEG_SEG_', delix, '.jpg', sep = '')
      file.remove(file.path(subdir, delfiles))          
    }
  }
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
# Run ABSOLUTE
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
  RunAbsolute(seg.dat.fn, sigma.p = 0, max.sigma.h = 0.02, min.ploidy = 0.95, max.ploidy = 6, 
              primary.disease = "Breast Cancer", 
              platform = 'SNP_6.0', sample.name, results.dir, max.as.seg.count = 1500, 
              max.non.clonal = 0, 
              max.neg.genome = 0, copy_num_type = 'allelic', verbose=TRUE, maf.fn = maf.fn, 
              min.mut.af = 0.03, output.fn.base = mysample)
  sink()
}



foreach (i = c(2, 3, 12:16, 20:27, 29, 30, 34, 36, 51, 53, 56), .combine=c) %dopar% {
  DoAbsolute(myarray = info$myarray[i], mysample = info$sample.name[i], i)
}


d = paste(getwd(), '/ABSOLUTE/RunAbsolute', sep = '')
obj.name = 'ALL'
CreateReviewObject(obj.name=obj.name, absolute.files=paste(d, '/', 
                    info$sample.name[c(2, 3, 12:16, 20:27, 29, 30, 34, 36, 51, 53, 56)], 
                                                           '.ABSOLUTE.RData', sep = ''), 
                   indv.results.dir = file.path(d, obj.name), copy_num_type = "allelic", 
                   plot.modes = TRUE, verbose=TRUE)



#Not done yet:
obj.name = 'ALL'
calls.path = file.path(d, obj.name, paste(obj.name, ".PP-calls_tab.txt", sep = ''))
modes.path = file.path(d, obj.name, paste(obj.name, ".PP-modes.data.RData", sep = ''))
output.path = file.path(d, obj.name)
ExtractReviewedResults(reviewed.pp.calls.fn = calls.path, analyst.id = "test", 
                       modes.fn = modes.path, out.dir.base = output.path, 
                       obj.name = "default", copy_num_type="allelic")



################################################################
################################################################
#
# Check absolute solution
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
  
  n = 4 ###MANUALLY INSERT THE NUMBER OF CLUSTERS YOU LIKE HERE:
  pos = locator(n)
  
  ###FOR EACH CLUSTER (CHANGE c BETWEEN 1, 2, 3, ..., n)
  c = 4
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
library(IRanges)
library(aroma.affymetrix)
library(IRanges)
col.transp <- function(col1, alpha=0.5){
  tmp=as.list(t(col2rgb(col1)))
  tmp$alpha=alpha*255
  tmp$maxColorValue=255
  do.call('rgb',tmp)
}
cols = c(2:7, 'gold' , colors()[c(96)])
for (c in 7:length(cols)) cols[c] = col.transp(cols[c])

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

#A1 vs A2
x = (mat$A2)
y = (mat$A1)
ix = !(y %in% c(Inf, -Inf) | x %in% c(Inf, -Inf))
size = mat$length/sum(mat$length)*100
Wmod = lm(size~length, data = mat)
ranx = quantile(x[ix], prob = c(.05, .95), na.rm = T)
ranx = ranx + c(-1, 1)*diff(ranx)/10
rany = quantile(y[ix], prob = c(.05, .95), na.rm = T)
rany = rany + c(-1, 1)*diff(rany)/10
plot(x, y, xlim = ranx, ylim = rany, cex = size, xlab = 'Homologue 2 copy value',
     ylab = 'Homologue 1 copy value', main = 'Homologue Specific Copy Values')
ix = which(clust == 'cluster het')
points(x[ix], y[ix], col =cols[8], pch = 16, cex = size[ix])
try({ix = which(clust == 'cluster hom')
     points(x[ix], y[ix], col =cols[7], pch = 16, cex = size[ix])}, silent = T)
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
d = file.path(getwd(), 'ABSOLUTE','ScaledAbsolute')
dev.print(png, file=file.path(d, paste(info$sample.name[i], '.............png', sep = '')), 
          width=640, height=640)
#This plot, too has the problem with CN 2 clusters not on a vertical line.

################################################################
######Plot HAPSEG HSCRs coloured by the chosen clusters
library(IRanges)
col.transp <- function(col1, alpha=0.5){
  tmp=as.list(t(col2rgb(col1)))
  tmp$alpha=alpha*255
  tmp$maxColorValue=255
  do.call('rgb',tmp)
}
cols = c(2:7, 'gold' , colors()[c(96)])
for (c in 7:length(cols)) cols[c] = col.transp(cols[c])

for (i in 1:56){
  
  d = file.path(getwd(), 'ABSOLUTE', 'Choose_solution2')
  pdf(paste(d, '/', info$sample.name[i],'_Solutions_output', '.pdf', sep=''), 
      paper = 'a4', height = 12, width = 8)  
  
  
  ### Load data for this sample
  d = file.path(getwd(), 'ABSOLUTE', 'HAPSEG', info$sample.name[i])
  f = load(file.path(d, paste(info$sample.name[i], '_', info$arrname[i], '_segdat.RData', sep = '')))
  ab = as.data.frame(seg.dat$allele.segs)
  ab$Chromosome = factor(ab$Chromosome, levels = as.character(1:22))
  names(ab)[6:7] = c('A1','A2')
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
  
  ##Index cluster overlaps with absolute segments
  for (c in seq_along(unique(cl$cluster))){
    cc = subset(cl, cluster == unique(cl$cluster)[c])
    ab[,paste('cluster', unique(cl$cluster)[c])] = 0
    for (ch in seq_along(levels(cc$Chromosome))){
      ab.ix = which(ab$Chromosome == levels(cc$Chromosome)[ch])
      cc.ix = which(cc$Chromosome == levels(cc$Chromosome)[ch])
      IRab = IRanges(start = ab$Start.bp[ab.ix], end = ab$End.bp[ab.ix])
      IRcc = IRanges(start = cc$Start[cc.ix], end = cc$End[cc.ix])
      if (length(IRcc)>0 & length(IRab)>0){
        cov = coverage(c(IRab, IRcc)) - 1
        ab[ab.ix,paste('cluster', unique(cl$cluster)[c])] = viewSums(Views(cov, IRab))
      } 
    }  
  }
  
  
  ### Plots of HAPSEG output
  ###
  par(mfrow = c(3,2))
  
  #imba by logratio from hscrs
  imba = 2*abs(ab$A1/(ab$A1+ab$A2) - 0.5) 
  lr = log(ab$A1+ab$A2)
  size = ab$length/sum(ab$length)*100
  Wmod = lm(size~length, data = ab)
  ran = quantile(lr, prob = c(.05, .95))
  ran = ran + c(-1, 1)*diff(ran)/10
  plot(lr, imba, xlim = ran, cex = size, xlab = 'log(Total Copy Number)',
       ylab = 'Alleleic imbalance 2*|BAF - 0.5|', main = 'HAPSEG Allelic Imbalance')
  ix = which(ab[,'cluster het']>0)
  points(lr[ix], imba[ix], col =cols[8], pch = 16, cex = size[ix])
  try({ix = which(ab[,'cluster hom']>0)
       points(lr[ix], imba[ix], col =cols[7], pch = 16, cex = size[ix])}
      , silent = T)
  for (c in seq_along(unique(cl$cluster))){
    clust = unique(cl$cluster)[c]
    if (!(clust %in% c('het','hom'))){
      ix = which(ab[,paste('cluster', clust)]>0)
      clust.W = predict(Wmod, data.frame(length=ab[,paste('cluster', clust)]))
      points(lr[ix], imba[ix], col =cols[c], 
             pch = 16, cex = clust.W[ix])
    }
  }
  grid()
  clusters = unique(cl$cluster)
  nc = length(clusters)
  clusters = clusters[c(nc-1, nc, 1:(nc-2))]
  lcols = cols[c(7, 8, 1:(nc-2))]
  legend('bottomright', col = lcols, pch = 16, legend = c('CN2 Hom','CN2 Het',
                                                          paste('Cluster', 1:(nc-2))), bty = 'n')
  
  
  
  #A1 vs A2
  y = (ab$A1)
  x = (ab$A2)
  ix = !(y %in% c(Inf, -Inf) | x %in% c(Inf, -Inf))
  ranx = quantile(x[ix], prob = c(.05, .95))
  ranx = ranx + c(-1, 1)*diff(ranx)/10
  rany = quantile(y[ix], prob = c(.05, .95))
  rany = rany + c(-1, 1)*diff(rany)/10
  plot(x, y, xlim = ranx, ylim = rany, cex = size, xlab = 'Allele 2 logvalue',
       ylab = 'Allele 1 logvalue', main = 'HAPSEG Copy Values')
  ix = which(ab[,'cluster het']>0)
  points(x[ix], y[ix], col =cols[8], pch = 16, cex = size[ix])
  try({ix = which(ab[,'cluster hom']>0)
       points(x[ix], y[ix], col =cols[7], pch = 16, cex = size[ix])}, silent = T)
  for (c in seq_along(unique(cl$cluster))){
    clust = unique(cl$cluster)[c]
    if (!(clust %in% c('het','hom'))){
      ix = which(ab[,paste('cluster', clust)]>0)
      clust.W = predict(Wmod, data.frame(length=ab[,paste('cluster', clust)]))
      points(x[ix], y[ix], col =cols[c], pch = 16,, cex = clust.W[ix])
    }
  }
  grid()
  clusters = unique(cl$cluster)
  nc = length(clusters)
  clusters = clusters[c(nc-1, nc, 1:(nc-2))]
  lcols = cols[c(7, 8, 1:(nc-2))]
  legend('bottomright', col = lcols, pch = 16, legend = c('CN2 Hom','CN2 Het',
                                                          paste('Cluster', 1:(nc-2))), bty = 'n')
  
  
  ### Plots of histograms
  myab = ab
  myab$W = size
  myab$color.by = 8
  nclust = length(unique(cl$cluster))
  colvec = c(as.numeric(cols[c(1:(nclust-2))]), 6, 7)
  for (c in seq_along(unique(cl$cluster))){
    clust = unique(cl$cluster)[c]
    mycl = myab[,paste('cluster', clust)]
    if (clust == 'het') clust = nclust - 1
    if (clust == 'hom') clust = nclust
    
    clust.W = predict(Wmod, data.frame(length=mycl))
    clust.W[mycl == 0] = 0
    myab$W = pmax(0, myab$W - clust.W)
    ix = mycl>0
    toadd = myab[ix,]
    toadd$W = clust.W[ix]
    toadd$color.by= colvec[c]
    if (clust == nclust-1) toadd$color.by = colvec[nclust-1]
    if (clust == nclust) toadd$color.by = colvec[nclust]
    myab = rbind(myab, toadd)
  }
  
  x = c(myab$A1, myab$A2)
  l = rep(myab$W, 2)
  colseg = rep(myab$color.by, 2)
  
  ABSOLUTE:::PlotSeglenHist(x, l, color.by = colseg, use.pal = 2:8, 
                            ylab = '', 
                            data.whiskers = F,
                            xlab = 'Pooled haplytype copy values', xlim = c(0, 3))
  grid()
  clusters = unique(cl$cluster)
  nc = length(clusters)
  clusters = clusters[c(nc-1, nc, 1:(nc-2))]
  lcols = colvec[c(nc, nc-1, 1:(nc-2))]
  legend('topright', col = lcols, pch = 16, legend = c('CN2 Hom','CN2 Het',
                                                       paste('Cluster', 1:(nc-2))), bty = 'n')
  
  x = myab$A1
  l = myab$W
  colseg = myab$color.by
  
  ABSOLUTE:::PlotSeglenHist(x, l, color.by = colseg, use.pal = 2:8, ylab = '', 
                            data.whiskers = F,
                            xlab = 'Allele 1 copy values', xlim = c(0, 3))
  grid()
  
  frame()
  x = myab$A2
  l = myab$W
  colseg = myab$color.by
  
  ABSOLUTE:::PlotSeglenHist(x, l, color.by = colseg, use.pal = 2:8, ylab = '', 
                            data.whiskers = F,
                            xlab = 'Allele 2 copy values', xlim = c(0, 3))
  grid() 
  title(main = info$sample.name[i], outer = T, line = -1)
  dev.off()
}




################################################################
################################################################
#
# Sample 15
#
################################################################
################################################################
info = read.delim(paste(prefix.ann, '../../Infokey.txt', sep = ''), stringsAsFactors = F)
info$sample.name = paste(info$ind, info$sample, sep = '.')
library(IRanges)
library(Rsolnp)
library(limma)
source(file.path(getwd(), 'programs', 'make_HAPSEG_seg_obj.R'))

i = 15
d = file.path(getwd(), 'ABSOLUTE', 'HAPSEG2')
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
d = file.path(getwd(), 'ABSOLUTE', 'HAPSEG', info$sample.name[i])
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

#Define data for optimization
y = (mat$A1)
x = (mat$A2)
xclust = rep(NA, nrow(mat))
yclust = rep(NA, nrow(mat))

ix = which(clust == 'cluster hom')
xclust[ix] = 2
yclust[ix] = 0

ix = which(clust == 'cluster het')
xclust[ix] = 1
yclust[ix] = 1

ix = which(clust == 'cluster 1')
xclust[ix] = 2
yclust[ix] = 1

ix = which(clust == 'cluster 2')
xclust[ix] = 2
yclust[ix] = 2

ix = which(clust == 'cluster 3')
yclust[ix] = 2

#Save the cluster data to use
dd = data.frame(x=x, y=y, xclust = xclust, yclust=yclust, xw=w/sum(w), yw = w/sum(w),
                myclust = clust)
dd = subset(dd, !(is.na(xclust) & is.na(yclust)))


##########################################################
#Optimize one issue at the time


### Scale x
###
minx = function(pars, X, Y, xclust, xw, myclust){
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
miny = function(pars, X, Y, yclust, yw, myclust){
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

minxy = function(pars, X, Y, xclust, yclust, xw, yw, mat, myclust){
  A = pars[1]
  B = pars[2]
  C = pars[3]
  D = pars[4]
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
eqfunc = function(pars, X, Y, xclust, yclust, xw, yw, mat, myclust){
  A = pars[1]
  B = pars[2]
  C = pars[3]
  D = pars[4]
  a2 = A + B*mat$A2new
  a1 = C + D*mat$A1new
  sum(c(a1, a2)* rep(mat$W, 2)/2)
}
e.cr1 = sum(c(mat$A1, mat$A2)* rep(mat$W, 2)/2)
pars = solnp(pars = c(0,1, 0, 1), fun = minxy, myclust = dd$myclust, mat = mat,
             eqfun = eqfunc, eqB = 1,
             X=dd$x, Y=dd$y, xclust = dd$xclust, xw = dd$xw, yclust = dd$yclust, yw = dd$yw)$pars
A = pars[1]
B = pars[2]
C = pars[3]
D = pars[4]
mat$A2last = A + B*mat$A2new
mat$A1last = C + D*mat$A1new
param = c(param, A, B, C, D)
names(param) = c('a','b','c','d', 'A','B','C','D')
dd$x = A + B*dd$x
dd$y = C + D*dd$y


#A different attempt without changeing constant:
minxy = function(pars, X, Y, xclust, yclust, xw, yw, mat, myclust){
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
eqfunc = function(pars, X, Y, xclust, yclust, xw, yw, mat, myclust){
  A = 0
  B = pars[1]
  C = 0
  D = pars[2]
  a2 = A + B*mat$A2new
  a1 = C + D*mat$A1new
  sum(c(a1, a2)* rep(mat$W, 2)/2)
}
e.cr1 = sum(c(mat$A1, mat$A2)* rep(mat$W, 2)/2)
pars = solnp(pars = c(1, 1), fun = minxy, myclust = dd$myclust, mat = mat,
             eqfun = eqfunc, eqB = 1,
             X=dd$x, Y=dd$y, xclust = dd$xclust, xw = dd$xw, yclust = dd$yclust, yw = dd$yw)$pars
B = pars[1]
D = pars[2]
mat$A2last =  B*mat$A2new
mat$A1last =  D*mat$A1new
param = c(param, A, B, C, D)
names(param) = c('a','b','c','d', 'B','D')
dd$x = B*dd$x
dd$y = D*dd$y



###########################################################
# Figures of outcome
col.transp <- function(col1, alpha=0.5){
  tmp=as.list(t(col2rgb(col1)))
  tmp$alpha=alpha*255
  tmp$maxColorValue=255
  do.call('rgb',tmp)
}
cols = c(2:7, 'gold' , colors()[c(96)])
for (c in 7:length(cols)) cols[c] = col.transp(cols[c])

#A1 vs A2
x = (mat$A2last)
y = (mat$A1last)
ix = !(y %in% c(Inf, -Inf) | x %in% c(Inf, -Inf))
size = mat$length/sum(mat$length)*100
Wmod = lm(size~length, data = mat)
ranx = quantile(x[ix], prob = c(.05, .95))
ranx = ranx + c(-1, 1)*diff(ranx)/10
rany = quantile(y[ix], prob = c(.05, .95))
rany = rany + c(-1, 1)*diff(rany)/10
plot(x, y, xlim = ranx, ylim = rany, cex = size, xlab = 'Homologue 2 copy value',
     ylab = 'Homologue 1 copy value', main = 'Homologue Specific Copy Values')
ix = which(clust == 'cluster het')
points(x[ix], y[ix], col =cols[8], pch = 16, cex = size[ix])
try({ix = which(clust == 'cluster hom')
     points(x[ix], y[ix], col =cols[7], pch = 16, cex = size[ix])}, silent = T)
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
d = file.path(getwd(), 'ABSOLUTE','ScaledAbsolute2')
dev.print(png, file=file.path(d, paste(info$sample.name[i], '_ScaleAll.png', sep = '')), 
          width=640, height=640)



#Histograms
myab = mat
myab$W = size
myab$color.by = 8
nclust = length(unique(cl$cluster))
colvec = c(as.numeric(cols[c(1:(nclust-2))]), 6, 7)
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
  toadd$color.by= colvec[c]
  if (myclust == nclust-1) toadd$color.by = colvec[nclust-1]
  if (myclust == nclust) toadd$color.by = colvec[nclust]
  myab = rbind(myab, toadd)
}

x = c(myab$A1last, myab$A2last)
l = rep(myab$W, 2)
colseg = rep(myab$color.by, 2)

ABSOLUTE:::PlotSeglenHist(x, l, color.by = colseg, use.pal = 2:8, 
                          ylab = '', 
                          data.whiskers = F,
                          xlab = 'Pooled haplytype copy values', xlim = c(0, 3))
grid()
clusters = unique(cl$cluster)
nc = length(clusters)
clusters = clusters[c(nc-1, nc, 1:(nc-2))]
lcols = colvec[c(nc, nc-1, 1:(nc-2))]
legend('topright', col = lcols, pch = 16, legend = c('CN2 Hom','CN2 Het',
                                                     paste('Cluster', 1:(nc-2))), bty = 'n')
d = file.path(getwd(), 'ABSOLUTE','ScaledAbsolute')
dev.print(png, file=file.path(d, paste(info$sample.name[i], '_HistScaleAll.png', sep = '')), 
          width=640, height=640)


x = myab$A1
l = myab$W
colseg = myab$color.by

ABSOLUTE:::PlotSeglenHist(x, l, color.by = colseg, use.pal = 2:8, ylab = '', 
                          data.whiskers = F,
                          xlab = 'Allele 1 copy values', xlim = c(0, 3))
grid()

frame()
x = myab$A2
l = myab$W
colseg = myab$color.by

ABSOLUTE:::PlotSeglenHist(x, l, color.by = colseg, use.pal = 2:8, ylab = '', 
                          data.whiskers = F,
                          xlab = 'Allele 2 copy values', xlim = c(0, 3))
grid() 

  

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
d = file.path(getwd(), 'ABSOLUTE', 'HAPSEG2', info$sample.name[i])
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
    mode.res = ABSOLUTE:::FitSample(seg.dat, mut, kQ, pi_theta_qz, sigma.h, 
                                    kTauDom, kSigmaHDom, chr.arms.dat, verbose = verbose)
    mode.res = ABSOLUTE:::apply_subclonal_scna_model(seg.dat, mode.res, 
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
      mode.res = ABSOLUTE:::ApplySomaticMutsModel(mode.res, obs_scna = seg.dat$obs_SCNA, 
                                                  mut.cn.dat, pi.som.theta.q = kPiSomThetaQ, 
                                                  mut_class_w, Q=kQ, verbose = verbose)
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


myMakeSegObj = function(seg.dat.fn, min_probes = NA, 
             max_sd = NA, filter_segs = FALSE, verbose = FALSE){
  load(seg.dat.fn)
  seg.dat
}


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



#Not done yet:
obj.name = 'ALL'
calls.path = file.path(d, obj.name, paste(obj.name, ".PP-calls_tab.txt", sep = ''))
modes.path = file.path(d, obj.name, paste(obj.name, ".PP-modes.data.RData", sep = ''))
output.path = file.path(d, obj.name)
ExtractReviewedResults(reviewed.pp.calls.fn = calls.path, analyst.id = "test", 
                       modes.fn = modes.path, out.dir.base = output.path, 
                       obj.name = "default", copy_num_type="allelic")

################################################################
################################################################
#
# Sample 15 rescale again
#
################################################################
################################################################
info = read.delim(paste(prefix.ann, '../../Infokey.txt', sep = ''), stringsAsFactors = F)
info$sample.name = paste(info$ind, info$sample, sep = '.')
library(IRanges)
library(Rsolnp)
library(limma)
myMakeSegObj = function(seg.dat.fn, min_probes = NA, 
                        max_sd = NA, filter_segs = FALSE, verbose = FALSE){
  load(seg.dat.fn)
  seg.dat
}

i = 15
d = file.path(getwd(), 'ABSOLUTE', 'ScaledAbsolute2')
d = paste(d,'/', 15, sep = '')
load(file.path(d, "15.PP-modes.data.RData"))
modpars = segobj.list[[1]]$mode.res$mode.tab[2,]
b = modpars['b']
delta = modpars['delta']
newalpha = .52
newtau = modpars['genome mass']
D = newalpha*newtau + 2*(1-newalpha)
newb = 2*(1-newalpha)/D
newdelta = newalpha/D

d = file.path(getwd(), 'ABSOLUTE', 'HAPSEG2')
d = paste(d,'/', info$sample.name[i], sep = '')
seg.dat.fn = file.path(d, paste(info$sample.name[i], '_', info$arrname[i], '_scaled_segdat.RData', sep = ''))
seg.dat = myMakeSegObj(seg.dat.fn, filter_segs=FALSE, min_probes=10, max_sd=100, 
                           verbose=TRUE)  
as.seg.dat = seg.dat$as.seg.dat
nr = as.integer(nrow(as.seg.dat))
mat = as.seg.dat[1:(nr/2),]
names(mat)[6] = 'A1'
mat$A2 = as.seg.dat$copy_num[(nr/2+1):(nr)]
mat$W = mat$W*2



#############################################################
#Rescale and print new dataset

seg.dat$allele.segs[,'A1.Seg.CN'] = seg.dat$allele.segs[,'A1.Seg.CN'] * newdelta/modpars['delta'] - 
  modpars['b']*newdelta/modpars['delta'] + newb
seg.dat$allele.segs[,'A2.Seg.CN'] = seg.dat$allele.segs[,'A2.Seg.CN'] * newdelta/modpars['delta'] - 
  modpars['b']*newdelta/modpars['delta'] + newb
ix = seg.dat$allele.segs[,'A1.Seg.CN']>0 & seg.dat$allele.segs[,'A2.Seg.CN']>0
seg.dat$allele.segs = seg.dat$allele.segs[ix,]
as.seg.dat$copy_num = as.seg.dat$copy_num * newdelta/modpars['delta'] - modpars['b']*newdelta/modpars['delta'] + newb
seg.dat$as.seg.dat = as.seg.dat[c(ix, ix),]
seg.dat$seg.info[,'copy_num'] = seg.dat$seg.info[,'copy_num'] * newdelta/modpars['delta'] - 
  modpars['b']*newdelta/modpars['delta'] + newb
seg.dat$seg.info = seg.dat$seg.info[seg.dat$seg.info[,'copy_num']>0,]



mat$A2new = mat$A2* newdelta/modpars['delta'] - modpars['b']*newdelta/modpars['delta'] + newb
mat$A1new = mat$A1* newdelta/modpars['delta'] - modpars['b']*newdelta/modpars['delta'] + newb

#Print the new object
d = file.path(getwd(), 'ABSOLUTE', 'HAPSEG2', info$sample.name[i])
save(seg.dat, file = file.path(d, paste(info$sample.name[i], '_', info$arrname[i], 
                                        '_scaled2_segdat.RData', sep = '')))


###########################################################
# Figures of outcome

###Read selected segments
d = file.path(getwd(), 'ABSOLUTE', 'HAPSEG', info$sample.name[i])
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




col.transp <- function(col1, alpha=0.5){
  tmp=as.list(t(col2rgb(col1)))
  tmp$alpha=alpha*255
  tmp$maxColorValue=255
  do.call('rgb',tmp)
}
cols = c(2:7, 'gold' , colors()[c(96)])
for (c in 7:length(cols)) cols[c] = col.transp(cols[c])

#A1 vs A2
x = (mat$A2new)
y = (mat$A1new)
ix = !(y %in% c(Inf, -Inf) | x %in% c(Inf, -Inf))
size = mat$length/sum(mat$length)*100
Wmod = lm(size~length, data = mat)
ranx = quantile(x[ix], prob = c(.05, .95))
ranx = ranx + c(-1, 1)*diff(ranx)/10
rany = quantile(y[ix], prob = c(.05, .95))
rany = rany + c(-1, 1)*diff(rany)/10
plot(x, y, xlim = ranx, ylim = rany, cex = size, xlab = 'Haplotype 2 copy value',
     ylab = 'Haplotype 1 copy value', main = 'HAPSEG Copy Values')
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
d = file.path(getwd(), 'ABSOLUTE','ScaledAbsolute2')
dev.print(png, file=file.path(d, paste(info$sample.name[i], '_ScaleALL3.png', sep = '')), 
          width=640, height=640)



#Histograms
myab = mat
myab$W = size
myab$color.by = 8
nclust = length(unique(cl$cluster))
colvec = c(as.numeric(cols[c(1:(nclust-2))]), 6, 7)
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
  toadd$color.by= colvec[c]
  if (myclust == nclust-1) toadd$color.by = colvec[nclust-1]
  if (myclust == nclust) toadd$color.by = colvec[nclust]
  myab = rbind(myab, toadd)
}

x = c(myab$A1new, myab$A2new)
ix = which(x>0)
x = x[ix]
l = rep(myab$W, 2)[ix]
colseg = rep(myab$color.by, 2)[ix]

ABSOLUTE:::PlotSeglenHist(x, l, color.by = colseg, use.pal = 2:8, 
                          ylab = '', 
                          data.whiskers = F,
                          xlab = 'Pooled haplytype copy values', xlim = c(0,3))
grid()
clusters = unique(cl$cluster)
nc = length(clusters)
clusters = clusters[c(nc-1, nc, 1:(nc-2))]
lcols = colvec[c(nc, nc-1, 1:(nc-2))]
legend('topright', col = lcols, pch = 16, legend = c('CN2 Hom','CN2 Het',
                                                    paste('Cluster', 1:(nc-2))), bty = 'n')
d = file.path(getwd(), 'ABSOLUTE','ScaledAbsolute2')
dev.print(png, file=file.path(d, paste(info$sample.name[i], '_HistScaledALL3.png', sep = '')), 
          width=640, height=640)



###Junk


minbdelta = function(pars, X, Y, xclust, yclust, xw, yw, mat, myclust,
                     newb, newdelta){
  val = c(X, Y)
  w = c(xw, yw)
  mcl = c(xclust, yclust)
  A = pars[1]
  B = pars[2]
  newval = A + B*val
  e = newb + mcl*2*newdelta
  ix = !is.na(e)
  sum((newval[ix]-e[ix])^2*w[ix])
}

eqfunc = function(pars, X, Y, xclust, yclust, xw, yw, mat, myclust, newb, newdelta){
  val = c(mat$A1, mat$A2)
  w = c(mat$W, mat$W)
  A = pars[1]
  B = pars[2]
  newval = A + B*val
  sum(newval* w/sum(w))
}
pars = solnp(pars = c(0,1), fun = minbdelta, myclust = dd$myclust, mat = mat,
             eqfun = NULL, eqB = 1, newb = newb, newdelta = newdelta,
             X=dd$x, Y=dd$y, xclust = dd$xclust, xw = dd$xw, yclust = dd$yclust, yw = dd$yw)$pars
A = pars[1]
B = pars[2]
mat$A2new = A + B*mat$A2
mat$A1new = A + B*mat$A1



################################################################
################################################################
#
# Run ABSOLUTE with scaled again haplotype copy values
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
    mode.res = ABSOLUTE:::FitSample(seg.dat, mut, kQ, pi_theta_qz, sigma.h, 
                                    kTauDom, kSigmaHDom, chr.arms.dat, verbose = verbose)
    mode.res = myapply_subclonal_scna_model(seg.dat, mode.res, 
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
      mode.res = ABSOLUTE:::ApplySomaticMutsModel(mode.res, obs_scna = seg.dat$obs_SCNA, 
                                                  mut.cn.dat, pi.som.theta.q = kPiSomThetaQ, 
                                                  mut_class_w, Q=kQ, verbose = verbose)
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


myMakeSegObj = function(seg.dat.fn, min_probes = NA, 
                        max_sd = NA, filter_segs = FALSE, verbose = FALSE){
  load(seg.dat.fn)
  seg.dat
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


DoAbsolute <- function(myarray, mysample, i) {
  registerDoSEQ()
  sample.name <- mysample
  maf.fn = paste(getwd(), 
                 '/AROMA/rawData/responsify/Exome/MAF/ALL_130516_Filter_CBS_MinDepth_20', '_', 
                 mysample, '.maf', sep = '')
  seg.dat.fn <- paste(getwd(), '/ABSOLUTE/HAPSEG2/', mysample,'/',mysample,'_',myarray,
                      "_scaled3_segdat.RData", sep="")
  results.dir <- paste(getwd(), '/ABSOLUTE/Scaled3Absolute', sep = '')
  print(paste("Starting myarray", mysample, "at", results.dir))
  log.dir <- paste(getwd(), '/ABSOLUTE/Scaled3Absolute/abs_logs', sep = '') 
  if (!file.exists(log.dir)) {
    dir.create(log.dir, recursive=TRUE)
  }
  if (!file.exists(results.dir)) {
    dir.create(results.dir, recursive=TRUE)
  }
  sink(file=file.path(log.dir, paste(mysample, ".abs.out.txt", sep="")))
  myRunAbsolute(seg.dat.fn, sigma.p = 0, max.sigma.h = 0.02, min.ploidy = 2.5, max.ploidy = 3.5, 
                primary.disease = "Breast Cancer", 
                platform = 'SNP_6.0', sample.name, results.dir, max.as.seg.count = 1500, 
                max.non.clonal = 0, 
                max.neg.genome = 0, copy_num_type = 'allelic', verbose=TRUE, maf.fn = maf.fn, 
                min.mut.af = 0.03, output.fn.base = mysample)
  sink()
}


i = 15
DoAbsolute(myarray = info$myarray[i], mysample = info$sample.name[i], i)


d = paste(getwd(), '/ABSOLUTE/Scaled3Absolute', sep = '')
obj.name = '15'
CreateReviewObject(obj.name=obj.name, absolute.files=paste(d, '/', 
                                                           info$sample.name[15], '.ABSOLUTE.RData', sep = ''), 
                   indv.results.dir = file.path(d, obj.name), copy_num_type = "allelic", 
                   plot.modes = TRUE, verbose=TRUE)

obj.name = '15'
d = paste(getwd(), '/ABSOLUTE/Scaled3Absolute', sep = '')
calls.path = file.path(d, obj.name, paste(obj.name, ".PP-calls_tab.txt", sep = ''))
modes.path = file.path(d, obj.name, paste(obj.name, ".PP-modes.data.RData", sep = ''))
output.path = file.path(d, obj.name)
ExtractReviewedResults(reviewed.pp.calls.fn = calls.path, analyst.id = "IL", 
                       modes.fn = modes.path, out.dir.base = output.path, 
                       obj.name = "AfterScaling3", copy_num_type="allelic")

