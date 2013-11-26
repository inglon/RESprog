################################################################
################################################################
# Script: ProgR_5.1_Heterogeneity.R
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

prefix.raw = paste(getwd(), "/rawData/responsify/GenomeWideSNP_6/", sep='')
prefix.in = paste(getwd(), "/ABSOLUTE/APT/", sep='')
prefix.ann = paste(getwd(), "/annotationData/chipTypes/GenomeWideSNP_6/", sep='')
prefix.out = paste(getwd(), "/ABSOLUTE/HAPSEG/", sep='')
prefix.out2 = paste(getwd(), "/ABSOLUTE/HAPSEG2/", sep='')

#prefix.prog = paste(getwd(), "/programs/", sep='')
#source(paste(prefix.prog, 'ProgR_0.1_Functions.R', sep=''))
date = format(Sys.Date())



################################################################
################################################################
#
# Birdseed with Affymetrix Power Tools
#
################################################################
################################################################


################################################################
##  Move files around to be able to run APT before HAPSEG before ABSOLUTE

dr = paste(getwd(), '/rawData/responsify/GenomeWideSNP_6/', sep = '')
info$tmp = paste(info$arrname,'.CEL', sep = '')
files = list.files(dr)
nc = nchar(files, type = 'c')
files = files[substr(files, nc-2, nc) == 'CEL']
dir.create(paste(dr, 'notPassed', sep = ''), recursive = T, showWarnings = F)
for (i in 1:nrow(info)){
  index = which(files == info$tmp[i])
  if (files[i] != 'notPassed') if (!info$pass[i]) file.rename(paste(dr, files[i], sep = ''), 
                                                              paste(dr, 'notPassed/',files[i], sep = ''))
}

#Now rename arrays so that they have '.', not '-':
dr = paste(getwd(), '/rawData/responsify/GenomeWideSNP_6/', sep = '')
files = list.files(dr)
nc = nchar(files, type = 'c')
files = files[substr(files, nc-2, nc) == 'CEL']
newfiles = unlist(lapply(strsplit(files, split = '-', fixed = T), paste, collapse = '.'))
for (i in 1:length(files)){
  file.rename(paste(dr, files[i], sep = ''), paste(dr, newfiles[i], sep = ''))
}


################################################################
##  Run APT from anywhere


#Run apt-probeset-summarize from RESPONSIFY/HAPSEG
apt-probeset-summarize --cdf-file  ~/RESPONSIFY/annotationData/chipTypes/GenomeWideSNP_6/GenomeWideSNP_6.cdf --analysis rma-bg,quant-norm.sketch=50000,pm-only,med-polish,expr.genotype=true --out-dir ~/RESPONSIFY/ABSOLUTE/APT/ ~/RESPONSIFY/rawData/responsify/GenomeWideSNP_6/*.CEL

#This file (or similar) will be your snp.fn argument to HAPSEG:
rma-bg.quant-norm.pm-only.med-polish.expr.summary.txt 

#Next you will need to run birdseed:
apt-probeset-genotype -o BirdseedOutput -c ~/RESPONSIFY/annotationData/chipTypes/GenomeWideSNP_6/GenomeWideSNP_6.cdf --set-gender-method cn-probe-chrXY-ratio --chrX-probes ~/RESPONSIFY/ABSOLUTE/APT/Birdseedfiles/GenomeWideSNP_6.chrXprobes --chrY-probes ~/RESPONSIFY/ABSOLUTE/APT/Birdseedfiles/GenomeWideSNP_6.chrYprobes --special-snps ~/RESPONSIFY/ABSOLUTE/APT/Birdseedfiles/GenomeWideSNP_6.specialSNPs --read-models-birdseed ~/RESPONSIFY/ABSOLUTE/APT/Birdseedfiles/GenomeWideSNP_6.v2.6.birdseed.models -a birdseed-dev --write-models ~/RESPONSIFY/rawData/responsify/GenomeWideSNP_6/*.CEL

#These will be input to HAPSEG:
birdseed-dev.calls.txt will be calls.fn
birdseed-dev.snp-models.txt which will be clusters.fn  

#Also let
snp.file.parser=AptSnpFileParser
cluster.file.parser=BirdseedClustersFileParser

################################################################
##  Run APT from anywhere ATTEMPT 2


#Run apt-probeset-summarize from RESPONSIFY/HAPSEG
apt-probeset-summarize --cdf-file  ~/RESPONSIFY/annotationData/chipTypes/GenomeWideSNP_6/GenomeWideSNP_6.cdf --analysis quant-norm.sketch=50000,pm-only,med-polish,expr.genotype=true --out-dir ~/RESPONSIFY/ABSOLUTE/APT2/ ~/RESPONSIFY/rawData/responsify/GenomeWideSNP_6/*.CEL

#This file (or similar) will be your snp.fn argument to HAPSEG:
quant-norm.pm-only.med-polish.expr.summary.txt 

#Next you will need to run birdseed:
apt-probeset-genotype -o BirdseedOutput -c ~/RESPONSIFY/annotationData/chipTypes/GenomeWideSNP_6/GenomeWideSNP_6.cdf --set-gender-method cn-probe-chrXY-ratio --chrX-probes ~/RESPONSIFY/ABSOLUTE/APT/Birdseedfiles/GenomeWideSNP_6.chrXprobes --chrY-probes ~/RESPONSIFY/ABSOLUTE/APT/Birdseedfiles/GenomeWideSNP_6.chrYprobes --special-snps ~/RESPONSIFY/ABSOLUTE/APT/Birdseedfiles/GenomeWideSNP_6.specialSNPs --read-models-birdseed ~/RESPONSIFY/ABSOLUTE/APT/Birdseedfiles/GenomeWideSNP_6.v2.6.birdseed.models -a birdseed-v1 --write-models ~/RESPONSIFY/rawData/responsify/GenomeWideSNP_6/*.CEL

#These will be input to HAPSEG:
birdseed-v1.calls.txt will be calls.fn
birdseed-v1.snp-models.txt which will be clusters.fn  

#Also let
snp.file.parser=AptSnpFileParser
cluster.file.parser=BirdseedClustersFileParser


################################################################
################################################################
#
# HAPSEG
#
################################################################
################################################################
library(HAPSEG)
source(paste(getwd(), '/programs/Hapseg.r', sep = ''))
library(foreach)
library(doMC)
registerDoMC(20)
info = read.delim(paste(prefix.ann, '../../Infokey.txt', sep = ''))
info$sample = substr(info$Sample.ID, 1, 4)
info$myarray = info$arrname

DoHapseg <- function(myarray, mysample) {
  registerDoSEQ()
  plate.name <- mysample
  clusters.fn <- paste(getwd(), "/ABSOLUTE/APT/BirdseedOutput/birdseed-dev.snp-models.txt", sep = '')  
  phased.bgl.dir <- "/usr/local/bioinfsoftware/R/current/lib64/R/library/ABSOLUTE/etc/phasedBGL/hg19"
 
  snp.fn <- paste(prefix.in,"rma-bg.quant-norm.pm-only.med-polish.expr.summary.txt", sep = '')
  results.dir <- paste(prefix.out2, myarray, sep = '')
  if (!file.exists(results.dir)) {
    dir.create(results.dir, recursive=TRUE)
  }
  print(paste("Starting array", myarray, "at", results.dir))
  log.dir <- paste(prefix.out2, 'log', sep = '')
  if (!file.exists(log.dir)) {
    dir.create(log.dir, recursive=TRUE)
  }
  sink(file=file.path(log.dir, paste(myarray, ".hapseg.out.txt", sep="")))
  myRunHapSeg(plate.name, myarray, seg.fn=NULL, snp.fn, genome = 'hg19', results.dir, platform = "SNP_6.0", 
            use.pop='CEPH', 
            impute.gt=TRUE, plot.segfit=TRUE, merge.small=TRUE, merge.close=TRUE, min.seg.size=5, 
            normal=FALSE, 
            out.p=0.001, seg.merge.thresh=1e-10, use.normal = TRUE, adj.atten = FALSE, 
            phased.bgl.dir, calibrate.data=FALSE, 
            drop.x=FALSE, verbose=TRUE, snp.file.parser=mySnpFileParser)
  sink()
}



i = 24

foreach (i=1:56, .combine=c) %dopar% {
  DoHapseg(info$myarray[i], paste(i,info$sample[i], sep = '.'))
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
prefix.in = paste(getwd(), "/Exome/", sep='')
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


###Split maf file into one file for each sample
mafbase = 'ALL_130516_Filter_CBS_MinDepth_20'
maf = read.delim(paste(prefix.in, mafbase,'.maf', sep = ''), 
                 row.names = NULL, stringsAsFactors = FALSE, 
                 check.names = FALSE, na.strings = c("NA", "---"), 
                 blank.lines.skip = TRUE, comment.char = "#")
info = read.delim(paste(prefix.ann, 'Infokey.txt', sep = ''))
info$sample = substr(as.character(info$Sample.ID), 1, 4)
for (i in 1:nrow(info)){
  samp = info$sample[i]
  out = subset(maf, Tumor_Sample_Barcode == paste(samp, 'T', sep = ''))
  if (nrow(out)>0) write.table(out, file = paste(prefix.in,'MAF/', 
              mafbase, '_', i,'.', samp, '.maf', sep = ''))
}

################################################################
################################################################
#
# Absolute TOTAL copy numbers in segments: 
# 
# 
################################################################
################################################################
library(ABSOLUTE)
library(foreach)
library(doMC)
prefix.ABS = paste(getwd(), '/ABSOLUTE/RunAbsoluteTCN/', sep = '')
info = read.delim(paste(prefix.ann, '../../Infokey.txt', sep = ''))
info$sample = substr(info$Sample.ID, 1, 4)
info$myarray = info$arrname


DoAbsoluteTCN <- function(i) {
  registerDoSEQ()
  sample.name = paste(i, info$sample[i], sep = '.')
  fitname = paste('fit', i, sep = '')
  seg.dat.fn = paste(prefix.ABS, 'seg.dat.files/', fitname, 'tcn.txt', sep = '');
  results.dir = prefix.ABS
  maf.fn = paste(getwd(), 
                 '/rawData/Exome/MAF/ALL_130516_Filter_CBS_MinDepth_20', '_', i, '.', 
                 sample.name, '.maf', sep = '')
  print(paste("Starting myarray", sample.name, "at", results.dir))
  log.dir <- paste(getwd(), '/ABSOLUTE/RunAbsoluteTCN/abs_logs', sep = '') 
  if (!file.exists(log.dir)) {
    dir.create(log.dir, recursive=TRUE)
  }
  if (!file.exists(results.dir)) {
    dir.create(results.dir, recursive=TRUE)
  }
  
  sink(file=file.path(log.dir, paste(i, '.', sample.name, ".abs.out.txt", sep="")))
  RunAbsolute(seg.dat.fn, sigma.p = 0, max.sigma.h = 0.02, min.ploidy = 2, max.ploidy = 4, 
              primary.disease = 'Breast Cancer', platform = "SNP_6.0", sample.name = sample.name, 
              results.dir = results.dir,
              max.as.seg.count = 1500, max.non.clonal = 0.5, max.neg.genome = 0, copy_num_type = 'total', 
              maf.fn=maf.fn, min.mut.af=NULL, output.fn.base=sample.name,
              verbose=T)
  #Note, min and max ploidy might have to be changed between samples
  #As would max.non.clonal and the rest...
  sink()
}

foreach (i = 1:56, .combine=c) %dopar% {
  DoAbsoluteTCN(i)
}

###Follow up sample 24
sample.name = paste(i, info$sample[i], sep = '.')
CreateReviewObject(obj.name=sample.name, absolute.files=paste(prefix.ABS,'24.4742.ABSOLUTE.RData', sep = ''), 
                   prefix.ABS, "total", verbose=TRUE)
calls.path = paste(prefix.ABS, sample.name,".PP-calls_tab.txt", sep = '')
modes.path = paste(prefix.ABS, sample.name, ".PP-modes.data.RData", sep = '')
output.path = paste(prefix.ABS, sample.name, "abs_extract", sep = '')
ExtractReviewedResults(calls.path, "IL", modes.path, output.path, "absolute", "total")

absdat = read.delim(paste(output.path,'/reviewed/SEG_MAF/24.4742.segtab.txt', sep = ''))
hist(absdat$cancer_cell_frac[!(absdat$cancer_cell_frac %in% c(0,1))], breaks = 100) #Cancer cell fractions
na.omit(absdat$cancer_cell_frac[!(absdat$cancer_cell_frac %in% c(0,1))])
[1] 0.89 0.59 0.66 0.64 0.68 0.49 0.49 0.56 0.63 0.52 0.57
na.omit(absdat$subclonal[!(absdat$subclonal %in% c(0,1))])



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

DoAbsolute <- function(myarray, mysample, i) {
  registerDoSEQ()
  sample.name <- mysample
  maf.fn = paste(getwd(), 
     '/rawData/Exome/MAF/ALL_130516_Filter_CBS_MinDepth_20', '_', i, '.', 
                 mysample, '.maf', sep = '')
  seg.dat.fn <- paste(prefix.out2, myarray,'/',i,'.',mysample,'_',myarray,
                      "_segdat.RData", sep="")
  results.dir <- paste(getwd(), '/ABSOLUTE/RunAbsolute2', sep = '')
  print(paste("Starting myarray", myarray, "at", results.dir))
  log.dir <- paste(getwd(), '/ABSOLUTE/RunAbsolute2/abs_logs', sep = '') 
  if (!file.exists(log.dir)) {
    dir.create(log.dir, recursive=TRUE)
  }
  if (!file.exists(results.dir)) {
    dir.create(results.dir, recursive=TRUE)
  }
  sink(file=file.path(log.dir, paste(i, '.', mysample, ".abs.out.txt", sep="")))
  RunAbsolute(seg.dat.fn, sigma.p = 0, max.sigma.h = 0.02, min.ploidy = 0.95, max.ploidy = 10, 
              primary.disease = "Breast Cancer", 
              platform = 'SNP_6.0', sample.name, results.dir, max.as.seg.count = 1500, 
              max.non.clonal = 0, 
              max.neg.genome = 0, copy_num_type = 'allelic', verbose=TRUE, maf.fn = maf.fn, 
              min.mut.af = min.mut.af, output.fn.base = paste(i, mysample, sep = '.'))
  sink()
}


foreach (i = 1:56, .combine=c) %dopar% {
  DoAbsolute(myarray = info$myarray[i], mysample = info$sample[i], i)
}



################################################################
################################################################
#
# Continue ABSOLUTE
#
################################################################
################################################################

obj.name <- "DRAWS_summary"
results.dir <- file.path(".", "output", "abs_summary")
absolute.files <- file.path(".", "output",
                            scans, "absolute",
                            paste(scans, ".ABSOLUTE.RData", sep=""))
library(ABSOLUTE)
CreateReviewObject(obj.name, absolute.files, results.dir, "allelic", verbose=TRUE)

## At this point you'd perform your manual review and mark up the file 
## output/abs_summary/DRAWS_summary.PP-calls_tab.txt by prepending a column with
## your desired solution calls. After that (or w/o doing that if you choose to accept
## the defaults, which is what running this code will do) run the following command:

calls.path = file.path("output", "abs_summary", "DRAWS_summary.PP-calls_tab.txt")
modes.path = file.path("output", "abs_summary", "DRAWS_summary.PP-modes.data.RData")
output.path = file.path("output", "abs_extract")
ExtractReviewedResults(calls.path, "test", modes.path, output.path, "absolute", "allelic")




################################################################
################################################################




