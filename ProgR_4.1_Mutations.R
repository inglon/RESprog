################################################################
################################################################
# Script: ProgR_4.1_Mutations.R
# Author: Ingrid Lonnstedt
# Date:  15/05/2013
# R version: R 2.15.1
# Details: Analyses of mutations and CN changes
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

prefix.in = paste(getwd(), "/AROMA/rawData/responsify/Exome/", sep='')
prefix.raw = paste(getwd(), "/AROMA/rawData/responsify/", sep='')
prefix.out = paste(getwd(), "/output/Exome/", sep='')
prefix.ann = paste(getwd(), "/AROMA/annotationData/", sep='')
date = format(Sys.Date())



################################################################
################################################################
#
# Clinical data
#
################################################################
################################################################
#

### Read clinical data
###
dat = read.delim(paste(prefix.raw, 'Clinical_data/spss working-Demos combined Leuven Bordet series-Jan2013.txt',
                       sep = ''), na.strings = c('', NA), check.names = F, dec = ',')[1:108,1:141]
index = c(94, 95, 97:99, 101, 102, 107)
tmpdate = dat$date_diagnosis[index]
dat$date_diagnosis = as.Date(as.character(dat$date_diagnosis),format='%d/%m/%Y')
dat$date_diagnosis[index] = as.Date(as.character(tmpdate),format='%Y-%m-%d')
index = c(98, 99, 103, 105) #106 set to 1951-07-01 manually
tmpdate = dat$Date_of_birth[index]
dat$Date_of_birth = as.Date(as.character(dat$Date_of_birth),format='%Y-%m-%d')
dat$Date_of_birth[index] = as.Date(as.character(tmpdate),format='%d/%m/%Y')
dat$Date_of_birth[106] = as.Date('1951-07-01',format='%Y-%m-%d')
dat$age = dat$date_diagnosis - dat$Date_of_birth 
dat$til = as.numeric(as.character(dat$Stomal_LI.))
#Warning above is OK.
dat$ID = as.character(dat$Frozen_tissue_BO_no)

tmp = read.delim(paste(prefix.raw, 'Clinical_data/HER2-lev-bordet-Aug2012-responsify.txt',
                       sep = ''), na.strings = c('', NA), check.names = F, dec = ',')
tmp$ID = as.character(dat$Frozen_tissue_BO_no)
tmp = subset(tmp, select = c('ID','IDFS_Y_N','IDFS_days'))
dat = merge(dat, tmp, sort = F, all.x = T)

#'IDFS_Y_N','IDFS_days'
################################################################
################################################################
#
# Which sample is deep and severe?
#
################################################################
################################################################

ex = read.csv(paste(prefix.in, 'ALL_130516_24samples_PERFORMANCE.csv', 
                    sep = ''), row.names = 1)
ex = as.data.frame(t(ex))
coverage = ex[,29]
hist(coverage)
summary(coverage)

data = data.frame(rownames = rownames(ex),coverage = ex[,29])
data$ID = substr(data$rownames, 2, 5)
data$TN = substr(data$rownames, 6, 6)
data= subset(data, TN == 'T')
data = merge(data, dat)
#data = subset(data, IDFS_Y_N ==1)
data$IDFS_days = as.numeric(data$IDFS_days)
#plot(coverage~IDFS_days, data = data)
#data$ID[data$coverage>120] 
#[1] "4787"
#data$coverage[data$ID == '4787']
#[1]  50.07 135.28

mut = read.delim(paste(prefix.in, 'Loi_130516_24samples/ALL_130516_Filter_CBS_MinDepth_20.txt', 
                       sep = ''))
tmp = as.data.frame(table(mut$SAMPLE))
names(tmp) = c('ID', 'nmut')
tmp$ID = substr(as.character(tmp$ID), 1, 4)
data = merge(data, tmp)

png(paste(prefix.out, 'N_mutations_and_coverage.png', sep = ''), width = 540, height = 400)
plot(nmut~coverage, data = data, type = 'n', ylab = 'Number of mutations',
     xlab = 'Average coverage', log = 'y')
text(nmut~coverage, data = data, labels = data$ID, cex = .8,
     col = data$IDFS_Y_N + 1)
text(160, 50, labels='Relapsed', col = 'red', adj = 1)
text(160, 40, labels='Not relapsed', adj = 1)
dev.off()


################################################################
################################################################
#
# CN diversity classification of samples
#
################################################################
################################################################

info = read.delim(paste(prefix.ann, 'Infokey.txt', sep = ''))
info$CNdiversity = NA
info$CNdiversity[c(11, 12, 15, 18, 21, 24, 25, 26, 27, 32, 33, 43, 44, 45, 52, 55, 9)] = 'High'
info$CNdiversity[c(1, 2, 5, 13, 17, 19, 20, 29, 31, 34, 35, 37, 38, 40, 41, 42, 46, 47,
               48, 49, 50, 51, 53, 56)] = 'Medium'
info$CNdiversity[c(3, 4, 10, 14, 16, 22, 28, 30, 36, 54, 6)] = 'Low'


################################################################
################################################################
#
# Mutation allele frequency plot
#
################################################################
################################################################
mut = read.delim(paste(prefix.in, 'AllSamplesMAF/ALL_130516_Filter_CBS_MinDepth_20.txt', 
                    sep = ''))
names(mut)[53:55] = c('Read_depth','Variant_depth','Variant_frequency')
library(lattice)
xyplot(Variant_frequency~Read_depth|SAMPLE, data = mut)
png(paste(prefix.out, 'VariantFrequency_readDepth.png', sep = ''), 
    width = 540, height = 900)
par(mfrow = c(3,1), mar = c(4,4,2,.2))
plot(Variant_frequency~Read_depth, data = mut,
     subset = SAMPLE == '4787T', main = '26.4787')
plot(Variant_frequency~Read_depth, data = mut,
     subset = SAMPLE == '4742T', main = '24.4742')
plot(Variant_frequency~Read_depth, data = mut,
     subset = SAMPLE == '4669T', main = '21.4669')
plot(Variant_frequency~Read_depth, data = mut,
     subset = SAMPLE == '4433T', main = '15.4433', log = 'x')
dev.off()
hist(mut$Variant_frequency[mut$SAMPLE == '4787T'], breaks = 1000)

hist(mut$Variant_frequency[mut$SAMPLE == '4433T'], breaks = 500, main = '15.4433',
     xlab = 'Variant_frequency')


################################################################
################################################################
#
# Sample 15
#
################################################################
################################################################
info = read.delim(paste(prefix.ann, 'Infokey.txt', sep = ''))

d = file.path(getwd(), 'ABSOLUTE', 'Scaled2Absolute', '15', 'reviewed', 'SEG_MAF')
maf = read.delim(file.path(d, "A1761.16_ABS_MAF.txt"))
seg = read.delim(file.path(d, "A1761.16.segtab.txt"))

hist(maf$cancer_cell_frac[maf$subclonal.ix == 1], breaks = 50)
ccf.a1 = seg$cancer.cell.frac.a1[seg$subclonal.a1 == 1]
ccf.a2 = seg$cancer.cell.frac.a2[seg$subclonal.a2 == 1]
hist(c(ccf.a1, ccf.a2), breaks = 50)
hist(c(ccf.a1, ccf.a2, maf$cancer_cell_frac[maf$subclonal.ix == 1]), breaks = 50)

sort(maf$Hugo_Symbol)
maf = maf[order(),]

###################################################################
## List of clonal and subclonal mutations
##

maf$Hugo_Symbol[maf$subclonal.ix == 1] #(seems to be those with prob(clonal)<prob(subclonal))
maf$Hugo_Symbol[maf$subclonal.ix == 0] #(seems to be those with prob(clonal)<prob(subclonal))

table(maf$subclonal.ix, maf$subclonal_SCNA)
boxplot(Pr_subclonal~subclonal_SCNA, data = maf)

###################################################################
## Mutation patterns
##

identical(maf$Reference_Allele, maf$Tumor_Seq_Allele1)

#Junk
maf = read.delim(paste(prefix.in, 'MAF/ALL_130516_Filter_CBS_MinDepth_20_', info$sample.name[15], '.maf', 
                       sep = ''))
maf$af = maf$t_alt_count/(maf$t_alt_count + maf$t_ref_count)
hist(maf$af, main = 'Fractions of reads with mutation', xlab = 'Allele frequency', xlim = c(0,1))
hist(maf$af, breaks = 50, main = 'Fractions of reads with mutation', xlab = 'Allele frequency', xlim = c(0,1))


hist(maf$af[maf$Chromosome == '17'], breaks = 50, main = 'Fractions of reads with mutation', 
     xlab = 'Allele frequency', xlim = c(0,1))
maf$Hugo_Symbol[maf$Chromosome == '17']

