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

prefix.raw = paste(getwd(), "/AROMA/rawData/responsify/", sep='')
prefix.norm = paste(getwd(), "/AROMA/reports/", sep='')
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
# Organize data: Clinical data and mutations from MuSic
#
################################################################
################################################################

###Read mutations
d = file.path(getwd(), 'AROMA', 'rawData', 'responsify', 'Exome', 
              'MuSic smg_MuSic sign mt.txt')
mut = read.delim(d, stringsAsFactors = F)
muts = mut$X.Gene
d = file.path(getwd(), 'AROMA', 'rawData', 'responsify','Exome', 'Mutsig_non-silent')
mut = read.delim(file.path(d,"AllVarintsMutectOncotatorPlusMutsigCategNo5214.tsv"))
ids = unique(mut$Tumor_Sample_Barcode)
ids = substr(as.character(ids), 2,5)
#setdiff(muts, mut$Hugo_Symbol)
mut = unique(subset(mut, Hugo_Symbol %in% muts, 
                    select = c('Tumor_Sample_Barcode','Hugo_Symbol')))
mut$Tumor_Sample_Barcode = as.character(mut$Tumor_Sample_Barcode)
mut$Hugo_Symbol = as.character(mut$Hugo_Symbol)
Xdata = table(mut$Tumor_Sample_Barcode, mut$Hugo_Symbol)
rownames(Xdata) = substr(rownames(Xdata), 2, 5)
toadd = setdiff(ids, rownames(Xdata))
Xdata = rbind(Xdata, matrix(0, ncol = ncol(Xdata), nrow = length(toadd), dimnames=
                list(toadd, colnames(Xdata))))
Xdata = Xdata[order(rownames(Xdata)),]

#Read clinical variables
data = data.frame(ID = substr(rownames(Xdata), 1, 4), stringsAsFactors = F)
data = merge(data, subset(dat, select = c('ID','age','tumor_size','nodes_pos','Nodal_status', 'ER_status', 'til', 'date_diagnosis')), 
             all.x = T, sort = F)
data$age_diag = as.numeric(data$age/365.25)
data$year_diag = as.numeric(substr(as.character(data$date_diagnosis), 1, 4))
identical(rownames(Xdata), data$ID)
#[1] TRUE

#Now, data and Xdata could be cbind() and put into linear models.

################################################################
################################################################
#
# Clinical data ~ mutations using MuSic mutations
#
################################################################
################################################################

#############################################
# Multiple linear models TIL ~ mutation
prefix.out = paste(getwd(),'/AROMA/reports/TILS_vs_mutations/', sep = '')

out = matrix(NA, ncol = 2, nrow = ncol(Xdata))
colnames(out) = c('Estimate', 'Pr(>|t|)')
for (i in 1:ncol(Xdata)){
  data$x = Xdata[,i]
  mod = lm(log(til+1) ~ x, data = data) 
  if ('x' %in% rownames(summary(mod)$coef))  out[i,] = summary(mod)$coef['x',colnames(out)]
}

par(mfrow = c(1,2))
hist(out[,2],breaks=c(0,0.001,0.01,0.05,0.1,0.2,0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9, 1),
     main = 'Mutations linear model P-values', xlab = 'p-values')
out = cbind(out, p.adjust(out[,2], method = 'BH'))
hist(out[,3],breaks=c(0,0.001,0.01,0.05,0.1,0.2,0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9, 1),
     main = 'Mutations linear model FDR values', xlab = 'FDR')
#dev.print(png, file=paste(prefix.out,'results/Pvalue2_hist_mutations.png', sep = ''), 
#          width=400, height=400)
qqnorm(out[,1])
#dev.print(png, file=paste(prefix.out,'results/TILcoefficients2_mutations.png', sep = ''), 
#          width=400, height=400)


#Toptable
toptable = data.frame('Hugo.Symbol' = colnames(Xdata), 
                      Coefficient = out[, 'Estimate'], 
                      p.value = out[, 'Pr(>|t|)'],
                      FDR = out[, 3])
toptable = toptable[order(toptable$p.value),]
toptable$rank = 1:nrow(toptable)
toptable = toptable[, c(ncol(toptable), 1:(ncol(toptable)-1))]
write.table(toptable, file = paste(prefix.out, 'TIL_by_mutation.txt', sep = ''), 
            row.names = F, quote = F, sep = '\t')

#############################################
# Multiple linear models TIL ~ mutation: ER+ samples only
prefix.out = paste(getwd(),'/AROMA/reports/TILS_vs_mutations/', sep = '')

out = matrix(NA, ncol = 2, nrow = ncol(Xdata))
colnames(out) = c('Estimate', 'Pr(>|t|)')
for (i in 1:ncol(Xdata)){
  data$x = Xdata[,i]
  mod = lm(log(til+1) ~ x, data = data, subset =  ER_status == 1)
  if ('x' %in% rownames(summary(mod)$coef))  out[i,] = summary(mod)$coef['x',colnames(out)]
}

par(mfrow = c(1,2))
hist(out[,2],breaks=c(0,0.001,0.01,0.05,0.1,0.2,0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9, 1),
     main = 'Mutations linear model P-values', xlab = 'p-values')
out = cbind(out, p.adjust(out[,2], method = 'BH'))
hist(out[,3],breaks=c(0,0.001,0.01,0.05,0.1,0.2,0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9, 1),
     main = 'Mutations linear model FDR values', xlab = 'FDR')
#dev.print(png, file=paste(prefix.out,'results/Pvalue2_hist_mutations.png', sep = ''), 
#          width=400, height=400)
qqnorm(out[,1])
#dev.print(png, file=paste(prefix.out,'results/TILcoefficients2_mutations.png', sep = ''), 
#          width=400, height=400)


#Toptable
toptable = data.frame('Hugo.Symbol' = colnames(Xdata), 
                      Coefficient = out[, 'Estimate'], 
                      p.value = out[, 'Pr(>|t|)'],
                      FDR = out[, 3])
toptable = toptable[order(toptable$p.value),]
toptable$rank = 1:nrow(toptable)
toptable = toptable[, c(ncol(toptable), 1:(ncol(toptable)-1))]
write.table(toptable, file = paste(prefix.out, 'TIL_by_mutation_ERpositives.txt', sep = ''), 
            row.names = F, quote = F, sep = '\t')

#############################################
# Multiple linear models TIL ~ mutation: ER- samples only
prefix.out = paste(getwd(),'/AROMA/reports/TILS_vs_mutations/', sep = '')

out = matrix(NA, ncol = 2, nrow = ncol(Xdata))
colnames(out) = c('Estimate', 'Pr(>|t|)')
for (i in 1:ncol(Xdata)){
  data$x = Xdata[,i]
  mod = lm(log(til+1) ~ x, data = data, subset =  ER_status == 0)
  if ('x' %in% rownames(summary(mod)$coef))  out[i,] = summary(mod)$coef['x',colnames(out)]
}

par(mfrow = c(1,2))
hist(out[,2],breaks=c(0,0.001,0.01,0.05,0.1,0.2,0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9, 1),
     main = 'Mutations linear model P-values', xlab = 'p-values')
out = cbind(out, p.adjust(out[,2], method = 'BH'))
hist(out[,3],breaks=c(0,0.001,0.01,0.05,0.1,0.2,0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9, 1),
     main = 'Mutations linear model FDR values', xlab = 'FDR')
#dev.print(png, file=paste(prefix.out,'results/Pvalue2_hist_mutations.png', sep = ''), 
#          width=400, height=400)
qqnorm(out[,1])
#dev.print(png, file=paste(prefix.out,'results/TILcoefficients2_mutations.png', sep = ''), 
#          width=400, height=400)


#Toptable
toptable = data.frame('Hugo.Symbol' = colnames(Xdata), 
                      Coefficient = out[, 'Estimate'], 
                      p.value = out[, 'Pr(>|t|)'],
                      FDR = out[, 3])
toptable = toptable[order(toptable$p.value),]
toptable$rank = 1:nrow(toptable)
toptable = toptable[, c(ncol(toptable), 1:(ncol(toptable)-1))]
write.table(toptable, file = paste(prefix.out, 'TIL_by_mutation_ERnegatives.txt', sep = ''), 
            row.names = F, quote = F, sep = '\t')

#############################################
# Mann Whitney U tests Nodes ~ mutation
prefix.out = paste(getwd(),'/AROMA/reports/Nodes_vs_mutations/', sep = '')
library('exactRankTests')

out = data.frame(ps = rep(NA, ncol(Xdata)), median.nodes.with.mut = NA, 
                 median.nodes.without.mut = NA, n.with.mut = NA, n.without.mut = NA)
for (i in 1:ncol(Xdata)){
  data$x = Xdata[,i]
  mod = wilcox.exact(data$nodes_pos[data$x==0], data$nodes_pos[data$x == 1])
  out$ps[i] = mod$p.value
  nmut = table(data$x)
  meds = tapply(data$nodes_pos, data$x, median, na.rm = T)
  out$median.nodes.with.mut[i] = meds['1']
  out$median.nodes.without.mut[i] = meds['0']
  out$n.with.mut[i] = nmut['1']
  out$n.without.mut[i] = nmut['0']
  
}


ps = rep(NA, ncol(Xdata))
for (i in 1:ncol(Xdata)){
  data$x = Xdata[,i]
  mod = wilcox.exact(data$nodes_pos[data$x==0], data$nodes_pos[data$x == 1])
  ps[i] = mod$p.value
}

out = matrix(NA, ncol = 2, nrow = ncol(Xdata))
colnames(out) = c('Estimate', 'Pr(>|t|)')
for (i in 1:ncol(Xdata)){
  data$x = Xdata[,i]
  mod = glm(log(nodes_pos+1) ~ x, data = data) 
  if ('x' %in% rownames(summary(mod)$coef))  out[i,] = summary(mod)$coef['x',colnames(out)]
}
plot(ps, out[,2])
abline(0,1)

par(mfrow = c(1,2))
hist(ps,breaks=c(0,0.001,0.01,0.05,0.1,0.2,0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9, 1),
     main = 'Nodes Wilcoxon test P-values', xlab = 'p-values')
fdr = p.adjust(ps, method = 'BH')
hist(fdr,breaks=c(0,0.001,0.01,0.05,0.1,0.2,0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9, 1),
     main = 'Nodes Wilcoxon test FDR values', xlab = 'FDR')
#dev.print(png, file=paste(prefix.out,'results/Pvalue2_hist_mutations.png', sep = ''), 
#          width=400, height=400)
qqnorm(out[,1])
#dev.print(png, file=paste(prefix.out,'results/TILcoefficients2_mutations.png', sep = ''), 
#          width=400, height=400)

fdr = p.adjust(out$ps, method = 'BH')

#Toptable
toptable = data.frame('Hugo.Symbol' = colnames(Xdata), 
                      'N samples with mutation' = out$n.with.mut,
                      'N samples without mutation' = out$n.without.mut,
                      'Median nodes with mutation' = out$n.with.mut,
                      'Median nodes without mutation' = out$n.without.mut,
                      p.value = ps,
                      FDR = fdr)
toptable = toptable[order(toptable$p.value),]
toptable$rank = 1:nrow(toptable)
toptable = toptable[, c(ncol(toptable), 1:(ncol(toptable)-1))]
write.table(toptable, file = paste(prefix.out, 'Nodes_by_mutation.txt', sep = ''), 
            row.names = F, quote = F, sep = '\t')

#############################################
# Fisher tests ER ~ mutation

prefix.out = paste(getwd(),'/AROMA/reports/ER_vs_mutations/', sep = '')

out = data.frame(ps = rep(NA, ncol(Xdata)), freq.in.ERneg = NA, freq.in.ERpos = NA)
for (i in 1:ncol(Xdata)){
  data$x = Xdata[,i]
  tab = table(data$x, data$ER_status)
  out$ps[i] = fisher.test(tab)$p.value
  out$freq.in.ERneg[i] = tab[2,1]
  out$freq.in.ERpos[i] = tab[2,2]
}

par(mfrow = c(1,2))
hist(out$ps,breaks=c(0,0.001,0.01,0.05,0.1,0.2,0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9, 1),
     main = 'ER-status Fisher test P-values', xlab = 'p-values')
fdr = p.adjust(out$ps, method = 'BH')
hist(fdr,breaks=c(0,0.001,0.01,0.05,0.1,0.2,0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9, 1),
     main = 'ER-status Fisher test FDR values', xlab = 'FDR')
#dev.print(png, file=paste(prefix.out,'results/Pvalue2_hist_mutations.png', sep = ''), 
#          width=400, height=400)
qqnorm(out$ps)
#dev.print(png, file=paste(prefix.out,'results/TILcoefficients2_mutations.png', sep = ''), 
#          width=400, height=400)

paste(out$freq.in.ERneg ,' (',round(out$freq.in.ERneg/table(data$ER_status)['0']*100), ' %)', sep = '')

#Toptable
toptable = data.frame('Hugo.Symbol' = colnames(Xdata), 
                      Frequency.in.ER.neg.samples = paste(out$freq.in.ERneg ,' (',
                            round(out$freq.in.ERneg/table(data$ER_status)['0']*100), ' %)', 
                                                               sep = ''),
                      Frequency.in.ER.pos.samples = paste(out$freq.in.ERpos ,' (',
                            round(out$freq.in.ERpos/table(data$ER_status)['1']*100), ' %)', 
                                                               sep = ''),
                      p.value = out$ps,
                      FDR = fdr)

toptable = toptable[order(toptable$p.value),]
toptable$rank = 1:nrow(toptable)
toptable = toptable[, c(ncol(toptable), 1:(ncol(toptable)-1))]
write.table(toptable, file = paste(prefix.out, 'ER_by_mutation.txt', sep = ''), 
            row.names = F, quote = F, sep = '\t')




################################################################
################################################################
#
# Expression data and mutations: data preparation
#
################################################################
################################################################

#First save mutation data as Mdata, to avoid name clash
Mdata = Xdata

################################################################
### Organize data

#Load expression data
load(paste(prefix.norm, 'Expression_estimates.RData', sep = ''))
Xdata = as.data.frame(t(as.matrix(expr[,6:ncol(expr)])))
names(Xdata) = paste('unit', expr$unit, sep = '')
unitNames = expr$unitName
rm(expr)

#Load gene annotations
mysymbols = c('CTLA4','CD274','PDCD1','CD8A','IFNG','IDO1','APOBEC3B','APOBEC3A', 'AICDA',
            'NT5E', 'ADORA2A', 'ENTPD1', 'ADORA2B')
d = file.path(getwd(), 'AROMA', 'annotationData', 'chipTypes','HG-U133_Plus_2','HG-U133_Plus_2.na33.annot.csv')
ann = read.csv(d, skip = 25)
myann = as.data.frame(subset(ann, Gene.Symbol %in% mysymbols))
myann$Gene.Symbol = as.character(myann$Gene.Symbol)
table(myann$Gene.Symbol)
#ADORA2A    AICDA APOBEC3A APOBEC3B    CD274     CD8A    CTLA4   ENTPD1     IDO1     IFNG     NT5E    PDCD1 
#     1        2        1        1        2        1        5        4        1        1        4        1 


#Extract expression data for selected genes
index = which(unitNames %in% myann$Probe.Set.ID)
myunits = unitNames[index]
Xdata = Xdata[,index]
colnames(Xdata) = myunits
rownames(Xdata) = substr(rownames(Xdata), 3, 6)

#Make sure data, Mdata and Xdata have the same samples
identical(rownames(Mdata), data$ID) 
#[1] TRUE #Should be true from above
noX = setdiff(data$ID, rownames(Xdata))
noX
#[1] "4644" #This sample must be taken away from data and Mdata, because no expression
data = data[data$ID != noX,]
Mdata = Mdata[rownames(Mdata) != noX,]
Xdata = Xdata[rownames(Mdata),]
identical(rownames(Xdata), rownames(Mdata))
#[1] TRUE
identical(rownames(Xdata), data$ID)
#[1] TRUE

#Now, data and Xdata could be cbind() and put into linear models.
save(data,Mdata,Xdata,myunits,myann, file = 
       "/wehisan/home/allstaff/l/lonnstedt/RESPONSIFY/AROMA/reports/Expression_vs_mutation/Expression_vs_mutation.RData")
load(file.path(getwd(),"AROMA/reports/Expression_vs_mutation/Expression_vs_mutation.RData"))


################################################################
################################################################
#
# Expression data and mutations: data preparation
#
################################################################
################################################################
#Uses Xdata, Mdata, data, myann and myunits from previous section


#############################################
# Multiple linear models expression ~ mutation
prefix.out = paste(getwd(),'/AROMA/reports/Expression_vs_mutation/', sep = '')

res = NULL
for (probe in 1:nrow(myann)){
  data$y = Xdata[,probe] 
  out = data.frame(Gene = myann$Gene.Symbol[probe], 
                   Probe = colnames(Xdata)[probe],Mutation = colnames(Mdata),
                   Estimate = NA, P.value = NA)
  for (i in 1:ncol(Mdata)){
    data$x = Mdata[,i]
    mod = lm(log2(y) ~ x + age_diag + tumor_size + Nodal_status + ER_status, data = data)
    if ('x' %in% rownames(summary(mod)$coef))  {
      out[i,'Estimate'] = summary(mod)$coef['x','Estimate']
      out[i,'P.value'] = summary(mod)$coef['x','Pr(>|t|)']
    }
  }
  res = rbind(res, out)
}
res$FDR = p.adjust(res[,'P.value'], method = 'BH')

par(mfrow = c(1,2))
hist(res[,'P.value'],breaks=c(0,0.001,0.01,0.05,0.1,0.2,0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9, 1),
     main = 'Expression ~ mutations lm P-values', xlab = 'p-values')
hist(res[,'FDR'],breaks=c(0,0.001,0.01,0.05,0.1,0.2,0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9, 1),
     main = 'Expression ~ mutations lm FDR', xlab = 'FDR')

#Reshape output res
res[,4] = format(round(res[,4],2), nsmall = 2)
res[,5] = format(round(res[,5],3), nsmall = 3)
res[,6] = format(round(res[,6],3), nsmall = 3)

resw = NULL
for (probe in 1:nrow(myann)){
  top = rbind(c(myann$Gene.Symbol[probe], '', ''), 
              c(colnames(Xdata)[probe],'',''),
              names(res)[4:6])
  tmp = as.matrix(res[((probe-1)*ncol(Mdata) + 1):(probe*ncol(Mdata)),4:6])
  resw = cbind(resw,rbind(top, tmp))
}
resw = cbind(c('GENE->','PROBESET->','Mutation', colnames(Mdata)), resw)

write.table(resw, file = paste(prefix.out, 'Expression_by_mutation.txt', sep = ''), 
            row.names = F, quote = F, sep = '\t')



#############################################
# Multiple linear models

out = matrix(NA, ncol = 2, nrow = ncol(Xdata))
colnames(out) = c('Estimate', 'Pr(>|t|)')
for (i in 1:ncol(Xdata)){
  data$x = log2(Xdata[,i])
  mod = lm(til ~ x + batch + age_diag + tumor_size + Nodal_status + ER_status, data = data)
  out[i,] = summary(mod)$coef['x',colnames(out)]
}
hist(out[,2],breaks=c(0,0.001,0.01,0.05,0.1,0.2,0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9, 1),
     main = 'Expression linear model P-values', xlab = 'p-values')
out = cbind(out, p.adjust(out[,2], method = 'BH'))
hist(out[,3],breaks=c(0,0.001,0.01,0.05,0.1,0.2,0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9, 1),
     main = 'Expression linear model FDR values', xlab = 'FDR')
mtext('Multiple models')
dev.print(png, file=paste(prefix.out,'results/Pvalue2_hist_expression.png', sep = ''), 
          width=400, height=400)
qqnorm(out[,1])
mtext('Multiple TIL ~ expression model coefficients')
dev.print(png, file=paste(prefix.out,'results/TILcoefficients2_expression.png', sep = ''), 
          width=400, height=400)

#Index of the probe sets with p-value<0.01:
indexLM = which(out[, 2]<0.001)   #Gives 1187 probe sets

#############################################
# Rearrangements of Multiple linear models results

#Gene symbol list for GOstat test of enrichment
highann = subset(ann, Probe.Set.ID %in% unitNames[indexLM])
highsymbols = unique(highann$Gene.Symbol)
highsymbols = highsymbols[highsymbols != '---']
write.table(highsymbols, file = paste(prefix.out, 'TIL_on_expr2_symbols.txt'), col.names=F,
            row.names = F, quote = F)

#Toptable
toptable = data.frame('Probe.Set.Number' = (1:nrow(out))[indexLM], 
                      'Probe.Set.ID' = unitNames[indexLM], 
                      log2exp = out[indexLM, 'Estimate'], 
                      p.value = out[indexLM, 'Pr(>|t|)'],
                      FDR = out[indexLM, 3])
highann = subset(ann, Probe.Set.ID %in% unitNames[indexRUV])
tmp = subset(highann, select = c('Probe.Set.ID','Gene.Symbol'))
toptable = merge(tmp, toptable)
toptable = toptable[order(toptable$FDR),]
toptable$rank = 1:nrow(toptable)
toptable = toptable[, c(ncol(toptable), 1:(ncol(toptable)-1))]
write.table(toptable, file = paste(prefix.out, 'TIL_by_exprLM_toptable.txt'), 
            row.names = F, quote = F)

#Heatmap
set.seed(567890)
set = paste('unit', c(sample(ncol(Xdata), size = 10), toptable$Probe.Set.Number[c(50:1)]), sep = '')
mat = Xdata[,set]
mat = t(mat)
mat = mat[, order(data$til)]
heatmap(log2(mat), Rowv = NA, Colv = NA, col = heat.colors(256), margin = c(10,7))



save.image(paste(prefix.out, 'TILs_on_expr.RData', sep = '')) 


#############################################
# Rearrangements of Naive RUV results, which will be considered primary
(the 835 probes that intersect with the linear models will be primary)

topindex = intersect(indexLM, indexNR)

#Gene symbol list for GOstat test of enrichment
highann = subset(ann, Probe.Set.ID %in% unitNames[topindex])
highsymbols = unique(highann$Gene.Symbol)
highsymbols = highsymbols[highsymbols != '---']
write.table(highsymbols, file = paste(prefix.out, 'TIL_on_expr_TOP_symbols.txt'), col.names=F,
            row.names = F, quote = F, sep = '\t')

#Toptable
toptable = data.frame('Probe.Set.Number' = (1:nrow(res))[topindex], 
                      'Probe.Set.ID' = unitNames[topindex], 
                      log2exp = res[topindex, 'Estimate'], 
                      p.value = res[topindex, 'Pr(>|t|)'],
                      FDR = p.adjust(res[,2], method = 'fdr')[topindex])
highann = subset(ann, Probe.Set.ID %in% unitNames[topindex])
tmp = subset(highann, select = c('Probe.Set.ID','Gene.Symbol'))
toptable = merge(tmp, toptable)
toptable = toptable[order(toptable$FDR),]
toptable$rank = 1:nrow(toptable)
toptable = toptable[, c(ncol(toptable), 1:(ncol(toptable)-1))]
write.table(toptable, file = paste(prefix.out, 'TIL_by_expression_toptable.txt'), 
            row.names = F, quote = F, sep = '\t')

#Heatmap
set.seed(12345)
set = paste('unit', c(sample(ncol(Xdata), size = 10), toptable$Probe.Set.Number[c(50:1)]), sep = '')
mat = Xdata[,set]
mat = t(mat)
mat = mat[, order(data$til)]
heatmap(log2(mat), Rowv = NA, Colv = NA, col = heat.colors(256), margin = c(10,7))

library(RColorBrewer)
hmcols<-colorRampPalette(c('blue','white','orange','red', 'black'))(256)

heatmap(log2(mat), Rowv = NA, Colv = NA, col = hmcols, margin = c(5,0), labRow = '', labCol = '',
        RowSideColors = c(rep('darkorange', 10), rep('darkblue', 50)))
#axis(1, line = 2, labels = c('', '','','',''), at = seq(0, 800, 200))
title(xlab = 'Samples ordered by TIL (maximum TIL to the right)', line = 2)
#text(x = 1, y = 10, labels = 'Random probe sets')
title(ylab = '10 random probe sets                                                                          ', line = 0, cex = .8)
title(ylab = '                                          Top 50 probe sets ordered by p-value', line = 0, cex = .8)
title(ylab = '                                              (minimum p at the top)', line = -1, cex = .8)
legend("right", fill = c('black','red','orange','white', 'blue'),
       legend = rep('', 5), bty = 'n', border = 'black', y.intersp = .5, cex = 2)
text(9.5,.4,'Max expression', cex = .8)
text(9.5,-.4,'Min expression', cex = .8)

#dev.print(png, file=paste(prefix.out,'results/TIL_expression_heatmap.png', sep = ''), 
#          width=640, height=500)


save.image(paste(prefix.out, 'TILs_on_expr.RData', sep = '')) 



