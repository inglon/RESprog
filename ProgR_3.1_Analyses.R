################################################################
################################################################
# Script: ProgR_3.1_Analyses.R
# Author: Ingrid Lonnstedt
# Date:  15/03/2013
# R version: R 2.15.2
# Details: Statistical Analyses with clinical data
################################################################
################################################################



################################################################
################################################################
#
# File paths, date
#
################################################################
################################################################
setwd(paste(getwd(), '/RESPONSIFY', sep=''))
#setwd('/Users/lonnstedt/Documents/RESPONSIFY')

prefix.norm = paste(getwd(), "/output/", sep='')
prefix.raw = paste(getwd(), "/rawData/responsify/", sep='')
prefix.ann = paste(getwd(), "/annotationData/chipTypes/", sep='')

date = format(Sys.Date())

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



################################################################
################################################################
#
# Expression data and tils: checks and data preparation
#
################################################################
################################################################
prefix.out = paste(getwd(), "/reports/report1/TILs_on_expression/", sep='')

################################################################
### Organize data

load(paste(prefix.norm, 'Expression/Expression_estimates.RData', sep = ''))

data = data.frame(array = names(expr)[6:ncol(expr)])
data$ID = substr(data$array, 3, 6)
data$index = 1:nrow(data)
data = merge(data, subset(dat, select = c('ID','age','tumor_size','Nodal_status','ER_status', 'til', 'date_diagnosis')), 
             all.x = T, sort = F)
data$batch = 1
data$array = as.character(data$array)
tmp = nchar(data$array, type = 'c')
data$batch[substr(data$array, tmp, tmp) == '2'] = 2
data$batch[substr(data$array, 7, 9) == 'bis'] = 3
data$batch = as.factor(data$batch)
data$age_diag = as.numeric(data$age/365.25)
data$year_diag = as.numeric(substr(as.character(data$date_diagnosis), 1, 4))
data = data[order(data$index),]

###Load array dates
tmp = read.table(paste(prefix.norm, 'Expression/Array_dates.txt', sep = ''), sep = '\t')
#All arrays are hybridized within 30 minutes, so that is not interesting

Xdata = as.data.frame(t(as.matrix(expr[,6:ncol(expr)])))
names(Xdata) = paste('unit', expr$unit, sep = '')

unitNames = expr$unitName
rm(expr)

#Finally remove samples with missing til:
index = which(!is.na(data$til))
data = data[index,]
Xdata = Xdata[index,]

#Now, data and Xdata could be cbind() and put into linear models.

################################################################
### Clinical data overview
png(paste(prefix.out,'results/Clinical_and_exprbatch_overview.png', sep = ''), 
    width=400, height=400)
plot(data[,c('tumor_size', 'til', 'age_diag')])
dev.off()

par(mfrow = c(1,2))
boxplot(til~Nodal_status, data = data, xlab = 'Nodal status', ylab = 'Stomal LI (%)', col = 'wheat')
boxplot(til~ER_status, data = data, xlab = 'ER status', ylab = 'Stomal LI (%)', col = 'wheat',
        xaxt = 'n')
axis(1, labels = c('ER-','ER+'), at = 1:2)
dev.print(png, file=paste(prefix.out,'results/Til_by_nodal_and_ER.png', sep = ''), 
          width=400, height=400)

summary(data[,c(5, 8, 11)])

################################################################
### Check unwanted variation through principal components plots

pc = prcomp(Xdata)
## Get the standard deviation and variance per principal component
sd.per.pc <- pc$sdev
var.per.pc <- sd.per.pc^2
## Display the percentage of total variance explained by each
sd.per.pc.percent <- sd.per.pc/sum(sd.per.pc)
var.per.pc.percent <- var.per.pc/sum(var.per.pc)


png(paste(prefix.out,'results/Expression perc variance per pc.png', sep = ''), 
    width=400, height=400)
barplot(var.per.pc.percent[1:10]*100, main='Gene expression, Percent of variance per component', 
        xlab='Component', ylab = '% variance')
dev.off()


newids = data$ID
#newids[c(19, 78)] = paste(newids[c(19, 78)], 'bis', sep = '')
#newids[c(25, 26, 46, 51, 80)] = paste(newids[c(25, 26, 46, 51, 80)], '_2', sep = '')
dimnames(pc$x)[[1]] = newids



## Plot components PC1 and PC2
plot(pc$x[,1:2],
     type='n',
     panel.first=grid(col='black'),
     main=paste('PCA gene expression ',
                nrow(Xdata), 'samples *', ncol(Xdata), 'genes', sep=' '),
     xlab='PC1', ylab='PC2')
text(pc$x[,1:2],labels=dimnames(pc$x)[[1]],col=as.numeric(data$batch),pch=0.5)
legend('bottomleft',col=as.numeric(unique(data$batch)),
       legend=paste('Batch', as.numeric(unique(data$batch))),pch=1,cex=0.7,bg='white',bty='o')


## Plot components PC1 and PC3
plot(pc$x[,c(1,3)],
     type='n',
     panel.first=grid(col='black'),
     main=paste('PCA gene expression ',
                nrow(Xdata), 'samples *', ncol(Xdata), 'genes', sep=' '),
     xlab='PC1', ylab='PC3')
text(pc$x[,c(1,3)],labels=dimnames(pc$x)[[1]],col=as.numeric(data$batch),pch=0.5)
legend('bottomleft',col=as.numeric(unique(data$batch)),
       legend=paste('Batch', as.numeric(unique(data$batch))),pch=1,cex=0.7,bg='white',bty='o')


## Plot components PC1 and PC4
plot(pc$x[,c(1,4)],
     type='n',
     panel.first=grid(col='black'),
     main=paste('PCA gene expression ',
                nrow(Xdata), 'samples *', ncol(Xdata), 'genes', sep=' '),
     xlab='PC1', ylab='PC4')
text(pc$x[,c(1,4)],labels=dimnames(pc$x)[[1]],col=as.numeric(data$batch),pch=0.5)
legend('bottomleft',col=as.numeric(unique(data$batch)),
       legend=paste('Batch', as.numeric(unique(data$batch))),pch=1,cex=0.7,bg='white',bty='o')

## Plot components PC2 and PC3
plot(pc$x[,c(2,3)],
     type='n',
     panel.first=grid(col='black'),
     main=paste('PCA gene expression ',
                nrow(Xdata), 'samples *', ncol(Xdata), 'genes', sep=' '),
     xlab='PC2', ylab='PC3')
text(pc$x[,c(2,3)],labels=dimnames(pc$x)[[1]],col=as.numeric(data$batch),pch=0.5)
legend('bottomleft',col=as.numeric(unique(data$batch)),
       legend=paste('Batch', as.numeric(unique(data$batch))),pch=1,cex=0.7,bg='white',bty='o')


## Plot components PC2 and PC4
plot(pc$x[,c(2,4)],
     type='n',
     panel.first=grid(col='black'),
     main=paste('PCA gene expression ',
                nrow(Xdata), 'samples *', ncol(Xdata), 'genes', sep=' '),
     xlab='PC2', ylab='PC4')
text(pc$x[,c(2,4)],labels=dimnames(pc$x)[[1]],col=as.numeric(data$batch),pch=0.5)
legend('bottomleft',col=as.numeric(unique(data$batch)),
       legend=paste('Batch', as.numeric(unique(data$batch))),pch=1,cex=0.7,bg='white',bty='o')

## Plot components PC3 and PC4
plot(pc$x[,c(3,4)],
     type='n',
     panel.first=grid(col='black'),
     main=paste('PCA gene expression ',
                nrow(Xdata), 'samples *', ncol(Xdata), 'genes', sep=' '),
     xlab='PC3', ylab='PC4')
text(pc$x[,c(3,4)],labels=dimnames(pc$x)[[1]],col=as.numeric(data$batch),pch=0.5)
legend('bottomleft',col=as.numeric(unique(data$batch)),
       legend=paste('Batch', as.numeric(unique(data$batch))),pch=1,cex=0.7,bg='white',bty='o')



## Plot components PC1 and PC2 coloured by ER status
plot(pc$x[,1:2],
     type='n',
     panel.first=grid(col='black'),
     main=paste('PCA gene expression ',
                nrow(Xdata), 'samples *', ncol(Xdata), 'genes', sep=' '),
     xlab='PC1', ylab='PC2')
text(pc$x[,1:2],labels=dimnames(pc$x)[[1]],col=as.numeric(data$ER_status)+1,pch=0.5)
legend('bottomleft',col=as.numeric(unique(data$ER_status)+1),
       legend=paste('ER', as.numeric(unique(data$ER_status))),pch=1,cex=0.7,bg='white',bty='o')



## Plot components PC1 and PC2 coloured by Nodal status
plot(pc$x[,1:2],
     type='n',
     panel.first=grid(col='black'),
     main=paste('PCA gene expression ',
                nrow(Xdata), 'samples *', ncol(Xdata), 'genes', sep=' '),
     xlab='PC1', ylab='PC2')
text(pc$x[,1:2],labels=dimnames(pc$x)[[1]],col=as.numeric(data$Nodal_status)+1,pch=0.5)
legend('bottomleft',col=as.numeric(unique(data$Nodal_status)+1),
       legend=paste('Nodal status', as.numeric(unique(data$Nodal_status))),pch=1,cex=0.7,bg='white',bty='o')

## Plot components PC1 and PC2 coloured by year of diagnosis
plot(pc$x[,1:2],
     type='n',
     panel.first=grid(col='black'),
     main=paste('PCA gene expression ',
                nrow(Xdata), 'samples *', ncol(Xdata), 'genes', sep=' '),
     xlab='PC1', ylab='PC2')
text(pc$x[,1:2],labels=dimnames(pc$x)[[1]],col=data$year_diag,pch=0.5)
legend('bottomleft',col=as.numeric(unique(data$year_diag)),
       legend=paste('Year of diagnosis', as.numeric(unique(data$year_diag))),pch=1,cex=0.7,bg='white',bty='o')

## Plot components PC1 and PC2 coloured by age_diag
colvec = 1
colvec[data$age_diag>55] = 2
plot(pc$x[,1:2],
     type='n',
     panel.first=grid(col='black'),
     main=paste('PCA gene expression ',
                nrow(Xdata), 'samples *', ncol(Xdata), 'genes', sep=' '),
     xlab='PC1', ylab='PC2')
text(pc$x[,1:2],labels=dimnames(pc$x)[[1]],col=unique(colvec),pch=0.5)
legend('bottomleft',col=unique(colvec),
       legend=paste('Age at diag', unique(colvec)),pch=1,cex=0.7,bg='white',bty='o')

## Plot components PC1 and PC2 coloured by tumor_size
colvec = 1
colvec[data$tumor_size>40] = 2
plot(pc$x[,1:2],
     type='n',
     panel.first=grid(col='black'),
     main=paste('PCA gene expression ',
                nrow(Xdata), 'samples *', ncol(Xdata), 'genes', sep=' '),
     xlab='PC1', ylab='PC2')
text(pc$x[,1:2],labels=dimnames(pc$x)[[1]],col=unique(colvec),pch=0.5)
legend('bottomleft',col=unique(colvec),
       legend=paste('Tumour size', unique(colvec)),pch=1,cex=0.7,bg='white',bty='o')


data$pc1 = pc$x[,1]
plot(data[,5:ncol(data)])

tmp = subset(data, select = c(ID, pc1))
dat = merge(dat, tmp, all.x = T)
par(ask = T)
for (i in 70:(ncol(dat)-1)){ plot(dat$pc1,dat[,i], ylab = names(dat)[i]) }

#We have dependendy between PC1 and some biology variables, but not with array dates or similar.





################################################################
################################################################
#
# Genes expressions assoxiated with TIL
#
################################################################
################################################################
prefix.out = paste(getwd(), "/reports/report1/TILs_on_expression/", sep='')

#############################################
# RUV model

library(ruv)

rinv<-RUVinv(Y=as.matrix(log2(Xdata)), X=as.matrix(data$til, ncol = 1), ctl = rep(T, ncol(Xdata)), Z = NULL)
#hist(rinv$p,breaks=c(0,0.001,0.01,0.05,0.1,0.2,0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9, 1))
radj<-variance_adjust(rinv)
#hist(radj$p.rsvar,breaks=c(0,0.001,0.01,0.05,0.1,0.2,0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9, 1))
#hist(radj$p.ebayes,breaks=c(0,0.001,0.01,0.05,0.1,0.2,0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9, 1))
hst = hist(radj$p.evar,breaks=c(0,0.001,0.01,0.05,0.1,0.2,0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9, 1),
           main = 'Empirical variances expression p-values', xlab = 'p-values')
dev.print(png, file=paste(prefix.out,'results/Pvalue_hist_expression.png', sep = ''), 
          width=400, height=400)
qqnorm(rinv$betahat)
mtext('Expression ~ TIL coefficients')
dev.print(png, file=paste(prefix.out,'results/TILcoefficients_expression.png', sep = ''), 
          width=400, height=400)


#Index of the probe sets with p-value<0.01:
indexRUV = which(radj$p.evar<0.01)   #Gives 785 probe sets

ann = read.csv(paste(prefix.ann, 'HG-U133_Plus_2/HG-U133_Plus_2.na33.annot.csv',
                     sep = ''), skip = 25)


#radj has p-values (p.evar)
#indexRUV has probe sets with low p-values
#ann has annotation data for the complete array
#hst has p-value histogram

#############################################
# Rearrangements of RUV model results

#Gene symbol list for GOstat test of enrichment
highann = subset(ann, Probe.Set.ID %in% unitNames[indexRUV])
highsymbols = unique(highann$Gene.Symbol)
highsymbols = highsymbols[highsymbols != '---']
write.table(highsymbols, file = paste(prefix.out, 'Expr_on_TIL_symbols.txt'), col.names=F,
            row.names = F, quote = F)

#FDR estimates and toptable
mean.freq = mean((hst$counts/(hst$breaks[-1] - hst$breaks[-length(hst$breaks)]))[3:length(hst$counts)])
toptable = data.frame('Probe.Set.Number' = (1:length(radj$betahat))[indexRUV], 'Probe.Set.ID' = unitNames[indexRUV], 
                      log2exp = radj$betahat[indexRUV], p.value = radj$p.evar[indexRUV])
Nle = numeric(nrow(toptable))
Nle[order(toptable$p.value)] = 1:nrow(toptable)
Ele = toptable$p.value * mean.freq
toptable$FDR = Ele/(Nle)

highann = subset(ann, Probe.Set.ID %in% unitNames[indexRUV])
tmp = subset(highann, select = c('Probe.Set.ID','Gene.Symbol'))
toptable = merge(tmp, toptable)
toptable = toptable[order(toptable$FDR),]
toptable$rank = 1:nrow(toptable)
toptable = toptable[, c(ncol(toptable), 1:(ncol(toptable)-1))]
write.table(toptable, file = paste(prefix.out, 'Expr_by_TIL_toptable.txt'), 
            row.names = F, quote = F)

#Heatmap
set = paste('unit', c(sample(ncol(Xdata), size = 10), toptable$Probe.Set.Number[c(50:1)]), sep = '')
mat = Xdata[,set]
mat = t(mat)
mat = mat[, order(data$til)]
heatmap(log2(mat), Rowv = NA, Colv = NA, col = heat.colors(256), margin = c(10,7))
abline(h = 11.5)
abline(v = 53.4)

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


################################################################
# Naive RUV adjustment model

## Y: expression matrix where the rows are the samples and the columns are the genes.
## cIdx: column index of the negative control genes in Y, for estimation of unwanted variation.
## nuCoeff: regularization parameter for the unwanted variation.
## k: rank of the unwanted variation term.


naiveRandRuv <- function(Y, cIdx, nuCoeff=1e-3, k=m){
  
  ## W is the square root of the empirical covariance on the control
  ## genes.
  
  svdYc <- svd(Y[, cIdx])
  W <- svdYc$u[, 1:k] %*% diag(svdYc$d[1:k]) #/ sqrt(length(cIdx)+1)
  
  ## Regularization heuristic: nu is a fraction of the largest eigenvalue of WW'
  
  nu <- nuCoeff*svdYc$d[1]^2 #/ (length(cIdx)+1)
  
  ## Naive correction: ridge regression of Y against W
  
  nY <- Y - W %*% solve(t(W)%*%W + nu*diag(k), t(W) %*% Y)
  
  return(nY)
}



rnaive10 = naiveRandRuv(Y=as.matrix(log2(Xdata)), cIdx = rep(T, ncol(Xdata)), nuCoeff=1e-3, k=10)
res10 = matrix(NA, ncol = 2, nrow = ncol(Xdata))
colnames(res10) = c('Estimate', 'Pr(>|t|)')
for (i in 1:ncol(Xdata)){
  data$x = rnaive10[,i]
  mod = lm(til ~ x + age_diag + tumor_size + Nodal_status + ER_status, data = data)
  res10[i,] = summary(mod)$coef['x',colnames(res10)]
}
hist(res10[,2],breaks=c(0,0.001,0.01,0.05,0.1,0.2,0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9, 1),
     main = 'Naive RUV P-values', xlab = 'p-values')
mtext('k = 10')

rnaive = naiveRandRuv(Y=as.matrix(log2(Xdata)), cIdx = rep(T, ncol(Xdata)), nuCoeff=1e-3, k=20)
res = matrix(NA, ncol = 2, nrow = ncol(Xdata))
colnames(res) = c('Estimate', 'Pr(>|t|)')
for (i in 1:ncol(Xdata)){
  data$x = rnaive[,i]
  mod = lm(til ~ x + age_diag + tumor_size + Nodal_status + ER_status, data = data)
  res[i,] = summary(mod)$coef['x',colnames(res)]
}
hist(res[,2],breaks=c(0,0.001,0.01,0.05,0.1,0.2,0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9, 1),
     main = 'Naive RUV P-values', xlab = 'p-values')
mtext('k = 20')


rnaive100 = naiveRandRuv(Y=as.matrix(log2(Xdata)), cIdx = rep(T, ncol(Xdata)), nuCoeff=1e-1, k=20)
res100 = matrix(NA, ncol = 2, nrow = ncol(Xdata))
colnames(res100) = c('Estimate', 'Pr(>|t|)')
for (i in 1:ncol(Xdata)){
  data$x = rnaive100[,i]
  mod = lm(til ~ x + age_diag + tumor_size + Nodal_status + ER_status, data = data)
  res100[i,] = summary(mod)$coef['x',colnames(res100)]
}
hist(res100[,2],breaks=c(0,0.001,0.01,0.05,0.1,0.2,0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9, 1),
     main = 'Naive RUV P-values', xlab = 'p-values')
mtext('Times 100')

rnaive.01 = naiveRandRuv(Y=as.matrix(log2(Xdata)), cIdx = rep(T, ncol(Xdata)), nuCoeff=1e-5, k=20)
res.01 = matrix(NA, ncol = 2, nrow = ncol(Xdata))
colnames(res.01) = c('Estimate', 'Pr(>|t|)')
for (i in 1:ncol(Xdata)){
  data$x = rnaive.01[,i]
  mod = lm(til ~ x + age_diag + tumor_size + Nodal_status + ER_status, data = data)
  res.01[i,] = summary(mod)$coef['x',colnames(res.01)]
}
hist(res.01[,2],breaks=c(0,0.001,0.01,0.05,0.1,0.2,0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9, 1),
     main = 'Naive RUV P-values', xlab = 'p-values')
mtext('Div 100')

#Index of the probe sets with p-value<0.01:
indexNR10 = which(res10[, 2]<0.001)   #Gives 862 probe sets, out of which 859 are shared with indexNR
indexNR = which(res[, 2]<0.001)   #Gives 889 probe sets, out of which 835 are shared with indexLM
indexNR100 = which(res100[, 2]<0.001)   #Gives 1231 probe sets, out of which 877 are shared with indexNR
#and 1021 are shared with indexLM
indexNR.01 = which(res.01[, 2]<0.001)   #Gives 26 probe sets, this histogram looks bad
indexNRx = which(res[, 2]<0.01)   #Gives 889 probe sets, out of which 835 are shared with indexNR


#Decided to use the intersect(indexNR, indexLM) as toptable
#Print p-value histogram of the naive RUV with k = 20 and nu = 1e-3 (indexNR model)
hstNR = hist(res[,2],breaks=c(0,0.001,0.01,0.05,0.1,0.2,0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9, 1),
             main = 'Expression naive RUV p-values', xlab = 'p-values')
mtext('Multiple models')
dev.print(png, file=paste(prefix.out,'results/Pvalue3_hist_expression.png', sep = ''), 
          width=400, height=400)


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







#####################################################################
#####################################################################
#
# IDFS relapse associated with gene expression
#
#
#####################################################################
#####################################################################
prefix.out = paste(getwd(), "/reports/report1/Relapse_on_expression/", sep='')

### Read IDFS relapse data and add to the other clinical data in dat
tmp = read.delim(paste(prefix.raw, 'HER2-lev-bordet-Aug2012-responsify.txt',
                       sep = ''), na.strings = c('', NA), check.names = F, dec = ',')
tmp$ID = as.character(dat$Frozen_tissue_BO_no)
tmp = subset(tmp, select = c('ID','IDFS_Y_N','IDFS_days'))
dat = merge(dat, tmp, sort = F, all.x = T)


################################################################
### Organize data

load(paste(prefix.norm, 'Expression/Expression_estimates.RData', sep = ''))

data = data.frame(array = names(expr)[6:ncol(expr)])
data$ID = substr(data$array, 3, 6)
data$index = 1:nrow(data)
data = merge(data, subset(dat, select = c('ID','age','tumor_size','Nodal_status',
                                          'ER_status', 'til', 'date_diagnosis',
                                          'IDFS_Y_N','IDFS_days')), 
             all.x = T, sort = F)
data$batch = 1
data$array = as.character(data$array)
tmp = nchar(data$array, type = 'c')
data$batch[substr(data$array, tmp, tmp) == '2'] = 2
data$batch[substr(data$array, 7, 9) == 'bis'] = 3
data$batch = as.factor(data$batch)
data$age_diag = as.numeric(data$age/365.25)
data$year_diag = as.numeric(substr(as.character(data$date_diagnosis), 1, 4))
data = data[order(data$index),]
data$IDFS_days = as.numeric(data$IDFS_days)

Xdata = as.data.frame(t(as.matrix(expr[,6:ncol(expr)])))
names(Xdata) = paste('unit', expr$unit, sep = '')

unitNames = expr$unitName
rm(expr)
gc()

#Finally remove samples with missing IDFS:

index = which(!is.na(data$IDFS_Y_N) & !is.na(data$IDFS_days))
data = data[index,]
Xdata = Xdata[index,]
#82 samples

#Now, data and Xdata could be cbind() and put into Cox models.


############################################################
# Cox models without naive RUV adjustment
library(survival)
out = matrix(NA, ncol = 2, nrow = ncol(Xdata))
colnames(out) = c('exp(coef)', 'Pr(>|z|)')
for (i in 1:ncol(Xdata)){
  data$x = log2(Xdata[,i])
  mod = coxph(Surv(IDFS_days, IDFS_Y_N) ~ x + age_diag + tumor_size + Nodal_status + ER_status, data = data)
  out[i,] = summary(mod)$coef['x',colnames(out)]
}
hist(out[,2],breaks=c(0,0.001,0.01,0.05,0.1,0.2,0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9, 1),
     main = 'Cox model P-values', xlab = 'p-values')
mtext('Multiple models')
dev.print(png, file=paste(prefix.out,'Pvalue_hist_relapse_expression.png', sep = ''), 
          width=400, height=400)
out = cbind(out, p.adjust(out[,2], method = 'BH'))
hist(out[,3],breaks=c(0,0.001,0.01,0.05,0.1,0.2,0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9, 1),
     main = 'Cox model FDR values', xlab = 'FDR')
qqnorm(out[,1])
mtext('Multiple Relapse time ~ expression model coefficients')

#Index of the probe sets with p-value<0.01:
indexLM = which(out[, 2]<0.01)   #Gives 352 probe sets
length(indexLM)


################################################################
# Naive RUV adjustment (Try 4 different levels of adjustment)

rnaive10 = naiveRandRuv(Y=as.matrix(log2(Xdata)), cIdx = rep(T, ncol(Xdata)), nuCoeff=1e-3, k=10)
res10 = matrix(NA, ncol = 2, nrow = ncol(Xdata))
colnames(res10) = c('exp(coef)', 'Pr(>|z|)')
for (i in 1:ncol(Xdata)){
  data$x = rnaive10[,i]
  mod = coxph(Surv(IDFS_days, IDFS_Y_N) ~ x + age_diag + tumor_size + Nodal_status + ER_status, data = data)
  res10[i,] = summary(mod)$coef['x',colnames(res10)]
}
hist(res10[,2],breaks=c(0,0.001,0.01,0.05,0.1,0.2,0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9, 1),
     main = 'Naive RUV P-values', xlab = 'p-values')
mtext('k = 10')

rnaive = naiveRandRuv(Y=as.matrix(log2(Xdata)), cIdx = rep(T, ncol(Xdata)), nuCoeff=1e-3, k=20)
res = matrix(NA, ncol = 2, nrow = ncol(Xdata))
colnames(res) = c('exp(coef)', 'Pr(>|z|)')
for (i in 1:ncol(Xdata)){
  data$x = rnaive[,i]
  mod = coxph(Surv(IDFS_days, IDFS_Y_N) ~ x + age_diag + tumor_size + Nodal_status + ER_status, data = data)
  res[i,] = summary(mod)$coef['x',colnames(res)]
}
hist(res[,2],breaks=c(0,0.001,0.01,0.05,0.1,0.2,0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9, 1),
     main = 'Naive RUV P-values', xlab = 'p-values')
mtext('k = 20')


rnaive100 = naiveRandRuv(Y=as.matrix(log2(Xdata)), cIdx = rep(T, ncol(Xdata)), nuCoeff=1e-1, k=20)
res100 = matrix(NA, ncol = 2, nrow = ncol(Xdata))
colnames(res100) = c('exp(coef)', 'Pr(>|z|)')
for (i in 1:ncol(Xdata)){
  data$x = rnaive100[,i]
  mod = coxph(Surv(IDFS_days, IDFS_Y_N) ~ x + age_diag + tumor_size + Nodal_status + ER_status, data = data)
  res100[i,] = summary(mod)$coef['x',colnames(res100)]
}
hist(res100[,2],breaks=c(0,0.001,0.01,0.05,0.1,0.2,0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9, 1),
     main = 'Naive RUV P-values', xlab = 'p-values')
mtext('Times 100')

rnaive.01 = naiveRandRuv(Y=as.matrix(log2(Xdata)), cIdx = rep(T, ncol(Xdata)), nuCoeff=1e-5, k=20)
res.01 = matrix(NA, ncol = 2, nrow = ncol(Xdata))
colnames(res.01) = c('exp(coef)', 'Pr(>|z|)')
for (i in 1:ncol(Xdata)){
  data$x = rnaive.01[,i]
  mod = coxph(Surv(IDFS_days, IDFS_Y_N) ~ x + age_diag + tumor_size + Nodal_status + ER_status, data = data)
  res.01[i,] = summary(mod)$coef['x',colnames(res.01)]
}
hist(res.01[,2],breaks=c(0,0.001,0.01,0.05,0.1,0.2,0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9, 1),
     main = 'Naive RUV P-values', xlab = 'p-values')
#mtext('Div 100')
mtext('Multiple models, k = 20, nu = 1e-5')
dev.print(png, file=paste(prefix.out,'Pvalue2_hist_relapse_expression.png', sep = ''), 
          width=400, height=400)

#No significances!!! Therefore we save no further results.
save.image(paste(prefix.out, 'IDFSs_on_expr.RData', sep = '')) 







################################################################
################################################################
#
# CN segments and tils: checks and data preparation
#
################################################################
################################################################
prefix.out = paste(getwd(), "/reports/report1/TILs_on_CN/", sep='')


library(aroma.affymetrix)
library(PSCBS)
library(IRanges)
library(parallel)

################################################################
### Fileinfo: which sample IDs correspond to which arrays?
### Resulting 'info' will be used further down

### Array names
###
filenames = list.files(paste(prefix.raw, 'GenomeWideSNP_6', sep = ''))
filenames = substr(filenames, 1, nchar(filenames, type = 'c')-4)
info = data.frame(arrname = filenames)
info$AROS.ID = substr(info$arrname, 1, 8)
tmp = read.delim(paste(prefix.ann, 'GenomeWideSNP_6/A1761_work sheet_ SNP 6.0_okay_Scan_info.txt', sep = ''),
                 sep = '\t')
info = merge(info, tmp, all = T)[,c('AROS.ID','arrname','Sample.ID')]
info$Sample.ID[info$AROS.ID %in% c('A1974-08','A1974-09','A1974-10')] = 
  c('3860 T','4541 T','4868 T')
info$pass = T
info$pass[info$arrname %in% c('A1761-02','A1761-08','A1761-08rh',
                              'A1761-20','A1761-20rh','A1761-26',
                              'A1761-30','A1761-30rh','A1761-50')] = F

################################################################
### Organize CN data

info$ntcn = T
info$ntcn[info$pass == F] = NA
segs = NULL
for (arr in 1:56){
  #Prepare data for this sample
  samp = as.character((info$AROS.ID[info$pass])[arr])
  fitname = paste('fit', arr, sep = '')
  filename <- paste(prefix.norm, 'SNP/CBSfits/', fitname, '.RData', sep = '');
  fit <- loadObject(filename);
  segtmp = getSegments(fit)
  segtmp$AROS.ID = samp
  if (is.null(segtmp$ntcnCall)) {
    segtmp$ntcnCall = NA
    info$ntcn[info$AROS.ID == samp] = F
  }
  segs = rbind(segs, segtmp)
}

tmp = subset(info, pass, select = c('AROS.ID','Sample.ID'))
tmp$ID = substr(as.character(tmp$Sample.ID), 1, 4)
tmp = tmp[,c('AROS.ID','ID')]
segs = merge(segs, tmp, all.x = T, sort = F) 
segs = subset(segs, !is.na(chromosome))
#segs = subset(segs, ! (AROS.ID %in% c('A1761-20rh','A1761-25','A1761-30rh','A1761-42'))) #Preliminary: exclude difficult samples
segs = subset(segs, chromosome <= 23 & !is.na(tcnMean))
ir = disjoin(IRanges(start = segs$tcnStart, end = segs$tcnEnd-1))

ir = mclapply(1:23, mc.cores = 2, function(k){
  index = which(segs$chromosome == k)
  disjoin(IRanges(start = segs$tcnStart[index], end = segs$tcnEnd[index]-1))
})
names(ir) = as.character(1:23)
ir = do.call('RangesList', ir) #ir defines the segments which will be correlated to TILs
seginfo = as.data.frame(ir)    #seginfo has the ir information in dataframe format
seginfo$mid = (seginfo$start + seginfo$end)/2

Xdata = matrix(NA, ncol = sum(mapply(length, ir)), nrow = length(unique(segs$ID)))
colnames(Xdata) = paste('seg', 1:ncol(Xdata), sep = '')
rownames(Xdata) = unique(segs$ID)
for (i in 1:nrow(Xdata)){
  index = which(segs$ID == unique(segs$ID)[i])
  rd = RangedData(IRanges(start = segs$tcnStart[index], end = segs$tcnEnd[index]-1), tcnMean=segs$tcnMean[index], 
                  space = segs$chromosome[index], universe = 'hg19') #Each interval is redefined to have an open right end
  rle = coverage(rd, weight = 'tcnMean', width = tapply(segs$tcnEnd, segs$chromosome, max))
  Xdata[i, ] = unlist(viewMeans(Views(rle, ir)))
} #Xdata holds segment tcn:s, one row for each sample
#There are 26075 segments

#Number of patients with changes in each segment
Ndata = matrix(NA, ncol = sum(mapply(length, ir)), nrow = length(unique(segs$ID)))
colnames(Ndata) = paste('seg', 1:ncol(Ndata), sep = '')
rownames(Ndata) = unique(segs$ID)
for (i in 1:nrow(Ndata)){
  index = which(segs$ID == unique(segs$ID)[i])
  rd = RangedData(IRanges(start = segs$tcnStart[index], end = segs$tcnEnd[index]-1), ntcnCall=segs$ntcnCall[index], 
                  space = segs$chromosome[index], universe = 'hg19') #Each interval is redefined to have an open right end
  ?
  ? Ndata[i, ] = unlist(viewApply(Views(rle, ir)))
} #Xdata holds segment tcn:s, one row for each sample


################################################################
#Now prepare clinical variables    
tmp = subset(info, pass, select = c('shortname','Sample.ID'))
tmp$ID = substr(as.character(tmp$Sample.ID), 1, 4)
tmp = tmp[,c('shortname','ID')]
names(tmp)[1] = 'array'
data = merge(tmp, subset(dat, select = c('ID','age','tumor_size','Nodal_status','ER_status', 'til', 'date_diagnosis')),
             sort = F, all.x = T)
data$batch = 1
data$batch[nchar(as.character(data$array))>8] = 2
data = subset(data, ID %in% segs$ID) #Limit the arrays to those for which we have CNs
data$index = 1:nrow(data)
data$age_diag = as.numeric(data$age/365.25)
data$year_diag = as.numeric(substr(as.character(data$date_diagnosis), 1, 4))
data = data[order(data$index),]

###Load array dates
#Array dates were checked and were all within 3 hours in the order of info$shortname, not informative

#Finally remove samples with missing til:
index = which(!is.na(data$til))
data = data[index,]
Xdata = Xdata[index,]

#Now, data and Xdata could be cbind() and put into linear models.


################################################################
### Clinical data overview
plot(data[,c('tumor_size', 'til', 'age_diag')]) #No essential difference from the plots of those 
#samples that have gene expression values

par(mfrow = c(1,2))
boxplot(til~Nodal_status, data = data, xlab = 'Nodal status', ylab = 'Stomal LI (%)', col = 'wheat')
boxplot(til~ER_status, data = data, xlab = 'ER status', ylab = 'Stomal LI (%)', col = 'wheat',
        xaxt = 'n')
axis(1, labels = c('ER-','ER+'), at = 1:2)
dev.print(png, file=paste(prefix.out,'results/Til_by_nodal_and_ER_CN.png', sep = ''), 
          width=400, height=400) #Shows that we do not have as many patients with very high TIL as we had with expression data,
#but the ones we have show a higher median TIL in ER+ than ER-.

summary(data[,c(5, 8, 11)])



################################################################
### Check unwanted variation through principal components plots

pc = prcomp(Xdata)
## Get the standard deviation and variance per principal component
sd.per.pc <- pc$sdev
var.per.pc <- sd.per.pc^2
## Display the percentage of total variance explained by each
sd.per.pc.percent <- sd.per.pc/sum(sd.per.pc)
var.per.pc.percent <- var.per.pc/sum(var.per.pc)


png(paste(prefix.out,'results/Expression perc variance per pc.png', sep = ''), 
    width=400, height=400)
barplot(var.per.pc.percent[1:10]*100, main='Gene expression, Percent of variance per component', 
        xlab='Component', ylab = '% variance')
dev.off()


newids = data$ID
#newids[c(19, 78)] = paste(newids[c(19, 78)], 'bis', sep = '')
#newids[c(25, 26, 46, 51, 80)] = paste(newids[c(25, 26, 46, 51, 80)], '_2', sep = '')
dimnames(pc$x)[[1]] = newids



## Plot components PC1 and PC2
plot(pc$x[,1:2],
     type='n',
     panel.first=grid(col='black'),
     main=paste('PCA gene expression ',
                nrow(Xdata), 'samples *', ncol(Xdata), 'genes', sep=' '),
     xlab='PC1', ylab='PC2')
text(pc$x[,1:2],labels=dimnames(pc$x)[[1]],col=as.numeric(data$batch),pch=0.5)
legend('bottomleft',col=as.numeric(unique(data$batch)),
       legend=paste('Batch', as.numeric(unique(data$batch))),pch=1,cex=0.7,bg='white',bty='o')


## Plot components PC1 and PC3
plot(pc$x[,c(1,3)],
     type='n',
     panel.first=grid(col='black'),
     main=paste('PCA gene expression ',
                nrow(Xdata), 'samples *', ncol(Xdata), 'genes', sep=' '),
     xlab='PC1', ylab='PC3')
text(pc$x[,c(1,3)],labels=dimnames(pc$x)[[1]],col=as.numeric(data$batch),pch=0.5)
legend('bottomleft',col=as.numeric(unique(data$batch)),
       legend=paste('Batch', as.numeric(unique(data$batch))),pch=1,cex=0.7,bg='white',bty='o')


## Plot components PC1 and PC4
plot(pc$x[,c(1,4)],
     type='n',
     panel.first=grid(col='black'),
     main=paste('PCA gene expression ',
                nrow(Xdata), 'samples *', ncol(Xdata), 'genes', sep=' '),
     xlab='PC1', ylab='PC4')
text(pc$x[,c(1,4)],labels=dimnames(pc$x)[[1]],col=as.numeric(data$batch),pch=0.5)
legend('bottomleft',col=as.numeric(unique(data$batch)),
       legend=paste('Batch', as.numeric(unique(data$batch))),pch=1,cex=0.7,bg='white',bty='o')

## Plot components PC2 and PC3
plot(pc$x[,c(2,3)],
     type='n',
     panel.first=grid(col='black'),
     main=paste('PCA gene expression ',
                nrow(Xdata), 'samples *', ncol(Xdata), 'genes', sep=' '),
     xlab='PC2', ylab='PC3')
text(pc$x[,c(2,3)],labels=dimnames(pc$x)[[1]],col=as.numeric(data$batch),pch=0.5)
legend('bottomleft',col=as.numeric(unique(data$batch)),
       legend=paste('Batch', as.numeric(unique(data$batch))),pch=1,cex=0.7,bg='white',bty='o')


## Plot components PC2 and PC4
plot(pc$x[,c(2,4)],
     type='n',
     panel.first=grid(col='black'),
     main=paste('PCA gene expression ',
                nrow(Xdata), 'samples *', ncol(Xdata), 'genes', sep=' '),
     xlab='PC2', ylab='PC4')
text(pc$x[,c(2,4)],labels=dimnames(pc$x)[[1]],col=as.numeric(data$batch),pch=0.5)
legend('bottomleft',col=as.numeric(unique(data$batch)),
       legend=paste('Batch', as.numeric(unique(data$batch))),pch=1,cex=0.7,bg='white',bty='o')

## Plot components PC3 and PC4
plot(pc$x[,c(3,4)],
     type='n',
     panel.first=grid(col='black'),
     main=paste('PCA gene expression ',
                nrow(Xdata), 'samples *', ncol(Xdata), 'genes', sep=' '),
     xlab='PC3', ylab='PC4')
text(pc$x[,c(3,4)],labels=dimnames(pc$x)[[1]],col=as.numeric(data$batch),pch=0.5)
legend('bottomleft',col=as.numeric(unique(data$batch)),
       legend=paste('Batch', as.numeric(unique(data$batch))),pch=1,cex=0.7,bg='white',bty='o')



## Plot components PC1 and PC2 coloured by ER status
plot(pc$x[,1:2],
     type='n',
     panel.first=grid(col='black'),
     main=paste('PCA gene expression ',
                nrow(Xdata), 'samples *', ncol(Xdata), 'genes', sep=' '),
     xlab='PC1', ylab='PC2')
text(pc$x[,1:2],labels=dimnames(pc$x)[[1]],col=as.numeric(data$ER_status)+1,pch=0.5)
legend('bottomleft',col=as.numeric(unique(data$ER_status)+1),
       legend=paste('ER', as.numeric(unique(data$ER_status))),pch=1,cex=0.7,bg='white',bty='o')



## Plot components PC1 and PC2 coloured by Nodal status
plot(pc$x[,1:2],
     type='n',
     panel.first=grid(col='black'),
     main=paste('PCA gene expression ',
                nrow(Xdata), 'samples *', ncol(Xdata), 'genes', sep=' '),
     xlab='PC1', ylab='PC2')
text(pc$x[,1:2],labels=dimnames(pc$x)[[1]],col=as.numeric(data$Nodal_status)+1,pch=0.5)
legend('bottomleft',col=as.numeric(unique(data$Nodal_status)+1),
       legend=paste('Nodal status', as.numeric(unique(data$Nodal_status))),pch=1,cex=0.7,bg='white',bty='o')

## Plot components PC1 and PC2 coloured by year of diagnosis
plot(pc$x[,1:2],
     type='n',
     panel.first=grid(col='black'),
     main=paste('PCA gene expression ',
                nrow(Xdata), 'samples *', ncol(Xdata), 'genes', sep=' '),
     xlab='PC1', ylab='PC2')
text(pc$x[,1:2],labels=dimnames(pc$x)[[1]],col=data$year_diag,pch=0.5)
legend('bottomleft',col=as.numeric(unique(data$year_diag)),
       legend=paste('Year of diagnosis', as.numeric(unique(data$year_diag))),pch=1,cex=0.7,bg='white',bty='o')

## Plot components PC1 and PC2 coloured by age_diag
colvec = 1
colvec[data$age_diag>55] = 2
plot(pc$x[,1:2],
     type='n',
     panel.first=grid(col='black'),
     main=paste('PCA gene expression ',
                nrow(Xdata), 'samples *', ncol(Xdata), 'genes', sep=' '),
     xlab='PC1', ylab='PC2')
text(pc$x[,1:2],labels=dimnames(pc$x)[[1]],col=unique(colvec),pch=0.5)
legend('bottomleft',col=unique(colvec),
       legend=paste('Age at diag', unique(colvec)),pch=1,cex=0.7,bg='white',bty='o')

## Plot components PC1 and PC2 coloured by tumor_size
colvec = 1
colvec[data$tumor_size>40] = 2
plot(pc$x[,1:2],
     type='n',
     panel.first=grid(col='black'),
     main=paste('PCA gene expression ',
                nrow(Xdata), 'samples *', ncol(Xdata), 'genes', sep=' '),
     xlab='PC1', ylab='PC2')
text(pc$x[,1:2],labels=dimnames(pc$x)[[1]],col=unique(colvec),pch=0.5)
legend('bottomleft',col=unique(colvec),
       legend=paste('Tumour size', unique(colvec)),pch=1,cex=0.7,bg='white',bty='o')


data$pc1 = pc$x[,1]
plot(data[,5:ncol(data)])

tmp = subset(data, select = c(ID, pc1))
dat = merge(dat, tmp, all.x = T)
par(ask = T, mfrow = c(1,1))
for (i in 70:(ncol(dat)-1)){ plot(dat$pc1,dat[,i], ylab = names(dat)[i]) }

#We have dependendy between PC1 and some biology variables, but not with array dates or similar.




################################################################
################################################################
#
# Genes expressions assoxiated with CN segments: models
#
################################################################
################################################################
require(parallel)


#############################################
# RUV model

library(ruv)

rinv<-RUVinv(Y=as.matrix((Xdata)), X=as.matrix(data$til, ncol = 1), ctl = rep(T, ncol(Xdata)), Z = NULL)
#hist(rinv$p,breaks=c(0,0.001,0.01,0.05,0.1,0.2,0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9, 1))
radj<-variance_adjust(rinv)
#hist(radj$p.rsvar,breaks=c(0,0.001,0.01,0.05,0.1,0.2,0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9, 1))
#hist(radj$p.ebayes,breaks=c(0,0.001,0.01,0.05,0.1,0.2,0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9, 1))
hst = hist(radj$p.evar,breaks=c(0,0.001,0.01,0.05,0.1,0.2,0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9, 1),
           main = 'Empirical variances expression p-values', xlab = 'p-values')
dev.print(png, file=paste(prefix.out,'results/Pvalue_hist_CN.png', sep = ''), 
          width=400, height=400)
qqnorm(rinv$betahat)
mtext('Expression ~ TIL coefficients')
dev.print(png, file=paste(prefix.out,'results/TILcoefficients_CN.png', sep = ''), 
          width=400, height=400)


#Index of the probe sets with p-value<0.01:
indexRUV = which(radj$p.evar<0.01)   
length(indexRUV) #[1] 454 segments



#radj has p-values (p.evar)
#indexRUV has probe sets with low p-values
#hst has p-value histogram

#############################################
# Multiple linear models

out = matrix(NA, ncol = 2, nrow = ncol(Xdata))
colnames(out) = c('Estimate', 'Pr(>|t|)')
for (i in 1:ncol(Xdata)){
  data$x = Xdata[,i]
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
indexLM = which(out[, 2]<0.01)  
length(indexLM) #[1] 231





################################################################
# Naive RUV adjustment, 4 different levels of adjustment


rnaive10 = naiveRandRuv(Y=as.matrix(Xdata), cIdx = rep(T, ncol(Xdata)), nuCoeff=1e-3, k=10)
res10 = matrix(NA, ncol = 2, nrow = ncol(Xdata))
colnames(res10) = c('Estimate', 'Pr(>|t|)')
for (i in 1:ncol(Xdata)){
  data$x = rnaive10[,i]
  mod = lm(til ~ x + age_diag + tumor_size + Nodal_status + ER_status, data = data)
  res10[i,] = summary(mod)$coef['x',colnames(res10)]
}
hist(res10[,2],breaks=c(0,0.001,0.01,0.05,0.1,0.2,0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9, 1),
     main = 'Naive RUV P-values', xlab = 'p-values')
mtext('k = 10')

rnaive = naiveRandRuv(Y=as.matrix((Xdata)), cIdx = rep(T, ncol(Xdata)), nuCoeff=1e-3, k=20)
res = matrix(NA, ncol = 2, nrow = ncol(Xdata))
colnames(res) = c('Estimate', 'Pr(>|t|)')
for (i in 1:ncol(Xdata)){
  data$x = rnaive[,i]
  mod = lm(til ~ x + age_diag + tumor_size + Nodal_status + ER_status, data = data)
  res[i,] = summary(mod)$coef['x',colnames(res)]
}
hist(res[,2],breaks=c(0,0.001,0.01,0.05,0.1,0.2,0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9, 1),
     main = 'Naive RUV P-values', xlab = 'p-values')
mtext('k = 20')


rnaive100 = naiveRandRuv(Y=as.matrix((Xdata)), cIdx = rep(T, ncol(Xdata)), nuCoeff=1e-1, k=20)
res100 = matrix(NA, ncol = 2, nrow = ncol(Xdata))
colnames(res100) = c('Estimate', 'Pr(>|t|)')
for (i in 1:ncol(Xdata)){
  data$x = rnaive100[,i]
  mod = lm(til ~ x + age_diag + tumor_size + Nodal_status + ER_status, data = data)
  res100[i,] = summary(mod)$coef['x',colnames(res100)]
}
hist(res100[,2],breaks=c(0,0.001,0.01,0.05,0.1,0.2,0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9, 1),
     main = 'Naive RUV P-values', xlab = 'p-values')
mtext('Times 100')

rnaive.01 = naiveRandRuv(Y=as.matrix((Xdata)), cIdx = rep(T, ncol(Xdata)), nuCoeff=1e-5, k=20)
res.01 = matrix(NA, ncol = 2, nrow = ncol(Xdata))
colnames(res.01) = c('Estimate', 'Pr(>|t|)')
for (i in 1:ncol(Xdata)){
  data$x = rnaive.01[,i]
  mod = lm(til ~ x + age_diag + tumor_size + Nodal_status + ER_status, data = data)
  res.01[i,] = summary(mod)$coef['x',colnames(res.01)]
}
hist(res.01[,2],breaks=c(0,0.001,0.01,0.05,0.1,0.2,0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9, 1),
     main = 'Naive RUV P-values', xlab = 'p-values')
mtext('Div 100')

#Index of the probe sets with p-value<0.01:
indexNR10 = which(res10[, 2]<0.01)  #BAD HISTOGRAM
length(indexNR10) #245 segments
length(intersect(indexNR10, indexRUV)) #112 shared with RUV

indexNR = which(res[, 2]<0.01)  #BAD HISTOGRAM
length(indexNR) #217
length(intersect(indexNR, indexRUV)) #130 shared with RUV

indexNR100 = which(res100[, 2]<0.01) 
length(indexNR100) #275
length(intersect(indexNR100, indexRUV)) #85 shared with RUV

indexNR.01 = which(res.01[, 2]<0.01)   #BAD HISTOGRAM
length(indexNR.01) #245
length(intersect(indexNR.01, indexRUV)) #73 shared with RUV

length(intersect(indexLM, indexNR100)) #155 shared with LM and 73 shared with RUV: Use the 155 as primary
length(intersect(indexLM, indexRUV)) #48, LM and RUV are very different


#Decided to use the intersect(indexNR100, indexLM) as toptable
#Print p-value histogram of the naive RUV with k = 20 and nu = 1e-1 (indexNR100 model)
hstNR = hist(res100[,2],breaks=c(0,0.001,0.01,0.05,0.1,0.2,0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9, 1),
             main = 'Segment CN naive RUV p-values', xlab = 'p-values')
mtext('Multiple models')
dev.print(png, file=paste(prefix.out,'results/Pvalue3_hist_CN.png', sep = ''), 
          width=400, height=400)


#############################################
# Rearrangements of Naive RUV 100 results, which will be considered primary
#(the 155 segments that intersect with the linear models will be preliminary primary)

topindex = intersect(indexLM, indexNR100)
topsegs = seginfo[topindex,]
topsegs$space = as.character(topsegs$space)
table(topsegs$space)
1  2  3  4  5  6  7  8  9 10 11 12 13 14 15 16 17 18 19 20 21 22 23 
23 43  0  2  0  6  2  3  2  0  3  5  0 30  0 13 12  0  0  5  0  0  6 

##### Read Refseq gene info
#Like to tie these segments to gene symbol, so we can find corresponding expression probes

plac = read.delim(paste(prefix.ann, "../refGene_hg19.txt",sep=''), header=F)
plac = subset(plac,  select=c('V3', 'V5','V6','V13'))
names(plac) = c('chrom','start','end','Gene.Symbol')
plac$chrom = as.character(plac$chrom)
plac$chrom[plac$chrom == 'chrX'] = 'chr23'
plac$chrom[plac$chrom == 'chrY'] = 'chr24'
plac$chrom[plac$chrom == 'chrM'] = 'chr25'
plac = subset(plac, !is.na(start) & !is.na(end) & !is.na(chrom))
plac$chrom = substr(plac$chrom, 4, nchar(plac$chrom, type = 'c'))
plac$width = plac$end - plac$start +1
oldind = mclapply(1:length(unique(plac$Gene.Symbol)), function(i){
  index = which(plac$Gene.Symbol == unique(plac$Gene.Symbol)[i])
  newindex = (index[order(-plac$width[index])])[1]
  setdiff(index, newindex)
})
oldind = unlist(oldind)
plac0 = plac[setdiff(which(!is.na(plac$chrom)), oldind),]   #Null set for GO Fisher       
ann = read.csv(paste(prefix.ann, 'HG-U133_Plus_2/HG-U133_Plus_2.na33.annot.csv',
                     sep = ''), skip = 25)
tmp = subset(ann, Gene.Symbol != '---', select = c('Gene.Symbol','Probe.Set.ID'))
plac = merge(plac0, tmp, all.x = T)


plac = RangedData(IRanges(start = plac$start, end = plac$end),  
                  space = plac$chrom, universe = 'hg19', 
                  Gene.Symbol = plac$Gene.Symbol, Probe.Set.ID= plac$Probe.Set.ID) 
topsegs = RangedData(IRanges(start = topsegs$start, end = topsegs$end),  
                     space = topsegs$space, universe = 'hg19')
hits = findOverlaps(query = plac, subject = topsegs)

mymax = 0
for (i in 1:length(names(hits))){
  nhits = countSubjectHits(hits[[names(hits)[i]]])
  m = max(nhits)
  if (m>mymax) mymax = m
}  #Warnings are OK
mymax  
[1] 31

load(paste(prefix.norm, 'Expression/Expression_estimates.RData', sep = ''))
GEdata = as.data.frame(t(as.matrix(expr[,6:ncol(expr)])))
index = which(!is.na(data$til))
GEdata = GEdata[index,] #The 51 CN samples' expression data only used
names(GEdata) = expr$unitName
rm(expr)
gc()

#############################################
#Toptable
toptable = data.frame('Segment.Number' = (1:nrow(res100))[topindex], 
                      'Chromosome' = seginfo$space[topindex],
                      'Start' = seginfo$start[topindex],
                      'End' = seginfo$end[topindex],
                      'Width' = seginfo$width[topindex],
                      'CN.Coef' = res100[topindex, 'Estimate'], 
                      p.value = res100[topindex, 'Pr(>|t|)'],
                      FDR = p.adjust(res100[,2], method = 'fdr')[topindex])

ge = matrix(NA, ncol = 4*mymax, nrow = nrow(toptable))
colnames(ge) = rep(c('Probe.Set.ID','Gene.Symbol','Mean.log2.expr','SD.log2.expr'), mymax)
ge = as.data.frame(ge)
for (i in 1:nrow(toptable)){
  chrom = toptable$Chromosome[i]
  x = hits[[as.character(chrom)]]
  myhits = queryHits(x)[subjectHits(x) == i]
  if (length(myhits) > 0){
    tmp = as.data.frame(plac[space(plac) == as.character(chrom),])[myhits,]
    tmp$Gene.Symbol = as.character(tmp$Gene.Symbol)
    tmp$Probe.Set.ID = as.character(tmp$Probe.Set.ID)
    ge[i, seq(1, 4*nrow(tmp), 4)] = tmp$Probe.Set.ID
    ge[i, seq(2, 4*nrow(tmp), 4)] = tmp$Gene.Symbol
    notna.ind  = which(!is.na(tmp$Probe.Set.ID))
    ex = log2(GEdata[,tmp$Probe.Set.ID[notna.ind]])
    if (is.null(dim(ex))) ex = as.matrix(ex, ncol = 1)
    ge[i, seq(3, 4*nrow(tmp), 4)[notna.ind]] = apply(ex, 2, mean)
    ge[i, seq(4, 4*nrow(tmp), 4)[notna.ind]] = apply(ex, 2, sd)
  }
}
plac1 = na.omit(unique(unlist(ge[,seq(2, ncol(ge), 4)]))) #Hit genes for GO analysis


toptable = cbind(toptable, ge)
toptable = toptable[order(toptable$FDR),]
toptable$rank = 1:nrow(toptable)
toptable = toptable[,c(ncol(toptable), 1:(ncol(toptable)-1))]
write.table(toptable, file = paste(prefix.out, 'TIL_by_CN_toptable.txt'), 
            row.names = F, quote = F, sep = '\t', na = '')

save.image(paste(prefix.out, 'TILs_on_CN.RData', sep = '')) 


#############################################
#Frequency plot
hist(Xdata, breaks = 1000)
quantile(Xdata, probs = c(.025, .05, .95, .975))
2.5%       5%      95%    97.5% 
  1.649629 1.734067 2.446999 2.787242 

cuts = 2 + c(-1, 1)*1.5


library(rtracklayer)
genomeinfo = SeqinfoForUCSCGenome(genome = 'hg19')
chrlens = genomeinfo@seqlengths[genomeinfo@seqnames %in% 
                                  paste('chr', c(as.character(1:22), 'X'), sep = '')]

###All 51 samples
nup = function(x){length(x[x>cuts[2]])}
ndown = function(x){length(x[x<cuts[1]])}
yup = apply(Xdata, 2, nup)
ydown = apply(Xdata, 2, ndown)
ymax = max(yup, ydown)


plot(1, 1, type = 'n', ylim = c(-ymax, ymax), xlim = c(0, sum(as.numeric(chrlens))),
     xaxt = 'n', ylab = 'Frequency', xlab = 'Genome position')
nsegs.chrom = table(seginfo$space)
xbase = 0
segbase = 0
abline(v = xbase)
for (i in 1:length(chrlens)){ 
  segind = (segbase + 1):(segbase + nsegs.chrom[i])
  lines(seginfo$mid[segind] + xbase, yup[segind], col = 'darkseagreen')
  lines(seginfo$mid[segind] + xbase, -ydown[segind], col = 'salmon')
  xbase = xbase + chrlens[i] 
  segbase = segbase + nsegs.chrom[i]
  abline(v = xbase)
}
abline(h = 0)
xmids = c(0,cumsum(as.numeric(chrlens)))
xmids = (xmids[-length(xmids)] + xmids[-1])/2
axis(1, at = xmids, labels = c(as.character(1:22), 'X'), tick = F, line = -1)

dev.print(png, file=paste(prefix.out,'results/CN_frequencies_all51samples', date,'.png', 
                          sep = ''), width=1200, height=500)


##########
#ER-

yup = apply(Xdata[data$ER_status == 0,], 2, nup)
ydown = apply(Xdata[data$ER_status == 0,], 2, ndown)
ymax = max(yup, ydown)

plot(1, 1, type = 'n', ylim = c(-11, 11), xlim = c(0, sum(as.numeric(chrlens))),
     xaxt = 'n', ylab = 'Frequency', xlab = 'Genome position')
nsegs.chrom = table(seginfo$space)
xbase = 0
segbase = 0
abline(v = xbase)
for (i in 1:length(chrlens)){ 
  segind = (segbase + 1):(segbase + nsegs.chrom[i])
  lines(seginfo$mid[segind] + xbase, yup[segind], col = 'darkseagreen')
  lines(seginfo$mid[segind] + xbase, -ydown[segind], col = 'salmon')
  xbase = xbase + chrlens[i] 
  segbase = segbase + nsegs.chrom[i]
  abline(v = xbase)
}
abline(h = 0)
xmids = c(0,cumsum(as.numeric(chrlens)))
xmids = (xmids[-length(xmids)] + xmids[-1])/2
axis(1, at = xmids, labels = c(as.character(1:22), 'X'), tick = F, line = -1)

dev.print(png, file=paste(prefix.out,'results/CN_frequencies_21ERnegsamples', date,'.png', 
                          sep = ''), width=1200, height=500)


##########
#ER+

yup = apply(Xdata[data$ER_status == 1,], 2, nup)
ydown = apply(Xdata[data$ER_status == 1,], 2, ndown)
ymax = max(yup, ydown)

plot(1, 1, type = 'n', ylim = c(-ymax, ymax), xlim = c(0, sum(as.numeric(chrlens))),
     xaxt = 'n', ylab = 'Frequency', xlab = 'Genome position')
nsegs.chrom = table(seginfo$space)
xbase = 0
segbase = 0
abline(v = xbase)
for (i in 1:length(chrlens)){ 
  segind = (segbase + 1):(segbase + nsegs.chrom[i])
  lines(seginfo$mid[segind] + xbase, yup[segind], col = 'darkseagreen')
  lines(seginfo$mid[segind] + xbase, -ydown[segind], col = 'salmon')
  xbase = xbase + chrlens[i] 
  segbase = segbase + nsegs.chrom[i]
  abline(v = xbase)
}
abline(h = 0)
xmids = c(0,cumsum(as.numeric(chrlens)))
xmids = (xmids[-length(xmids)] + xmids[-1])/2
axis(1, at = xmids, labels = c(as.character(1:22), 'X'), tick = F, line = -1)

dev.print(png, file=paste(prefix.out,'results/CN_frequencies_30ERpossamples', date,'.png', 
                          sep = ''), width=1200, height=500)


#Genomic positions of probe sets in expression arrays
plac = ann$Alignments
plac = strsplit(as.character(plac), split = '//', fixed = T)
first = function(x) x[1]
plac = mapply(first, plac)
plac = strsplit(plac, split = ' (', fixed = T)
plac = mapply(first, plac)
plac = strsplit(plac, split = ':', fixed = T)
second = function(x) x[2]
tmp = data.frame(chrom = mapply(first, plac), pos = mapply(second, plac),
                 stringsAsFactors = F)
tmp$chrom[tmp$chrom == 'chrX'] = 'chr23'
tmp$chrom[tmp$chrom == 'chrY'] = 'chr24'
tmp$chrom[tmp$chrom == 'chrM'] = 'chr25'
tmp$clen = nchar(tmp$chrom, type = 'c')
plac = data.frame(chrom = substr(tmp$chrom, 4, tmp$clen))
tmp = strsplit(tmp$pos, split = '-', fixed = T)
plac$start = as.numeric(mapply(first, tmp))
plac$end = as.numeric(mapply(second, tmp))
plac$chrom = as.character(plac$chrom)
plac$Gene.Symbol = ann$Gene.Symbol
plac$Probe.Set.ID = ann$Probe.Set.ID
plac = subset(plac, !is.na(start) & !is.na(end))

plac = RangedData(IRanges(start = plac$start, end = plac$end),  
                  space = plac$chrom, universe = 'hg19', Gene.Symbol = plac$Gene.Symbol,
                  Probe.Set.ID = plac$Probe.Set.ID) 










################################################################
################################################################
#
# CIN associated with TIL? Not done yet (requires integer CNs)
#
################################################################
################################################################

load("seg.mat.copy.RData")
source("Sherene.wGII.score.R")

# determine ploidy for segmented matrix
ploidy.samples <- sapply(unique(seg.mat.copy[,1]), fun.ploidy, seg.mat.copy)

# determine weighted GII for samples
wGII.samples  <- sapply(unique(seg.mat.copy[,1]), gii.weight.chrom, seg.mat.copy, ploidy.samples)


for (arr in c(4, 10, 14, 22)){
  samp = as.character((info$arrname[info$pass])[arr])
  fitname = paste('fit', arr, sep = '')
  filename <- paste(prefix.out, 'CBSfits/', fitname, '.RData', sep = '');
  fit = loadObject(filename)
  debug(estimateKappaByC1Density)
  kappa <- estimateKappa(fit, verbose=-10)
  
}
