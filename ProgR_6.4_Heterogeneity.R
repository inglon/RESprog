################################################################
################################################################
# Script: ProgR_6.4_Heterogeneity.R
# Author: Ingrid Lonnstedt
# Date:  30/08/2013
# R version: R 3.0.0
# Details: Our own absolute CNs and heterogeneity
# Simplified code to run more samples
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

#Samples todo

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
info = read.delim(file.path(getwd(), 'AROMA','annotationData','Infokey.txt'), stringsAsFactors = F)
info$sample.name = format(info$sample.name, nsmall = 4, trim = T)

#Samples todo
#inds = c(1,2,3,4,5,6,8,9,10,11,13,14,18,19,20,21,24,26,30,32,33,34,35,36,38,39,40,41,42,43,45,46,48,49,50,52,53,55,56)

###Load data
#obj = loadSample(i, info = info)

###Read/create data (specify i)
i = 2
obj = list(i = i, d = file.path(getwd(), 'HET', info$sample.name[i]),
              info = info)
obj$h = ifelse(obj$i %in% c(1:3, 6, 10, 22, 28, 35, 42, 46, 48), 1, 2)
obj = readHAPSEG(obj, printfolder = obj$d)

obj = readCBS(obj, printfolder = obj$d)

####################
###Rotation
####################
obj$mystart = startRotation(obj)

obj = setWeights(obj)
#w = 1-exp(-mat$length/5000000) #Sample 56, 3

#You may have to adjust weights and rerun this:
obj = optimizeRotation(obj)

####################
### Choose purity with help from mutations
####################

obj = readMutations(obj)
plotAfByCN(obj, ns = 1:4) 
plotAfByPosition(obj, ns = 2:3, restrict = FALSE)
obj = setRotation(obj, n = 3)

#Choose n
n = 3 #Sample 24
n = 5 #Sample 30
n = 9 #Sample 3
n = 7 #Sample 56

####################
### Nonlinearity fix
### and final CNs, multiplicity
####################

linearize(obj, 'a1r', output = FALSE, logis = TRUE, asymporig = TRUE, plot = T)
linearize(obj, 'a2r', output = FALSE, logis = TRUE, asymporig = TRUE, plot = T)

obj$mat$a1rl = linearize(obj, 'a1r', output = TRUE, none = TRUE)
obj$linfunc1 = 'none'
obj$mat$a1rl = linearize(obj, 'a1r', output = TRUE, logis = TRUE, none = FALSE)
obj$linfunc1 = 'logis'

obj$mat$a2rl = linearize(obj, 'a2r', output = TRUE, logis = TRUE)
obj$linfunc2 = 'logis'
obj$mat$a2rl = linearize(obj, 'a2r', output = TRUE, asymporig = TRUE)
obj$linfunc2 = 'asymporig'

plotLin(obj)
obj = setCNs(obj, linearize = FALSE)
obj = setMultiplicity(obj)


####################
# Subclonality
####################

obj = subclonalMutations(obj)
obj = subclonalCNdist(obj, cut = 3) #Choose cut
obj = ITH(obj, cn.cut = 3.5, cn.cut0 = -.5)
plotCNmutsGrid(obj)
plotCNmutsGenome(obj)
plotCNmutsChrom(obj)




####################
### Output dataset details
####################

out = data.frame(ind = obj$i, sample = obj$info$sample[obj$i],
                 purity = obj$alpha, ploidy = obj$ploidy, ith = obj$ith,
                 n.mutations = nrow(obj$muts), n.nonsilent.mutations = 
                   nrow(obj$muts[obj$muts$effect == 'nonsilent',]))
write.table(out, file = paste(obj$d, '/',obj$info$sample.name[obj$i], 
                              '_details.tsv', sep = ''), sep = '\t', col.names = T,
            row.names = FALSE)
tmp = obj$muts[obj$muts$effect == 'nonsilent',]
mutout = data.frame(Hugo_Symbol = tmp$Hugo_Symbol, subclonal = tmp$subclonal)
write.table(mutout, file = paste(obj$d, '/',obj$info$sample.name[obj$i], 
                              '_nonsilent_mutations.tsv', sep = ''), sep = '\t', col.names = T,
            row.names = FALSE)

####################
### Output workspace
####################

#Save workspace
save(obj = obj,file = file.path(obj$d, paste(obj$info$sample.name[obj$i], '_obj.RData', sep = '')))

#Load
#d = file.path(getwd(), 'HET', info$sample.name[i])
#load(file = file.path(d, paste(info$sample.name[i], '.RData', sep = '')))

####################
### Output HAPSEG transformed version dataset
####################

seg.dat = hapseg.format(obj$seg.dat, obj$mat)
#Print the new object
d = file.path(getwd(), 'ABSOLUTE', 'WithNormals', 'HAPSEG2')
d = paste(d,'/', obj$info$sample.name[obj$i], sep = '')
save(seg.dat, file = file.path(d, paste(obj$info$sample.name[obj$i], '_', 
                                        obj$info$arrname[obj$i], 
                                        '_scaled_segdat.RData', sep = '')))



################################################################
################################################################
#End of file
