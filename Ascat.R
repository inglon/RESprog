ascat.loadData
function(Tumor_LogR_file, Tumor_BAF_file, Germline_LogR_file = NULL, 
         Germline_BAF_file = NULL, chrs = c(1:22,"X")) {
  
  # read in SNP array data files
  print.noquote("Reading Tumor LogR data...")
  Tumor_LogR <- read.table(Tumor_LogR_file, header=T, row.names=1, comment.char="", sep = "\t")
  print.noquote("Reading Tumor BAF data...")
  Tumor_BAF <- read.table(Tumor_BAF_file, header=T, row.names=1, comment.char="", sep = "\t")
  
  #infinite values are a problem - change those
  Tumor_LogR[Tumor_LogR==-Inf]=NA
  Tumor_LogR[Tumor_LogR==Inf]=NA
  
  Germline_LogR = NULL
  Germline_BAF = NULL
  if(!is.null(Germline_LogR_file)) {
    print.noquote("Reading Germline LogR data...")
    Germline_LogR <- read.table(Germline_LogR_file, header=T, row.names=1, comment.char="", sep = "\t")
    print.noquote("Reading Germline BAF data...")
    Germline_BAF <- read.table(Germline_BAF_file, header=T, row.names=1, comment.char="", sep = "\t")
    
    #infinite values are a problem - change those
    Germline_LogR[Germline_LogR==-Inf]=NA
    Germline_LogR[Germline_LogR==Inf]=NA
  }
  
  # make SNPpos vector that contains genomic position for all SNPs and remove all data not on chromosome 1-22,X
  print.noquote("Registering SNP locations...")
  SNPpos <- Tumor_LogR[,1:2]
  SNPpos = SNPpos[SNPpos[,1]%in%chrs,]
  
  Tumor_LogR = Tumor_LogR[rownames(SNPpos),c(-1,-2),drop=F]
  Tumor_BAF = Tumor_BAF[rownames(SNPpos),c(-1,-2),drop=F]
  
  if(!is.null(Germline_LogR_file)) {
    Germline_LogR = Germline_LogR[rownames(SNPpos),c(-1,-2),drop=F]
    Germline_BAF = Germline_BAF[rownames(SNPpos),c(-1,-2),drop=F]
  }
  
  # sort all data by genomic position
  last = 0;
  ch = list();
  SNPorder = vector(length=dim(SNPpos)[1])
  for (i in 1:length(chrs)) {
    chrke = SNPpos[SNPpos[,1]==chrs[i],]
    chrpos = chrke[,2]
    names(chrpos) = rownames(chrke)
    chrpos = sort(chrpos)
    ch[[i]] = (last+1):(last+length(chrpos))  
    SNPorder[ch[[i]]] = names(chrpos)
    last = last+length(chrpos)
  }
  SNPpos = SNPpos[SNPorder,]
  Tumor_LogR=Tumor_LogR[SNPorder,,drop=F]
  Tumor_BAF=Tumor_BAF[SNPorder,,drop=F]
  
  if(!is.null(Germline_LogR_file)) {
    Germline_LogR = Germline_LogR[SNPorder,,drop=F]
    Germline_BAF = Germline_BAF[SNPorder,,drop=F]
  }
  
  # write final SNP positions to file
  #write.table(SNPpos,"SNPpos.txt", sep="\t")
  
  # split the genome into distinct parts to be used for segmentation (e.g. chromosome arms, parts of genome between gaps in array design)
  print.noquote("Splitting genome in distinct chuncks...")
  chr = split_genome(SNPpos)
  
  return(list(Tumor_LogR = Tumor_LogR, Tumor_BAF = Tumor_BAF, 
              Tumor_LogR_segmented = NULL, Tumor_BAF_segmented = NULL, 
              Germline_LogR = Germline_LogR, Germline_BAF = Germline_BAF, 
              SNPpos = SNPpos, ch = ch, chr = chr, chrs = chrs, 
              samples = colnames(Tumor_LogR)))
}






ascat.predictGermlineGenotypes
function(ASCATobj, platform = "AffySNP6") {
  attach(ASCATobj)
  
  Homozygous = matrix(nrow = dim(Tumor_LogR)[1], ncol = dim(Tumor_LogR)[2])
  colnames(Homozygous) = colnames(Tumor_LogR)
  rownames(Homozygous) = rownames(Tumor_LogR)
  
  if (platform=="Custom10k") {
    maxHomozygous = 0.05
    proportionHetero = 0.59
    proportionHomo = 0.38
    proportionOpen = 0.02
    segmentLength = 20
  }
  else if (platform=="Illumina109k") {
    maxHomozygous = 0.05
    proportionHetero = 0.35
    proportionHomo = 0.60
    proportionOpen = 0.02
    segmentLength = 20
  }
  else if (platform=="IlluminaCytoSNP") {
    maxHomozygous = 0.05
    proportionHetero = 0.28
    proportionHomo = 0.62
    proportionOpen = 0.03
    segmentLength = 100
  }
  else if (platform=="Illumina610k") {
    maxHomozygous = 0.05
    proportionHetero = 0.295
    proportionHomo = 0.67
    proportionOpen = 0.015
    segmentLength = 30
  }
  else if (platform=="Illumina660k") {
    maxHomozygous = 0.05
    proportionHetero = 0.295
    proportionHomo = 0.67
    proportionOpen = 0.015
    segmentLength = 30
  }
  else if (platform=="Illumina700k") {
    maxHomozygous = 0.05
    proportionHetero = 0.295
    proportionHomo = 0.67
    proportionOpen = 0.015
    segmentLength = 30
  }
  else if (platform=="Illumina1M") {
    maxHomozygous = 0.05
    proportionHetero = 0.22
    proportionHomo = 0.74
    proportionOpen = 0.02
    segmentLength = 100
    #previousvalues:
    #proportionHetero = 0.24
    #proportionOpen = 0.01
    #segmentLength = 60
  }
  else if (platform=="Illumina2.5M") {
    maxHomozygous = 0.05
    proportionHetero = 0.21
    proportionHomo = 0.745
    proportionOpen = 0.03
    segmentLength = 100
  }
  else if (platform=="Affy100k") {
    maxHomozygous = 0.05
    proportionHetero = 0.27
    proportionHomo = 0.62
    proportionOpen = 0.09
    segmentLength = 30
  }
  else if (platform=="Affy250k_sty") {
    maxHomozygous = 0.05
    proportionHetero = 0.26
    proportionHomo = 0.66
    proportionOpen = 0.05
    segmentLength = 50
  }
  else if (platform=="Affy250k_nsp") {
    maxHomozygous = 0.05
    proportionHetero = 0.26
    proportionHomo = 0.66
    proportionOpen = 0.05
    segmentLength = 50
  }
  else if (platform=="Affy500k") {
    maxHomozygous = 0.05
    proportionHetero = 0.26
    proportionHomo = 0.66
    proportionOpen = 0.05
    segmentLength = 50
  }
  else if (platform=="AffySNP6") {
    maxHomozygous = 0.05
    proportionHetero = 0.25
    proportionHomo = 0.67
    proportionOpen = 0.04
    segmentLength = 100
  }
  else if (platform=="AffyOncoScan") {
    maxHomozygous = 0.05
    proportionHetero = 0.24
    proportionHomo = 0.69
    proportionOpen = 0.04
    segmentLength = 30
  }
  else {
    print("Error: platform unknown")
  }
  
  for (i in 1:dim(Tumor_LogR)[2]) {
    png(filename = paste("tumorSep",colnames(Tumor_LogR)[i],".png",sep=""), width = 2000, height = 500, res = 200)
    par(mar = c(0.5,5,5,0.5), cex = 0.4, cex.main=3, cex.axis = 2, pch = ifelse(dim(Tumor_LogR)[1]>100000,".",20))
    
    Tumor_BAF_noNA = Tumor_BAF[!is.na(Tumor_BAF[,i]),i]
    names(Tumor_BAF_noNA) = rownames(Tumor_BAF)[!is.na(Tumor_BAF[,i])]
    Tumor_LogR_noNA = Tumor_LogR[names(Tumor_BAF_noNA),i]
    names(Tumor_LogR_noNA) = names(Tumor_BAF_noNA)
    
    chr_noNA = list()
    prev = 0
    for(j in 1:length(chr)) {
      chrke = chr[[j]]
      next2 = prev + sum(!is.na(Tumor_BAF[chrke,i]))
      chr_noNA[[j]] = (prev+1):next2
      prev = next2
    }
    
    ch_noNA = list()
    prev = 0
    for(j in 1:length(ch)) {
      chrke = ch[[j]]
      next2 = prev + sum(!is.na(Tumor_BAF[chrke,i]))
      ch_noNA[[j]] = (prev+1):next2
      prev = next2
    }
    
    tbsam = Tumor_BAF_noNA
    # sample, mirrored
    bsm = ifelse(tbsam<0.5, tbsam, 1-tbsam)
    
    homoLimit = max(sort(bsm)[round(length(bsm)*proportionHomo)],maxHomozygous)
    
    Hom = ifelse(bsm<homoLimit,T,NA)
    
    Hetero = sum(Hom==F, na.rm=T)
    Homo = sum(Hom==T, na.rm=T)
    Undecided = sum(is.na(Hom))
    
    extraHetero = round(min(proportionHetero * length(Tumor_BAF_noNA) - Hetero,Undecided-proportionOpen*length(Tumor_BAF_noNA)))
    
    if(extraHetero>0) {
      
      allProbes=1:length(Tumor_BAF_noNA)
      nonHomoProbes = allProbes[is.na(Hom)|Hom==F]
      
      lowestDist = NULL
      
      # bsm, with homozygous replaced by NA
      bsmHNA=bsm
      bsmHNA[!is.na(Hom)&Hom]=NA
      
      for (chrke in chr_noNA) {
        
        chrNonHomoProbes = intersect(nonHomoProbes,chrke)
        
        # there must be a minimum number of probes on the chromosome, otherwise these are called homozygous anyway
        if (length(chrNonHomoProbes)>5) {        
          
          #make sure we're not going over any borders..
          segmentLength2 = min(length(chrNonHomoProbes)-1,segmentLength)
          
          chrNonHomoProbesStartWindowLeft = c(rep(NA,segmentLength2),chrNonHomoProbes[1:(length(chrNonHomoProbes)-segmentLength2)])
          chrNonHomoProbesEndWindowLeft = c(NA,chrNonHomoProbes[1:(length(chrNonHomoProbes)-1)])
          chrNonHomoProbesStartWindowRight = c(chrNonHomoProbes[2:length(chrNonHomoProbes)],NA)
          chrNonHomoProbesEndWindowRight = c(chrNonHomoProbes[(segmentLength2+1):length(chrNonHomoProbes)],rep(NA,segmentLength2))
          chrNonHomoProbesStartWindowMiddle = c(rep(NA,segmentLength2/2),chrNonHomoProbes[1:(length(chrNonHomoProbes)-segmentLength2/2)])
          chrNonHomoProbesEndWindowMiddle = c(chrNonHomoProbes[(segmentLength2/2+1):length(chrNonHomoProbes)],rep(NA,segmentLength2/2))
          
          chrLowestDist = NULL
          
          for (probeNr in 1:length(chrNonHomoProbes)) {
            probe = chrNonHomoProbes[probeNr]
            if(!is.na(chrNonHomoProbesStartWindowLeft[probeNr])&!is.na(chrNonHomoProbesEndWindowLeft[probeNr])) {
              medianLeft = median(bsmHNA[chrNonHomoProbesStartWindowLeft[probeNr]:chrNonHomoProbesEndWindowLeft[probeNr]], na.rm=T)
            }
            else {
              medianLeft = NA
            }
            if(!is.na(chrNonHomoProbesStartWindowRight[probeNr])&!is.na(chrNonHomoProbesEndWindowRight[probeNr])) {
              medianRight = median(bsmHNA[chrNonHomoProbesStartWindowRight[probeNr]:chrNonHomoProbesEndWindowRight[probeNr]], na.rm=T)
            }
            else {
              medianRight = NA
            }
            
            if(!is.na(chrNonHomoProbesStartWindowMiddle[probeNr])&!is.na(chrNonHomoProbesEndWindowMiddle[probeNr])) {
              medianMiddle = median(c(bsmHNA[chrNonHomoProbesStartWindowMiddle[probeNr]:chrNonHomoProbesEndWindowLeft[probeNr]],
                                      bsmHNA[chrNonHomoProbesStartWindowRight[probeNr]:chrNonHomoProbesEndWindowMiddle[probeNr]]), na.rm=T)
            }
            else {
              medianMiddle = NA
            }
            
            chrLowestDist[probeNr] = min(abs(medianLeft-bsm[probe]),abs(medianRight-bsm[probe]),abs(medianMiddle-bsm[probe]),Inf,na.rm=T)
          }
        }
        
        # if too few probes on the chromosome
        else {
          chrLowestDist = NULL
          if (length(chrNonHomoProbes)>0) {
            # 1 is higher than any practical distance
            chrLowestDist[1:length(chrNonHomoProbes)] = 1
          }
        }
        
        lowestDist = c(lowestDist,chrLowestDist)
      }
      
      lowestDistUndecided = lowestDist[is.na(Hom[nonHomoProbes])]
      names(lowestDistUndecided)=names(Tumor_LogR_noNA)[nonHomoProbes[is.na(Hom[nonHomoProbes])]]
      
      sorted = sort(lowestDistUndecided)
      Hom[names(sorted[1:min(length(sorted),extraHetero)])] = F
      
      Hetero = sum(Hom==F, na.rm=T)
      Homo = sum(Hom==T, na.rm=T)
      Undecided = sum(is.na(Hom))
      
    }
    
    title = paste(paste(colnames(Tumor_BAF)[i], Hetero), Homo)
    plot(c(1,length(Tumor_BAF_noNA)), c(0,1), type = "n", xaxt = "n", main = title, xlab = "", ylab = "")
    points(Tumor_BAF_noNA,col=ifelse(is.na(Hom),"green",ifelse(Hom,"blue","red")))
    
    abline(v=0.5,lty=1,col="lightgrey")
    chrk_tot_len = 0
    for (j in 1:length(ch_noNA)) {
      chrk = ch_noNA[[j]];
      chrk_tot_len_prev = chrk_tot_len
      chrk_tot_len = chrk_tot_len + length(chrk)
      vpos = chrk_tot_len;
      tpos = (chrk_tot_len+chrk_tot_len_prev)/2;
      text(tpos,1,chrs[j], pos = 1, cex = 2)
      abline(v=vpos+0.5,lty=1,col="lightgrey")
    }
    
    # set all Undecided to homozygous
    Hom[is.na(Hom)] = T
    
    dev.off()
    
    Homozygous[names(Hom),i] = Hom
  }
  
  detach(ASCATobj)
  
  return(Homozygous)
}
> 
  
  
  
  
  ascat.plotSegmentedData
function(ASCATobj, filenames = "ASPCF") {
  attach(ASCATobj)
  for (arraynr in 1:dim(Tumor_LogR)[2]) {
    Select_nonNAs = rownames(Tumor_BAF_segmented[[arraynr]])
    AllIDs = 1:dim(Tumor_LogR)[1]
    names(AllIDs) = rownames(Tumor_LogR)
    HetIDs = AllIDs[Select_nonNAs]
    png(filename = paste(filenames,samples[arraynr],".png",sep=""), width = 2000, height = 1000, res = 200)
    par(mar = c(0.5,5,5,0.5), mfrow = c(2,1), cex = 0.4, cex.main=3, cex.axis = 2)
    r = Tumor_LogR_segmented[rownames(Tumor_BAF_segmented[[arraynr]]),arraynr]
    beta = Tumor_BAF_segmented[[arraynr]][,]
    plot(c(1,length(r)), c(-1,1), type = "n", xaxt = "n", main = paste(colnames(Tumor_BAF)[arraynr],", LogR",sep=""), xlab = "", ylab = "")
    points(Tumor_LogR[rownames(Tumor_BAF_segmented[[arraynr]]),arraynr], col = "red", pch=ifelse(dim(Tumor_LogR)[1]>20000,".",20))
    points(r,col="green")
    abline(v=0.5,lty=1,col="lightgrey")
    chrk_tot_len = 0
    for (j in 1:length(ch)) {
      chrk = intersect(ch[[j]],HetIDs);
      chrk_tot_len_prev = chrk_tot_len
      chrk_tot_len = chrk_tot_len + length(chrk)
      vpos = chrk_tot_len;
      tpos = (chrk_tot_len+chrk_tot_len_prev)/2;
      text(tpos,1,chrs[j], pos = 1, cex = 2)
      abline(v=vpos+0.5,lty=1,col="lightgrey")
    }
    plot(c(1,length(beta)), c(0,1), type = "n", xaxt = "n", main = paste(colnames(Tumor_BAF)[arraynr],", BAF",sep=""), xlab = "", ylab = "")
    points(Tumor_BAF[rownames(Tumor_BAF_segmented[[arraynr]]),arraynr], col = "red", pch=ifelse(dim(Tumor_LogR)[1]>20000,".",20))
    points(beta, col = "green")
    points(1-beta, col = "green")
    abline(v=0.5,lty=1,col="lightgrey")
    chrk_tot_len = 0
    for (j in 1:length(ch)) {
      chrk = intersect(ch[[j]],HetIDs);
      chrk_tot_len_prev = chrk_tot_len
      chrk_tot_len = chrk_tot_len + length(chrk)
      vpos = chrk_tot_len;
      tpos = (chrk_tot_len+chrk_tot_len_prev)/2;
      text(tpos,1,chrs[j], pos = 1, cex = 2)
      abline(v=vpos+0.5,lty=1,col="lightgrey")
    }
    dev.off()
  }
  detach(ASCATobj)
}
> 
  plot(c(1,length(r)), c(-1,1), type = "n", xaxt = "n", main = paste(colnames(Tumor_BAF)[arraynr],", LogR",sep=""), xlab = "", ylab = "")
points(Tumor_LogR[rownames(Tumor_BAF_segmented[[arraynr]]),arraynr], col = "red", pch=ifelse(dim(Tumor_LogR)[1]>20000,".",20))                               points(r,col="green")                                                           abline(v=0.5,lty=1,col="lightgrey")                                             chrk_tot_len = 0                                                                for (j in 1:length(ch)) {
  chrk = intersect(ch[[j]],HetIDs);
  chrk_tot_len_prev = chrk_tot_len
  chrk_tot_len = chrk_tot_len + length(chrk)
  vpos = chrk_tot_len;
  tpos = (chrk_tot_len+chrk_tot_len_prev)/2;
  text(tpos,1,chrs[j], pos = 1, cex = 2)
  abline(v=vpos+0.5,lty=1,col="lightgrey")
}
plot(c(1,length(beta)), c(0,1), type = "n", xaxt = "n", main = paste(colnames(Tumor_BAF)[arraynr],", BAF",sep=""), xlab = "", ylab = "")                     points(Tumor_BAF[rownames(Tumor_BAF_segmented[[arraynr]]),arraynr], col = "red", pch=ifelse(dim(Tumor_LogR)[1]>20000,".",20))                                points(beta, col = "green")                                                     points(1-beta, col = "green")                                                   abline(v=0.5,lty=1,col="lightgrey")                                             chrk_tot_len = 0                                                                for (j in 1:length(ch)) {                                                         chrk = intersect(ch[[j]],HetIDs);                                               chrk_tot_len_prev = chrk_tot_len                                                chrk_tot_len = chrk_tot_len + length(chrk)                                      vpos = chrk_tot_len;                                                            tpos = (chrk_tot_len+chrk_tot_len_prev)/2;                                      text(tpos,1,chrs[j], pos = 1, cex = 2)                                          abline(v=vpos+0.5,lty=1,col="lightgrey")                                      }
dev.of  )
}
detach(ASCATobj)
}
> 
  plot(c(1,length(r)), c(-1,1), type = "n", xaxt = "n", main = paste(colnames(Tumor_BAF)[arraynr],", LogR",sep=""), xlab = "", ylab = "")
points(Tumor_LogR[rownames(Tumor_BAF_segmented[[arraynr]]),arraynr], col = "red", pch=ifelse(dim(Tumor_LogR)[1]>20000,".",20))
points(r,col="green")
abline(v=0.5,lty=1,col="lightgrey")
chrk_tot_len = 0
for (j in 1:length(ch)) {
  chrk = intersect(ch[[j]],HetIDs);
  chrk_tot_len_prev = chrk_tot_len
  chrk_tot_len = chrk_tot_len + length(chrk)
  vpos = chrk_tot_len;
  tpos = (chrk_tot_len+chrk_tot_len_prev)/2;
  text(tpos,1,chrs[j], pos = 1, cex = 2)
  abline(v=vpos+0.5,lty=1,col="lightgrey")
}
plot(c(1,length(beta)), c(0,1), type = "n", xaxt = "n", main = paste(colnames(Tumor_BAF)[arraynr],", BAF",sep=""), xlab = "", ylab = "")
points(Tumor_BAF[rownames(Tumor_BAF_segmented[[arraynr]]),arraynr], col = "red", pch=ifelse(dim(Tumor_LogR)[1]>20000,".",20))
points(beta, col = "green")
points(1-beta, col = "green")
abline(v=0.5,lty=1,col="lightgrey")
chrk_tot_len = 0
for (j in 1:length(ch)) {
  chrk = intersect(ch[[j]],HetIDs);
  chrk_tot_len_prev = chrk_tot_len
  chrk_tot_len = chrk_tot_len + length(chrk)
  vpos = chrk_tot_len;
  tpos = (chrk_tot_len+chrk_tot_len_prev)/2;
  text(tpos,1,chrs[j], pos = 1, cex = 2)
  abline(v=vpos+0.5,lty=1,col="lightgrey")
}
dev.off()
}
detach(ASCATobj)
}
> 
  plot(c(1,length(r)), c(-1,1), type = "n", xaxt = "n", main = paste(colnames(Tumor_BAF)[arraynr],", LogR",sep=""), xlab = "", ylab = "")
points(Tumor_LogR[rownames(Tumor_BAF_segmented[[arraynr]]),arraynr], col = "red", pch=ifelse(dim(Tumor_LogR)[1]>20000,".",20))
points(r,col="green")
abline(v=0.5,lty=1,col="lightgrey")
chrk_tot_len = 0
for (j in 1:length(ch)) {
  chrk = intersect(ch[[j]],HetIDs);
  chrk_tot_len_prev = chrk_tot_len
  chrk_tot_len = chrk_tot_len + length(chrk)
  vpos = chrk_tot_len;
  tpos = (chrk_tot_len+chrk_tot_len_prev)/2;
  text(tpos,1,chrs[j], pos = 1, cex = 2)
  abline(v=vpos+0.5,lty=1,col="lightgrey")
}
plot(c(1,length(beta)), c(0,1), type = "n", xaxt = "n", main = paste(colnames(Tumor_BAF)[arraynr],", BAF",sep=""), xlab = "", ylab = "")
points(Tumor_BAF[rownames(Tumor_BAF_segmented[[arraynr]]),arraynr], col = "red", pch=ifelse(dim(Tumor_LogR)[1]>20000,".",20))
points(beta, col = "green")
points(1-beta, col = "green")
abline(v=0.5,lty=1,col="lightgrey")
chrk_tot_len = 0
for (j in 1:length(ch)) {
  chrk = intersect(ch[[j]],HetIDs);
  chrk_tot_len_prev = chrk_tot_len
  chrk_tot_len = chrk_tot_len + length(chrk)
  vpos = chrk_tot_len;
  tpos = (chrk_tot_len+chrk_tot_len_prev)/2;
  text(tpos,1,chrs[j], pos = 1, cex = 2)
  abline(v=vpos+0.5,lty=1,col="lightgrey")
}
dev.off()
}
detach(ASCATobj)
}
> 
  plot(c(1,length(r)), c(-1,1), type = "n", xaxt = "n", main = paste(colnames(Tumor_BAF)[arraynr],", LogR",sep=""), xlab = "", ylab = "")
points(Tumor_LogR[rownames(Tumor_BAF_segmented[[arraynr]]),arraynr], col = "red", pch=ifelse(dim(Tumor_LogR)[1]>20000,".",20))
points(r,col="green")
abline(v=0.5,lty=1,col="lightgrey")
chrk_tot_len = 0
for (j in 1:length(ch)) {
  chrk = intersect(ch[[j]],HetIDs);
  chrk_tot_len_prev = chrk_tot_len
  chrk_tot_len = chrk_tot_len + length(chrk)
  vpos = chrk_tot_len;
  tpos = (chrk_tot_len+chrk_tot_len_prev)/2;
  text(tpos,1,chrs[j], pos = 1, cex = 2)
  abline(v=vpos+0.5,lty=1,col="lightgrey")
}
plot(c(1,length(beta)), c(0,1), type = "n", xaxt = "n", main = paste(colnames(Tumor_BAF)[arraynr],", BAF",sep=""), xlab = "", ylab = "")
points(Tumor_BAF[rownames(Tumor_BAF_segmented[[arraynr]]),arraynr], col = "red", pch=ifelse(dim(Tumor_LogR)[1]>20000,".",20))
points(beta, col = "green")
points(1-beta, col = "green")
abline(v=0.5,lty=1,col="lightgrey")
chrk_tot_len = 0
for (j in 1:length(ch)) {
  chrk = intersect(ch[[j]],HetIDs);
  chrk_tot_len_prev = chrk_tot_len
  chrk_tot_len = chrk_tot_len + length(chrk)
  vpos = chrk_tot_len;
  tpos = (chrk_tot_len+chrk_tot_len_prev)/2;
  text(tpos,1,chrs[j], pos = 1, cex = 2)
  abline(v=vpos+0.5,lty=1,col="lightgrey")
}
dev.off()
}
detach(ASCATobj)
}
> 

  
  
  ascat.runAscat
function(ASCATobj, gamma = 0.55, sunrisefiles = "sunrise", profilefiles = "ASCATprofile", rawprofilefiles = "rawprofile", allow100percent = F) {
  attach(ASCATobj)
  goodarrays=NULL
  res = vector("list",dim(Tumor_LogR)[2])
  for (arraynr in 1:dim(Tumor_LogR)[2]) {
    print.noquote(paste("Sample ", samples[arraynr], " (",arraynr,"/",length(samples),")",sep=""))
    lrr=Tumor_LogR[,arraynr]
    names(lrr)=rownames(Tumor_LogR)
    baf=Tumor_BAF[,arraynr]
    names(baf)=rownames(Tumor_BAF)
    lrrsegm = Tumor_LogR_segmented[,arraynr]
    names(lrrsegm) = rownames(Tumor_LogR_segmented)
    bafsegm = Tumor_BAF_segmented[[arraynr]][,]
    names(bafsegm) = rownames(Tumor_BAF_segmented[[arraynr]])
    res[[arraynr]] = runASCAT(lrr,baf,lrrsegm,bafsegm,ch,paste(sunrisefiles,samples[arraynr],".png",sep=""),paste(profilefiles,samples[arraynr],".png",sep=""),paste(rawprofilefiles,samples[arraynr],".png",sep=""),gamma,allow100percent)
    if(!is.na(res[[arraynr]]$rho)) {
      goodarrays[length(goodarrays)+1] = arraynr
    }
  }
  
  if(length(goodarrays)>0) {
    n1 = matrix(nrow = dim(Tumor_LogR)[1], ncol = length(goodarrays))
    n2 = matrix(nrow = dim(Tumor_LogR)[1], ncol = length(goodarrays))
    rownames(n1) = rownames(Tumor_LogR)
    rownames(n2) = rownames(Tumor_LogR)
    colnames(n1) = colnames(Tumor_LogR)[goodarrays]
    colnames(n2) = colnames(Tumor_LogR)[goodarrays]
    for (i in 1:length(goodarrays)) {
      n1[,i] = res[[goodarrays[i]]]$nA
      n2[,i] = res[[goodarrays[i]]]$nB
    }
    
    tp = vector(length=length(goodarrays))
    psi = vector(length=length(goodarrays))
    ploidy = vector(length=length(goodarrays))
    for (i in 1:length(goodarrays)) {
      tp[i] = res[[goodarrays[i]]]$rho
      psi[i] = res[[goodarrays[i]]]$psi
      ploidy[i] = mean(res[[goodarrays[i]]]$nA+res[[goodarrays[i]]]$nB,na.rm=T)
    }
    fa = colnames(Tumor_LogR)[-goodarrays]
  }
  else {
    n1 = NULL
    n2 = NULL
    tp = NULL
    ploidy = NULL
    psi = NULL
    fa = colnames(Tumor_LogR)
  }
  
  seg = list()
  for (i in 1:length(samples)) {
    seg[[i]] = res[[i]]$seg
  }
  
  detach(ASCATobj)
  
  return(list(nA = n1, nB = n2, aberrantcellfraction = tp, ploidy = ploidy, psi = psi, failedarrays = fa, segments = seg))
}
