## The Broad Institute
## SOFTWARE COPYRIGHT NOTICE AGREEMENT
## This software and its documentation are copyright (2012) by the
## Broad Institute/Massachusetts Institute of Technology. All rights are
## reserved.
##
## This software is supplied without any warranty or guaranteed support
## whatsoever. Neither the Broad Institute nor MIT can be responsible for its
## use, misuse, or functionality.


AllelicMakeSegObj <- function(seg.dat.fn, gender, filter_segs=FALSE, min_probes=NA, max_sd=NA, verbose=FALSE) 
{   
  ## provides seg.dat
  if (!file.exists(seg.dat.fn)) {
    stop("seg.dat.fn does not exist")
  }

  seg.dat = read.delim( seg.dat.fn, row.names=NULL, stringsAsFactors=FALSE, check.names=FALSE)
  
 
  print("Gender not supported yet!!")
  gender = NA

  X.ix = seg.dat[,"Chromosome"]==23
  seg.dat[X.ix, "Chromosome"] = "X"

  nix= is.na(seg.dat[,"mu.minor"]) | is.na(seg.dat[,"mu.major"])
  seg.dat = seg.dat[!nix,]

  if( !is.na(gender) && gender %in% c( "Male", "Female") )
  {
     nix = seg.dat[,"Chromosome"] %in% c("Y", "M", "chrY", "chrM")
  }
  else
  {
     if( verbose ) {
       print("No or invalid gender specified - dropping X chromosome segs")
     }
     nix = seg.dat[,"Chromosome"] %in% c("X", "Y", "M", "chrX", "chrY", "chrM")
  }
  seg.dat = seg.dat[!nix, ]




 ## Construct allele.segs object for genome HSCR plotting
  allele.segs = seg.dat[, c("Chromosome", "Start.bp", "End.bp", "n_probes", "length", "mu.minor", "mu.major")]
  colnames(allele.segs) = c("Chromosome", "Start.bp", "End.bp", "n_probes", "length", "A1.Seg.CN", "A2.Seg.CN")
 

 ## reformat into concatenated seg-tab for each allele...

  as.1 <- seg.dat[, c("Chromosome", "Start.bp", "End.bp",
                                       "n_probes", "length", "mu.minor", "sigma.minor")]

  as.2 <- seg.dat[, c("Chromosome", "Start.bp", "End.bp",
                                       "n_probes", "length", "mu.major", "sigma.major")]
  
  colnames(as.1)[6] <- "copy_num"
  colnames(as.2)[6] <- "copy_num"
  
  colnames(as.1)[7] <- "seg_sigma"
  colnames(as.2)[7] <- "seg_sigma"

  ## remember each segment's mate
  seg.ix <- c(1:nrow(as.1))
  as.1 <- cbind(as.1, seg.ix)
  as.2 <- cbind(as.2, seg.ix)
  
  as.seg.dat <- rbind(as.1, as.2)
#  colnames(as.seg.dat)[7] <- "seg_sigma"
  colnames(as.seg.dat)[8] <- "seg.ix"
 

  if (filter_segs) {
    as.seg.dat = FilterSegs(as.seg.dat, min_probes=min_probes, max_sd=max_sd,
                                verbose=verbose)$seg.info
  }
  
  W <- as.numeric(as.seg.dat[, "length"])
  W <- W / sum(W)
  as.seg.dat <- cbind(as.seg.dat, W)
  colnames(as.seg.dat)[ncol(as.seg.dat)] <- "W"


  seg.obj = list()
  seg.obj$gender = gender
  seg.obj$as.seg.dat = as.seg.dat
  seg.obj$segtab = as.seg.dat
  seg.obj$allele.segs = allele.segs

  return(seg.obj)
}


AllelicExtractSampleObs <- function(seg.obj) 
{
  seg.dat <- seg.obj[["as.seg.dat"]]
  d <- seg.dat[, "copy_num"]
  stderr <- seg.dat[, "seg_sigma"]
  W <- seg.dat[, "W"]

  if("seg.ix" %in% colnames(seg.dat)) {
    seg.ix <- seg.dat[,"seg.ix"]
  } else {
    seg.ix <- NA
  }
  
  e.cr = sum(W*d )
  gender=seg.obj$gender
  
  if( !is.na(gender) && gender == "Male" ) 
  {
     male_X = (seg.dat[,"Chromosome"]=="X")
  }
  else
  {
     male_X = rep(FALSE, length(d) )
  }
#   print( paste(sum(male_X), " male_X segs, gender = ", gender, sep=""))

## TODO: purge error.model from all obs.   Get it from the parent seg.dat instead
  obs <- list(d=d, d.tx=d, d.stderr=stderr, W=W, seg.ix=seg.ix,
              error.model=list(), n_probes=seg.dat[,"n_probes"], e.cr=e.cr,
	      platform=seg.obj[["platform"]], data.type="ALLELIC", male_X=male_X )
  
  return(obs)
}


## In fact - this is the right function whenever segs have input std errors
Allelic_get_seg_sigma = function(SCNA_model, obs)
{
  seg_sigma =  exp(SCNA_model[["Theta"]]["sigma.A"]) * obs[["d.stderr"]] 
  sigma.h = SCNA_model[["sigma.h"]]

  seg_sigma = sqrt(seg_sigma^2 + sigma.h^2)
}





