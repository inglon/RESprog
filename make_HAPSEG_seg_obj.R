## The Broad Institute
## SOFTWARE COPYRIGHT NOTICE AGREEMENT
## This software and its documentation are copyright (2012) by the
## Broad Institute/Massachusetts Institute of Technology. All rights are
## reserved.
##
## This software is supplied without any warranty or guaranteed support
## whatsoever. Neither the Broad Institute nor MIT can be responsible for its
## use, misuse, or functionality.


## Note: this code should be deprecated by standardizing the output of HAPSEG (SNP6) to an
## allelic segtab TSV format.   Once this is done, make_allelic_seg_obj.R should replace this## file.


HAPSEGMakeSegObj <- function(seg.dat.fn, gender, filter_segs=FALSE, min_probes=NA, max_sd=NA, verbose=FALSE) 
{   
 AllelicInvAtten <- function(r, at) {
    return(r / (1 + at - (at * r)))
  }

  ## provides seg.dat
  if (!file.exists(seg.dat.fn)) {
    stop("seg.dat.fn does not exist")
  }
  load(seg.dat.fn)
  if (!exists("seg.dat")) {
    stop("Invalid object contained in seg.dat.fn")
  }

  ## Detect pre-1.0 HAPSEG inputs and convert
  if ("EM_res" %in% names(seg.dat)) {
    if (verbose) {
      print("Detected a pre-1.0 HAPSEG input, converting to 1.0 format ....")
    }
    seg.dat = RenameSegDatFields(seg.dat)
  }
  
  new.seg.info <- seg.dat[["allele.segs"]][, c("Chromosome",
                                               "Start.bp",
                                               "End.bp", "n_probes",
                                               "length", "A1.Seg.CN",
                                               "tCN.Seg.sd")]
  new.seg.info[, "A1.Seg.CN"] <- apply(seg.dat[["allele.segs"]][, c("A1.Seg.CN", "A2.Seg.CN")], 1, sum) / 2
  colnames(new.seg.info)[c(6, 7)] <- c("copy_num", "seg_sigma")
  seg.info <- new.seg.info

## tCR based filter
  seg.dat[["seg.info"]] <- seg.info
  seg.dat[["seg.info"]] <- ABSOLUTE:::FilterSegs(seg.dat[["seg.info"]],
                                      verbose=verbose)[["seg.info"]]
  
  W <- as.numeric(seg.dat[["seg.info"]][, "length"])
  W <- W / sum(W)
  seg.dat[["seg.info"]] <- cbind(seg.dat[["seg.info"]], W)
  colnames(seg.dat[["seg.info"]])[ncol(seg.dat[["seg.info"]])] <- "W"
  
  ## also create AS seg.info
  AS.Seg.sd = sqrt( seg.dat[["allele.segs"]][,"AS.Seg.sd"]^2 + seg.dat[["allele.segs"]][,"tCN.Seg.sd"]^2 ) / 2
#  AS.Seg.sd = seg.dat[["allele.segs"]][,"AS.Seg.sd"]




  as.1 <- seg.dat[["allele.segs"]][, c("Chromosome", "Start.bp", "End.bp",
                                       "n_probes", "length", "A1.Seg.CN")]

  as.2 <- seg.dat[["allele.segs"]][, c("Chromosome", "Start.bp", "End.bp",
                                       "n_probes", "length", "A2.Seg.CN")]
  
  as.1 = data.frame( as.1, "AS.Seg.sd"=AS.Seg.sd, check.names=FALSE, stringsAsFactors=FALSE)
  as.2 = data.frame( as.2, "AS.Seg.sd"=AS.Seg.sd, check.names=FALSE, stringsAsFactors=FALSE)

  colnames(as.1)[6] <- "copy_num"
  colnames(as.2)[6] <- "copy_num"
  
  ## remember each segment's mate
  seg.ix <- c(1:nrow(as.1))
  as.1 <- cbind(as.1, seg.ix)
  as.2 <- cbind(as.2, seg.ix)
  
  as.seg.dat <- rbind(as.1, as.2)
  colnames(as.seg.dat)[7] <- "seg_sigma"
  colnames(as.seg.dat)[8] <- "seg.ix"
  

  if( any(is.na(as.seg.dat[,"seg_sigma"])) )
  {
     print( paste("Capping ", sum(is.na(as.seg.dat[,"seg_sigma"])), " seg_sigma values at 10", sep=""))

     as.seg.dat[   is.na(as.seg.dat[,"seg_sigma"]), "seg_sigma" ] = 10
  }

  if (filter_segs) {
    as.seg.dat = ABSOLUTE:::FilterSegs(as.seg.dat, min_probes=min_probes, max_sd=max_sd,
                                verbose=verbose)$seg.info
  }
  
  W <- as.numeric(as.seg.dat[, "length"])
  W <- W / sum(W)
  as.seg.dat <- cbind(as.seg.dat, W)
  colnames(as.seg.dat)[ncol(as.seg.dat)] <- "W"

  if( !is.na(gender) && gender %in% c( "Male", "Female") )
  {
     nix = as.seg.dat[,"Chromosome"] %in% c("Y", "M", "chrY", "chrM")
  }
  else
  {
     if( verbose ) {
       print("No or invalid gender specified - dropping X chromosome segs")
     }
     nix = as.seg.dat[,"Chromosome"] %in% c("X", "Y", "M", "chrX", "chrY", "chrM")
  }


  seg.dat$gender = gender
  as.seg.dat = as.seg.dat[!nix, ]

  
  ## add error-model parameters from EM-fit
  em <- seg.dat[["em.res"]]
  
  if (is.null(seg.dat[["em.res"]][["theta"]])) {
    ## This indicates an older hapseg output, these values aren't
    ## contained by theta
    seg.dat[["error.model"]] <- list()
    seg.dat[["error.model"]][["sigma.eta"]] <- em[["sigma.eta"]]
    seg.dat[["error.model"]][["sigma.nu"]] <- em[["sigma.nu"]]
    seg.dat[["error.model"]][["het.cov"]] <- em[["het.cov"]]
    seg.dat[["error.model"]][["nu"]] <- em[["nu"]]
    seg.dat[["error.model"]][["AT"]] <- em[["at"]]
    seg.dat[["error.model"]][["BG"]]<- em[["bg"]]
  } else {
    seg.dat[["error.model"]] <- seg.dat[["em.res"]][["theta"]]
    seg.dat[["error.model"]][["AT"]] <- seg.dat[["em.res"]][["theta"]][["at"]]
    seg.dat[["error.model"]][["at"]] <- NULL
  }
  
  seg.dat[["error.model"]][["loglik"]] <- em[["theta"]][["loglik"]]
#  seg.dat[["error.model"]][["fit.at"]] <- em[["at"]]


## un-attenuate data
  AT <- seg.dat[["error.model"]][["AT"]]
  d = as.seg.dat[,"copy_num"]
  d.tx = AllelicInvAtten( d, AT )
  as.seg.dat[,"copy_num"] = d.tx

  seg.dat[["as.seg.dat"]] <- as.seg.dat
  seg.dat[["segtab"]] <- as.seg.dat


  snp.clust.p = seg.dat$em.res$snp.clust.p
  seg.dat[["BEAGLE.seg.expected.phase"]] = BuildSegExpectedPhase(length(snp.clust.p), snp.clust.p) 
  
  switches = unlist(lapply(seg.dat$h.switch.ix, length))
  seg.dat[["BEAGLE.seg.expected.phase"]][ switches == 0, ] =
  seg.dat[["seg.expected.phase"]][ switches == 0, ]


  ## Remove multiple fields from the resulting list as we no longer
  ## need them
  remove.fields <- c("em.res", "c.hat.1", "g.0",
                     "c.hat.phased.unfixed", "h.switch.ix",
                     "merged.loci", "final.merge.prob")
  seg.dat <- seg.dat[setdiff(names(seg.dat), remove.fields)]
  
  return(seg.dat)
}



BuildSegExpectedPhase <- function(n.segs, h.snp.clust.p) 
{
  seg.expected.phase <- matrix(NA, nrow=n.segs, ncol=2)

  for (i in seq(n.segs)) 
  {
    snp.clust.p <- h.snp.clust.p[[i]]
    state.mat <- matrix(c(0, 1, 2, 3), ncol=4, nrow=nrow(snp.clust.p),
                        byrow=TRUE)
    snp.e.state <- rowSums(state.mat * snp.clust.p[, c(3, 2, 4, 1)])
    ## now take expectation of phase over segment, weight by P(het)
    seg.expected.phase[i, 1] <- (sum(snp.clust.p[, 4] * snp.e.state) /
      sum(as.vector(snp.clust.p[, 4])))
    seg.expected.phase[i, 2] <- (sum(snp.clust.p[, 2] * snp.e.state) /
      sum(as.vector(snp.clust.p[, 2])))
  }

  return( seg.expected.phase )
}


