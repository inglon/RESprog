##############################################################################
gii.weight.chrom <- function(x, seg, ploidy){
  # x: sample Name
  # seg: segmentation matrix, cols: Sample Name, Chromosome, Start Region, End Region, Number Probes, Copy Number Value.
  # ploidy: Vector containing ploidy of all samples (calculated using fun.ploidy)
  # The names of the ploidy vector must be the sample names
  # Example for function call: sapply(seg.mat[, 1], gii.weight.chrom, seg = seg.mat, ploidy = ploidy)
  gii.weight <- 0
  ploidy.x <- ploidy[x]
  sub <- subset(seg, seg[, 1] == x)
  
  for(cc in unique(seg[, 2])){
    sub.chrom <- subset(sub, sub[, 2] == cc)
    
    a <- sub.chrom[as.numeric(sub.chrom[, 6]) > ploidy.x | as.numeric(sub.chrom[, 6]) < ploidy.x, ]
    b <- sum(as.numeric(sub.chrom[, 4]) - as.numeric(sub.chrom[, 3])) 
    if(!is.null(dim(a)))
      gii.chrom <- sum(as.numeric(a[, 4]) - as.numeric(a[, 3]))/b
    else
      gii.chrom <-  sum(as.numeric(a[4]) - as.numeric(a[3]))/b 
    
    gii.weight <- gii.weight + gii.chrom	
  }
  return(gii.weight/length(unique(seg[, 2])))
}	


###########################################################################
fun.mode <- function(x, show_all_modes = FALSE)
{
  x_freq <- table(x)
  mode_location <- which.max(x_freq)
  The_mode <- names(x_freq)[mode_location]
  Number_of_modes <- length(mode_location)
  #
  if(show_all_modes) {
    if(Number_of_modes >1) {
      warning(paste("Multiple modes exist - returning all",Number_of_modes,"of them"))}
    return(The_mode)
  } 
  
  else {
    if(Number_of_modes >1) {
      warning(paste("Multiple modes exist - returning only the first one out of", Number_of_modes))}
    return(The_mode[1])
  }
}

###################################################################


fun.ploidy <- function(x, seg){
  # This function estimates the median ploidy from a segmentation matrix of PICNIC results (integer copy number)
  # by calculating the weighted median ploidy. The weights are equal to the segment length.
  # The variabel seg needs to me a segmentation matrix with columns in the follwing order:
  # Sample Name, Chromosome, Start Region, End Region, Number Probes, Copy Number Value.
  # The matrix for seg must be a character matrix and mustn't be a dataframe.
  # The variable x is should in include the unique sample names in the first column of seg .
  # The function returns the median ploidy for each sample.
  sub = subset(seg, seg[, 1] == x)
  require(limma)
  
  w.meds <- c()
  for (chr in unique(seg[,2]))
  {
    sub.chr <- subset(sub, sub[,2]==chr)
    w.med   <- weighted.median(as.numeric(sub.chr[, 6]), w = as.numeric(sub.chr[, 4]) - as.numeric(sub.chr[, 3])    )
    w.meds  <- c(w.meds,w.med)  
  }
  
  return(as.numeric(fun.mode(w.meds)))
  }