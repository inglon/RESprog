> RunAbsolute = function (seg.dat.fn, sigma.p, max.sigma.h, min.ploidy, max.ploidy, 
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
  seg.dat = MakeSegObj(seg.dat.fn, min_probes = min_probes, 
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
    mode.res = myapply_subclonal_scna_model(segobj = seg.dat, mode_res = mode.res, 
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
      mode.res = ABSOLUTE:::ApplySomaticMutsModel(mode.res, seg.dat$obs_SCNA, 
                                       mut.cn.dat, kPiSomThetaQ, mut_class_w, kQ, verbose = verbose)
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
<environment: namespace:ABSOLUTE>
  
  
  
  
  > ABSOLUTE:::CreateMutCnDat
function (maf, seg.dat, min.mut.af, verbose = FALSE) 
{
  mut.cn.dat <- maf
  if ("total_normals_called" %in% colnames(mut.cn.dat)) {
    ix <- mut.cn.dat[, "total_normals_called"] > 1
    if (verbose) {
      print(paste("Removing ", sum(ix), " of ", length(ix), 
                  " mutations due to seen in > 1 normals", sep = ""))
    }
    mut.cn.dat <- mut.cn.dat[!ix, ]
  }
  if ("dbSNP_Val_Status" %in% colnames(mut.cn.dat)) {
    mut.cn.dat[["dbSNP_Val_Status"]][is.na(mut.cn.dat[["dbSNP_Val_Status"]])] <- ""
  }
  cols <- colnames(mut.cn.dat)
  cix <- which(cols %in% c("i_t_ref_count", "t_ref_count"))
  if (length(cix) > 0) {
    colnames(mut.cn.dat)[cix] <- "ref"
  }
  else {
    stop("Malformed MAF file, no ref column supplied")
  }
  cix <- which(cols %in% c("i_t_alt_count", "t_alt_count"))
  if (length(cix) > 0) {
    colnames(mut.cn.dat)[cix] <- "alt"
  }
  else {
    stop("Malformed MAF file, no alt column supplied")
  }
  cix <- which(cols == "dbSNP_Val_Status")
  if (length(cix) > 0) {
    colnames(mut.cn.dat)[cix] <- "dbSNP"
  }
  else {
    stop("Malformed MAF file, no dbSNP_Val_Status column supplied")
  }
  cix <- which(cols == "Tumor_Sample_Barcode")
  if (length(cix) > 0) {
    colnames(mut.cn.dat)[cix] <- "sample"
  }
  else {
    stop("Malformed MAF file, no Tumor_Sample_Barcode column supplied")
  }
  na.ix <- apply(is.na(mut.cn.dat[, c("ref", "alt")]), 1, sum) > 
    0
  if (verbose) {
    print(paste("Removing ", sum(na.ix), " of ", length(na.ix), 
                " mutations with NA coverage", sep = ""))
  }
  mut.cn.dat <- mut.cn.dat[!na.ix, ]
  af <- mut.cn.dat[, "alt"]/(mut.cn.dat[, "alt"] + mut.cn.dat[, 
                                                              "ref"])
  ix <- af < min.mut.af
  if (verbose) {
    print(paste("Removing ", sum(ix), " of ", length(ix), 
                " mutations due to allelic fraction < ", min.mut.af, 
                sep = ""))
  }
  mut.cn.dat <- mut.cn.dat[!ix, ]
  ix = mut.cn.dat[, "Hugo_Symbol"] %in% c("ADAM6")
  if (verbose) {
    print(paste("Removing ", sum(ix), " mutations in IG regions", 
                sep = ""))
  }
  if (sum(!ix) == 0) {
    stop("no mutations left!")
  }
  mut.cn.dat = mut.cn.dat[!ix, , drop = FALSE]
  mut.seg.ix <- GetMutSegIx(mut.cn.dat, seg.dat[["obs.scna"]][["segtab"]])
  ix <- apply(is.na(mut.seg.ix), 1, sum) == 0
  if (verbose && (sum(!ix) > 0)) {
    print(paste("Removing ", sum(!ix), " unmapped mutations on Chrs: ", 
                sep = ""))
    print(mut.cn.dat[!ix, "Chromosome"])
  }
  if (sum(ix) == 0) {
    stop("No mutations left")
  }
  mut.cn.dat <- mut.cn.dat[ix, ]
  mut.seg.ix <- mut.seg.ix[ix, , drop = FALSE]
  mut.cn.dat <- cbind(mut.cn.dat, mut.seg.ix)
  return(mut.cn.dat)
}
<environment: namespace:ABSOLUTE>
  
  > ABSOLUTE:::ApplySomaticMutsModel
function (mode.res, obs_scna, mut.cn.dat, pi.som.theta.q, mut_class_w, 
          Q, verbose = FALSE) 
{
  q <- dim(mode.res[["theta.q.tab"]])[2]
  if (verbose) {
    print(paste("Evaluating ", nrow(mut.cn.dat), " mutations over ", 
                nrow(mode.res[["mode.tab"]]), " purity/ploidy modes: ", 
                sep = ""))
  }
  for (j in 1:nrow(mode.res[["mode.tab"]])) {
    alpha <- mode.res[["mode.tab"]][j, "alpha"]
    seg.qz.hat <- mode.res[["seg.qz.tab"]][j, , ]
    seg.q.hat <- mode.res[["seg.q.tab"]][j, , ]
    subclonal_scna_tab = mode.res$subclonal_SCNA_res$subclonal_SCNA_tab[j, 
                                                                        , ]
    ccf_dens = mode.res$subclonal_SCNA_res$CCF_dens[j, , 
                                                    ]
    res = FitSomaticMuts(mut.cn.dat, mode.res$mode.tab[j, 
                                                       ], subclonal_scna_tab, ccf_dens, obs_scna, seg.qz.hat, 
                         seg.q.hat, pi.som.theta.q, mut_class_w, Q, verbose = verbose)
    muts.post.prs <- res[["muts.post.prs"]]
    som.theta.q.map <- res[["som.theta.q.map"]]
    if (j == 1) {
      mode.res[["muts.post.prs"]] <- array(NA, dim = c(nrow(mut.cn.dat), 
                                                       ncol(muts.post.prs) + 1, nrow(mode.res[["mode.tab"]])))
      dimnames(mode.res[["muts.post.prs"]])[[2]] <- c(colnames(muts.post.prs), 
                                                      "purity")
      dimnames(mode.res[["muts.post.prs"]])[[1]] <- rownames(muts.post.prs)
    }
    mode.res[["muts.post.prs"]][, , j] <- cbind(muts.post.prs, 
                                                purity = alpha)
    mode.res[["mode.tab"]][j, "somatic.mut.ll"] <- sum(muts.post.prs[, 
                                                                     "LL"], na.rm = TRUE) + dDirichlet(som.theta.q.map, 
                                                                                                       pi.som.theta.q, log.p = TRUE)
    if (verbose) {
      cat(".")
    }
  }
  if (verbose) {
    cat("\n")
  }
  new.ll <- mode.res[["mode.tab"]][, "clust_LL"] + mode.res[["mode.tab"]][, 
                                                                          "log.partition"] + mode.res[["mode.tab"]][, "somatic.mut.ll"]
  mode.res[["mode.tab"]][, "post_LL"] <- new.ll
  return(mode.res)
}
<environment: namespace:ABSOLUTE>
  
  >
  
  > ABSOLUTE:::FitSomaticMuts
function (mut.cn.dat, mode_info, subclonal_scna_tab, scna_ccf_dens, 
          obs_scna, seg.qz.hat, seg.q.hat, pi.som.theta.q, mut_class_w, 
          Q, verbose = FALSE) 
{
  mut.cn.dat = get_muts_nearest_clonal_scna(mut.cn.dat, seg.qz.hat, 
                                            seg.q.hat, Q)
  mult_res = dir_post_fit_somatic_multiplicity(mut.cn.dat, 
                                               mode_info, mut_class_w, pi.som.theta.q, print_fit = verbose)
  som_theta_q_map = mult_res$som_theta_Q_MAP
  post_prs = mult_res$post_Prs
  clonal_scna_mut_ix = !get_subclonal_scna_mut_ix(mut.cn.dat, 
                                                  subclonal_scna_tab)
  if (FALSE) {
    clonal_scna_mult_res = dir_post_fit_somatic_multiplicity(mut.cn.dat[clonal_scna_mut_ix, 
                                                                        ], mode_info, mut_class_w, pi.som.theta.q, print_fit = verbose)
    som_theta_q_map = clonal_scna_mult_res$som_theta_Q_MAP
    post_prs = matrix(NA, nrow = nrow(mut.cn.dat), ncol = ncol(clonal_scna_mult_res$post_prs))
    colnames(post_prs) = colnames(clonal_scna_mult_res$post_Prs)
    post_prs[clonal_scna_mut_ix, ] = clonal_scna_mult_res$post_Prs
    subclonal_scna_mult_res = calc_sample_muts_on_subclonal_scna(mut.cn.dat[!clonal_scna_mut_ix, 
                                                                            ], subclonal_scna_tab, scna_ccf_dens, mode_info, 
                                                                 obs_scna, mut_class_w, som_theta_q_map)
    post_prs[!clonal_scna_mut_ix, ] = subclonal_scna_mult_res$post_Prs
  }
  var_classes = ClassifySomaticVariants(post_prs, 0.5)
  q_s = post_prs[, "modal_q_s"]
  mut_ccf_res = calc_snv_ccf(mut.cn.dat, mode_info, q_s, post_prs[, 
                                                                  "Pr_somatic_clonal"])
  muts.post.prs <- cbind(as.matrix(mut.cn.dat[, c("q_hat", 
                                                  "HS_q_hat_1", "HS_q_hat_2"), drop = FALSE]), post_prs, 
                         var_classes, mut_ccf_res, subclonal_SCNA = !clonal_scna_mut_ix)
  return(list(muts.post.prs = muts.post.prs, som.theta.q.map = som_theta_q_map))
}
<environment: namespace:ABSOLUTE>
  
  
>  ABSOLUTE:::calc_snv_ccf  
  function (mut_dat, mode_info, modal_q_s, pr_somatic_clonal) 
  {
    muts = mut_dat
    alpha = mode_info["alpha"]
    cov = rowSums(muts[, c("alt", "ref")])
    f = muts[, "alt"]/cov
    Q = muts[, "q_hat"]
    q_s = rep(NA, length(pr_somatic_clonal))
    c.ix = pr_somatic_clonal > 0.5
    q_s[c.ix] = modal_q_s[c.ix]
    q_s[!c.ix] = 1
    som_delta = alpha/(2 * (1 - alpha) + alpha * Q)
    f_s = q_s * som_delta
    mut_scale = 1/f_s
    cell_frac = f * mut_scale
    cell_mult = f/som_delta
    ccf_grid = seq(0, 1, by = 0.01)
    ccf_dens = matrix(NA, nrow = nrow(mut_dat), ncol = length(ccf_grid))
    ccf_ci95 = matrix(NA, nrow = nrow(ccf_dens), ncol = 2)
    ccf_hat = rep(NA, nrow(mut_dat))
    for (i in seq_len(nrow(mut_dat))) {
      if (mut_dat[i, "q_hat"] == 0) {
        next
      }
      ccf_dens[i, ] = calc_ccf_posterior_grid(mut_dat[i, "alt"], 
                                              mut_dat[i, "ref"], alpha, mut_dat[i, "q_hat"], ccf_grid)
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
    ix = mut_dat[i, "q_hat"] == 0
    ccf_ci95[ix, ] = NA
    res = cbind(cell_mult, cell_frac, ccf_hat, ccf_ci95)
    colnames(res) = c("cell_mult", "old_cancer_cell_frac", "cancer_cell_frac", 
                      "ccf_CI95_low", "ccf_CI95_high")
    return(res)
  }
<environment: namespace:ABSOLUTE>
  
  
  > ABSOLUTE:::dir_post_fit_somatic_multiplicity
function (mut_cn_dat, mode_info, mut_class_w, pi_som_theta_q, 
          print_fit = TRUE) 
{
  som_theta_Q_mode = (pi_som_theta_q - 1)/(sum(pi_som_theta_q) - 
                                             length(pi_som_theta_q))
  mcv = c(mut_class_w$Pi_SM, mut_class_w$Pi_SC)
  mcw = (mcv - 1)/(sum(mcv) - length(mcv))
  mut_class_w$SM = mcw[1]
  mut_class_w$SC = mcw[2]
  mult_res = CalcSampleMutsPostPr(mut_cn_dat, mode_info, mut_class_w, 
                                  som_theta_Q_mode)
  mut_mat = mult_res$som_mut_Q_tab * matrix(mult_res$post_Prs[, 
                                                              "Pr_somatic_clonal"], nrow = nrow(mult_res$som_mut_Q_tab), 
                                            ncol = ncol(mult_res$som_mut_Q_tab), byrow = TRUE)
  nix = mult_res$post_Prs[, "Pr_somatic_clonal"] == 0
  mut_mat = mut_mat[!nix, , drop = FALSE]
  if (nrow(mut_mat) > 0) {
    som_q = colSums(mut_mat, na.rm = TRUE)
  }
  else {
    som_q = rep(0, length(pi_som_theta_q))
  }
  som_theta_q_map = (pi_som_theta_q + som_q - 1)/(sum(pi_som_theta_q + 
                                                        som_q) - length(pi_som_theta_q))
  if (print_fit) {
    cat("som_theta_Q_MAP: ")
    print(round(som_theta_q_map, 5))
  }
  mut_w = mult_res$post_Prs[, c("Pr_somatic_clonal", "Pr_subclonal"), 
                            drop = FALSE]
  mut_w = mut_w/rowSums(mut_w)
  mut_pr = colSums(mut_w, na.rm = TRUE)
  mcw = (mcv + mut_pr - 1)/(sum(mcv + mut_pr) - length(mcv))
  mut_class_w$SM = mcw[1]
  mut_class_w$SC = mcw[2]
  if (print_fit) {
    print(mut_class_w)
  }
  mult_res = CalcSampleMutsPostPr(mut_cn_dat, mode_info, mut_class_w, 
                                  som_theta_q_map)
  mult_res$som_theta_Q_MAP = som_theta_q_map
  return(mult_res)
}
<environment: namespace:ABSOLUTE>
  
  
  > ABSOLUTE:::ApplySomaticMutsModel
function (mode.res, obs_scna, mut.cn.dat, pi.som.theta.q, mut_class_w, 
          Q, verbose = FALSE) 
{
  q <- dim(mode.res[["theta.q.tab"]])[2]
  if (verbose) {
    print(paste("Evaluating ", nrow(mut.cn.dat), " mutations over ", 
                nrow(mode.res[["mode.tab"]]), " purity/ploidy modes: ", 
                sep = ""))
  }
  for (j in 1:nrow(mode.res[["mode.tab"]])) {
    alpha <- mode.res[["mode.tab"]][j, "alpha"]
    seg.qz.hat <- mode.res[["seg.qz.tab"]][j, , ]
    seg.q.hat <- mode.res[["seg.q.tab"]][j, , ]
    subclonal_scna_tab = mode.res$subclonal_SCNA_res$subclonal_SCNA_tab[j, 
                                                                        , ]
    ccf_dens = mode.res$subclonal_SCNA_res$CCF_dens[j, , 
                                                    ]
    res = FitSomaticMuts(mut.cn.dat, mode.res$mode.tab[j, 
                                                       ], subclonal_scna_tab, ccf_dens, obs_scna, seg.qz.hat, 
                         seg.q.hat, pi.som.theta.q, mut_class_w, Q, verbose = verbose)
    muts.post.prs <- res[["muts.post.prs"]]
    som.theta.q.map <- res[["som.theta.q.map"]]
    if (j == 1) {
      mode.res[["muts.post.prs"]] <- array(NA, dim = c(nrow(mut.cn.dat), 
                                                       ncol(muts.post.prs) + 1, nrow(mode.res[["mode.tab"]])))
      dimnames(mode.res[["muts.post.prs"]])[[2]] <- c(colnames(muts.post.prs), 
                                                      "purity")
      dimnames(mode.res[["muts.post.prs"]])[[1]] <- rownames(muts.post.prs)
    }
    mode.res[["muts.post.prs"]][, , j] <- cbind(muts.post.prs, 
                                                purity = alpha)
    mode.res[["mode.tab"]][j, "somatic.mut.ll"] <- sum(muts.post.prs[, 
                                                                     "LL"], na.rm = TRUE) + dDirichlet(som.theta.q.map, 
                                                                                                       pi.som.theta.q, log.p = TRUE)
    if (verbose) {
      cat(".")
    }
  }
  if (verbose) {
    cat("\n")
  }
  new.ll <- mode.res[["mode.tab"]][, "clust_LL"] + mode.res[["mode.tab"]][, 
                                                                          "log.partition"] + mode.res[["mode.tab"]][, "somatic.mut.ll"]
  mode.res[["mode.tab"]][, "post_LL"] <- new.ll
  return(mode.res)
}
<environment: namespace:ABSOLUTE>
  
  >
  
  > ABSOLUTE:::GetCallStatus
function (mode.res, seg.w) 
{
  status = "called"
  q.tab = mode.res[["seg.qz.tab"]][1, , ]
  q.tab = q.tab[, c(1:(ncol(q.tab)))]
  Q = ncol(q.tab) - 1
  max.q = (apply(q.tab, 1, which.max))
  peak_masses = rep(0, Q)
  for (i in 1:Q) {
    ix = which(max.q == i)
    peak_masses[i] = sum(seg.w[ix])
  }
  b = mode.res[["mode.tab"]][1, "b"]
  if (peak_masses[1] < 0.01 & b < 0.15) {
    peak_masses[1] = 0
  }
  six = order(peak_masses, decreasing = TRUE)
  if (peak_masses[six[3]] < 1e-04) {
    if (peak_masses[six[2]] < 1e-04) {
      if (mode.res[["mode.tab"]][1, "sigma.h.hat"] < 0.02 & 
            !is.null(mode.res[["muts.post.prs"]])) {
        status = "non-aneuploid"
      }
      else {
        status = "low purity"
      }
    }
  }
  if (mode.res[["mode.tab"]][1, "entropy"] > 0.2) {
    status = "high entropy"
  }
  if (mode.res[["mode.tab"]][1, "frac.het"] > 0.2) {
    status = "high non-clonal"
  }
  return(status)
}
<environment: namespace:ABSOLUTE>
  
  
  
  
  
myCreateReviewObject = function (obj.name, absolute.files, indv.results.dir, copy_num_type, 
          plot.modes = TRUE, verbose = FALSE, myModes = NULL) 
{
  if (copy_num_type == "total") {
    ABSOLUTE:::set_total_funcs()
  } else if (copy_num_type == "allelic") {
    ABSOLUTE:::set_allelic_funcs()
  } else {
    stop("Unsupported copy number type: ", copy_num_type)
  }
  dir.create(indv.results.dir, recursive = TRUE)
  nm = file.path(indv.results.dir, paste(obj.name, ".PP-modes", 
                                         sep = ""))
  modesegs.fn = paste(nm, ".data.RData", sep = "")
  pdf.fn = paste(nm, ".plots.pdf", sep = "")
  failed.pdf.fn = paste(nm, "FAILED_plots.pdf", sep = ".")
  failed.tab.fn = paste(nm, "FAILED_tab.txt", sep = ".")
  call.tab.fn = file.path(indv.results.dir, paste(obj.name, 
                                                  "PP-calls_tab.txt", sep = "."))
  segobj.list = vector(mode = "list", length = length(absolute.files))
  failed.list = vector(mode = "list", length = length(absolute.files))
  so.ix = 1
  fa.ix = 1
  sids = character()
  for (i in seq_along(absolute.files)) {
    seg.out.fn = absolute.files[i]
    if (!file.exists(seg.out.fn)) {
      if (verbose) {
        cat("\n")
        print(paste("sample #", i, " result not found", 
                    sep = ""))
      }
      next
    }
    load(absolute.files[i])
    SID = seg.dat$sample.name
    if (SID %in% sids) {
      stop("Sample names must be unique")
    }
    sids = c(sids, SID)
    if (is.null(seg.dat$array.name)) {
      seg.dat$array.name = SID
    }
    if (is.na(seg.dat[["mode.res"]][["mode.flag"]])) {
      segobj.list[[so.ix]] = seg.dat
      names(segobj.list)[so.ix] = SID
      so.ix = so.ix + 1
      if (verbose) {
        cat(".")
      }
    }    else {
      if (verbose) {
        print(paste("Failed list got updated", i))
      }
      failed.list[[fa.ix]] = seg.dat
      names(failed.list)[fa.ix] = SID
      failed.list[[fa.ix]][["sample.name"]] = seg.dat[["sample.name"]]
      fa.ix = fa.ix + 1
      if (verbose) {
        cat("-")
      }
    }
  }
  segobj.list = segobj.list[c(1:(so.ix - 1))]
  failed.list = failed.list[c(1:(fa.ix - 1))]
  mode.ent = rep(NA, length(segobj.list))
  names(mode.ent) = names(segobj.list)
  for (i in seq_along(segobj.list)) {
    mtab = segobj.list[[i]][["mode.res"]][["mode.tab"]]
    ix = which.max(mtab[, "post_LL"])
    mode.ent[i] = mtab[ix, "entropy"]
  }
  samples = names(sort(mode.ent))
  segobj.list = segobj.list[samples]
  save(segobj.list, file = modesegs.fn)
  ABSOLUTE:::PrintPpCallTable(segobj.list, call.tab.fn)
  if (plot.modes) {
    n.print = length(myModes)
    myPlotModes(segobj.list, pdf.fn, n.print = n.print, myModes = myModes)
  }
  if (!is.null(failed.list[[1]])) {
    ABSOLUTE:::PlotFailedSamples(failed.list, failed.pdf.fn)
    ABSOLUTE:::PrintFailedTable(failed.list, failed.tab.fn)
  }
  return(TRUE)
}
  
  

myPlotModes = function (segobj.list, pdf.fn, n.print = NA, write.tab = TRUE, 
          debug.info = TRUE, add = FALSE, fig.mode = FALSE, sideways = FALSE, 
          low.w = FALSE, called.mode.ix = NA, verbose = FALSE, myModes = myModes) 
{
  pal <- colorRampPalette(c("chocolate4", "cyan4"))
  if (fig.mode) {
    add <- TRUE
    debug.info <- FALSE
  }
  Q <- 15
  alpha.dom <- c(0, 1)
  tau.dom <- c(1, 10)
  mode.colors <- ABSOLUTE:::GetModeColors()
  x.max <- 2.25
  binW <- 0.025
  if (!add) {
    pdf(pdf.fn, 3.5 * (1 + 2), 15)
    par(mfrow = c(5, 4))
  }
  for (s in seq_len(length(segobj.list))) {
    mode.tab <- segobj.list[[s]][["mode.res"]][["mode.tab"]]
    seg.z.tab <- segobj.list[[s]][["mode.res"]][["seg.z.tab"]]
    if (nrow(mode.tab) == 0) {
      for (i in seq_len(n.print + 2)) {
        frame()
      }
      next
    }
    obs <- ExtractSampleObs(segobj.list[[s]])
    SN <- segobj.list[[s]][["sample.name"]]
    if (is.na(n.print)) {
      s.n.print <- nrow(mode.tab)
    }
    else {
      s.n.print <- min(nrow(mode.tab), n.print)
    }
    s.mix <- s.n.print
    if (!fig.mode) {
      model.id <- segobj.list[[s]][["group"]]
      ABSOLUTE:::PostPlot(mode.tab[myModes,], mode.colors, alpha.dom, tau.dom, 
               sample.name = SN, obs, called.mode.ix, debug.info, call.status = 
                            segobj.list[[s]][["mode.res"]][["call.status"]], 
               model.id = model.id)
      ABSOLUTE:::PpModeScorerBarplot(mode.tab[myModes,], mode.colors, obs, n.print)
      if (length(segobj.list) == 1) {
        frame()
        frame()
      }
    }
    n.plot <- min(n.print, NROW(mode.tab))
     for (i in seq_len(n.plot)) {
      tau <- mode.tab[myModes[i], "tau"]
      alpha <- mode.tab[myModes[i], "alpha"]
      delta <- alpha/(2 * (1 - alpha) + alpha * tau)
      b <- 2 * (1 - alpha)/(2 * (1 - alpha) + alpha * tau)
      mode.info <- mode.tab[myModes[i], ]
      obs[["error.model"]][["fit.at"]] <- mode.tab[myModes[i], "AT"]
      comb <- InvTxData(obs[["error.model"]], GetCopyRatioComb(Q, 
                                            delta, b, obs[["error.model"]]))
      seg.dat <- InvTxData(obs[["error.model"]], obs[["d.tx"]])
      last.plot <- ifelse(i == n.plot, TRUE, FALSE)
      mid.plot <- ifelse(i == ceiling(n.plot/2), TRUE, 
                         FALSE)
      first.plot <- ifelse(i == 1, TRUE, FALSE)
      if ((length(segobj.list) > 1) && (!is.null(segobj.list[[s]][["mut.cn.dat"]])) && 
            (i > 1)) {
        frame()
        frame()
      }
      if (!is.null(segobj.list[[s]][["mode.res"]][["ab.tab"]])) {
        comb.ab <- segobj.list[[s]][["mode.res"]][["ab.tab"]][myModes[i], ]
      }
      else {
        comb.ab <- NA
      }
      ABSOLUTE:::ModeCombPlot(seg.dat, obs[["W"]], obs[["d.stderr"]], 
                   mode.info, seg.z.tab[myModes[i], ], comb.ab, mode.colors[i], 
                   comb, last.plot, mid.plot, first.plot, x.max, 
                   debug.info, fig.mode, sideways, pal)
      if (!is.na(called.mode.ix)) {
        if (i == called.mode.ix) {
          mtext(text = "*", side = 3, at = 0, col = "red", 
                line = -1, cex = 3 * par("cex") * par("cex.axis"))
        }
      }
      if (!is.null(segobj.list[[s]][["mut.cn.dat"]])) {
        mut.cn.dat <- segobj.list[[s]][["mut.cn.dat"]]
        modeled <- segobj.list[[s]][["mode.res"]][["muts.post.prs"]][, 
                                                     , myModes[i], drop = TRUE]
        if (is.null(dim(modeled))) {
          modeled <- matrix(modeled, nrow = 1)
          colnames(modeled) <- dimnames(segobj.list[[s]][["mode.res"]][["muts.post_prs"]])[[2]]
        }
        modeled.mut.dat <- cbind(mut.cn.dat, modeled)
        ABSOLUTE:::PlotSomaticMutDensities(modeled.mut.dat, segobj.list[[s]], 
                                1, mode.colors[i], min.cov = 3, verbose = verbose)
      }
    }
    if ((length(segobj.list) > 1) && (is.null(segobj.list[[s]][["mut.cn.dat"]]))) {
      if (i < n.print) {
        for (j in seq_len(n.print - i)) {
          frame()
        }
      }
      frame()
    }
  }
  if (!add) {
    dev.off()
  }
}


ABSOLUTE:::ModeCombPlot
function (d, W, seg.stderr, mode.info, seg.z, comb.ab, mode.color, 
          comb, last.plot, mid.plot, first.plot, x.max, debug.info, 
          fig.mode, sideways, pal) 
{
  ix <- (d >= 0) & (d < x.max)
  colpal <- pal(1000)
  col <- FALSE
  old.mar <- par("mar")
  xlab = "Copy ratio"
  if (debug.info) {
    par(mar = old.mar + c(0, 0, 1, 0))
  }
  if (fig.mode) {
    x.max <- 2
    ix <- (d >= 0) & (d < x.max)
    if (col) {
      x.ax.labs <- ifelse(last.plot, TRUE, FALSE)
    }
    else {
      x.ax.labs <- TRUE
    }
    if (col) {
      xlab <- ifelse(last.plot, "Allelic copy ratio", "")
    }
    else {
      xlab <- ifelse(first.plot, "Allelic copy ratio", 
                     "")
    }
    if (col) {
      ylab <- "Genomic fraction"
    }
    else {
      ylab <- ifelse(mid.plot, "Genomic fraction", "")
    }
    if (!col) {
      y.ax.labs <- first.plot
    }
    else {
      y.ax.labs <- last.plot
    }
    if (sideways) {
      tmp <- xlab
      xlab <- ylab
      ylab <- tmp
    }
    if (sideways) {
      tmp <- x.ax.labs
      x.ax.labs <- y.ax.labs
      y.ax.labs <- tmp
    }
    if (!fig.mode) {
      order.by <- W[ix]
    }
    else {
      order.by <- seg.z[ix]
    }
    PlotSeglenHist(d[ix], W[ix], color.by = seg.z[ix], order.by = order.by, 
                   color.range = c(0, 1), use.pal = colpal, bin.w = 0.04, 
                   x.max = x.max, data.whiskers = debug.info, ylab = ylab, 
                   xlab = xlab, x.ax.labs = x.ax.labs, y.ax.labs = y.ax.labs, 
                   sideways = sideways, border = !fig.mode)
    if (!col) {
      if (first.plot) {
        mtext("Candidate interpretations of copy profile", 
              side = 3, line = 0, adj = 0, cex = par("cex.axis"))
      }
    }
  }
  else {
    PlotSeglenHist(d[ix], W[ix], color.by = seg.z[ix], color.range = c(0, 
                                                                       1), use.pal = colpal, bin.w = 0.025, x.max = x.max, 
                   data.whiskers = FALSE, xlab = xlab, sideways = sideways)
  }
  msg1 <- substitute(paste(purity(hat(alpha)) == x1, ", ", 
                           "ploidy: ", hat(tau) == x2, ", ", hat(tau)[g] == x3, 
                           sep = ""), list(x1 = round(mode.info["alpha"], 2), x2 = round(mode.info["tau"], 
                                                                                         2), x3 = round(mode.info["genome mass"], 2)))
  msg2 <- substitute(paste(hat(sigma[H]) == x1, ", ", hat(theta[Z]) == 
                             x2, ", ", "%", Het == x3, sep = ""), list(x1 = round(mode.info["sigma.h.hat"], 
                                                                                  3), x2 = round(mode.info["theta.z.hat"], 2), x3 = round(100 * 
                                                                                                                                            mode.info["frac.het"], 2)))
  msg3 <- substitute(paste("SCNAs" == x6, ", ", "Kar" == x7, 
                           ", ", "SSNVs" == x8, sep = ""), list(x6 = round(mode.info["log.partition"], 
                                                                           2), x7 = round(mode.info["clust_LL"], 2), x8 = round(mode.info["somatic.mut.ll"], 
                                                                                                                                2)))
  if (debug.info) {
    mtext(msg1, line = 4, cex = 0.75, adj = 0)
    mtext(msg2, line = 3, cex = 0.75, adj = 0)
    mtext(msg3, line = 2, cex = 0.75, adj = 0)
  }
  extra <- (comb[2] - comb[1])/2
  for (q in seq_along(comb)) {
    if (comb[q] > x.max) {
      break
    }
    q_ix <- ifelse(q > 1, q - 1, q)
    seg.ix <- d < comb[q_ix] - extra
    cum_gen <- sum(W[seg.ix])
    if (cum_gen > 0.99) {
      break
    }
  }
  max.q <- q - 1
  top.axis <- FALSE
  if (top.axis) {
    axis(side = 3, at = comb[c(1:max.q)], labels = c(1:max.q) - 
           1, col = mode.color, col.ticks = mode.color)
  }
  for (q in seq_along(comb)) {
    if (q > max.q) {
      next
    }
    if (!sideways) {
      abline(v = comb[q], lwd = 1.5, lty = 3, col = mode.color)
    }
    else {
      abline(h = comb[q], lwd = 1.5, lty = 3, col = mode.color)
    }
    side <- ifelse(!sideways, 3, 4)
    if (!top.axis) {
      if (q - 1 < 10 | (q - 1)%%2 == 1) {
        mtext(text = (q - 1), side = side, at = comb[q], 
              col = mode.color, line = 0.2, cex = par("cex") * 
                par("cex.axis"))
      }
    }
    if (debug.info & !is.na(comb.ab)) {
      q.ab <- round(comb.ab[q], 2)
      if (is.finite(q.ab) & q.ab > 0) {
        mtext(text = paste("(", q.ab, ")", sep = ""), 
              at = comb[q], col = 1, line = 1, cex = 0.5)
      }
    }
  }
  par(mar = old.mar)
}


PlotSeglenHist = function (D, W, color.by, color.range = NA, x.max = NA, y.max = NA, 
          bin.w = 0.05, order.by = NA, use.pal = NA, data.whiskers = TRUE, 
          xlab = "Observed allelic copy-ratio", x.axis = TRUE, y.axis = TRUE, 
          x.ax.labs = TRUE, y.ax.labs = TRUE, xlim = NA, ylab = "Genomic fraction", 
          border = TRUE, sideways = FALSE, add = FALSE) 
{
  P <- length(D)
  if (is.na(x.max)) {
    x.max <- max(D)
  }
  if (length(bin.w) == 1) {
    res <- hist(D, breaks = seq(0, x.max + bin.w, bin.w), 
                plot = FALSE)
    breaks <- res[["breaks"]]
  }
  else {
    breaks <- bin.w
  }
  heights <- matrix(0, ncol = length(breaks), nrow = length(D))
  seg.colors <- matrix(0, ncol = length(breaks), nrow = length(D))
  if (is.na(use.pal)) {
    col.scale <- 1000
    use.pal <- heat.colors(col.scale)
  }
  else {
    col.scale <- length(use.pal)
  }
  if (is.na(color.range)) {
    color.range <- range(color.by, na.rm = TRUE)
  }
  pal.idx <- floor((color.by - color.range[1])/(color.range[2] - 
                                                  color.range[1]) * (col.scale - 1))
  for (p in seq_len(P)) {
    bin <- max(which(breaks <= D[p]))
    heights[p, bin] <- W[p]
    seg.colors[p, bin] <- use.pal[pal.idx[p] + 1]
  }
  seg.colors[is.na(seg.colors)] <- "grey"
  for (i in seq_len(ncol(heights))) {
    if (is.na(order.by)) {
      res <- sort(heights[, i], decreasing = TRUE, index.return = TRUE)
    }
    else {
      res <- sort(order.by, index.return = TRUE)
    }
    heights[, i] <- heights[res[["ix"]], i]
    seg.colors[, i] <- seg.colors[res[["ix"]], i]
  }
  seg.colors[heights == 0] <- NA
  heights[heights == 0] <- NA
  MyBarplot(heights, seg.colors, ylab, xlab, x.max, y.max, 
            add = add, breaks = breaks, xlim = xlim, border = border, 
            sideways = sideways)
  if (data.whiskers) {
    side <- ifelse(!sideways, 1, 2)
    axis(side = side, at = D, labels = FALSE, line = -0.5, 
         col = "grey")
    axis(side = side, at = D[!is.na(color.by)], labels = FALSE, 
         line = -0.5, col = "red")
    x.axis.line <- 1
  }
  else {
    x.axis.line <- NA
  }
  if (x.axis) {
    if (!sideways) {
      axis(side = 1, line = x.axis.line, labels = x.ax.labs)
    }
    else {
      axis(side = 2, line = x.axis.line, labels = x.ax.labs)
    }
  }
  if (y.axis) {
    if (!sideways) {
      par(ylog = FALSE)
      axis(side = 2, labels = y.ax.labs, las = 2)
      par(ylog = FALSE)
    }
    else {
      par(xlog = FALSE)
      axis(side = 1, labels = y.ax.labs, las = 1)
      par(xlog = FALSE)
    }
  }
}



ABSOLUTE:::FitSample
function (seg.obj, mut, Q, pi_theta_qz, sigma.h, tau.dom, sigma.h.dom, 
          chr.arms.dat, verbose = FALSE) 
{
  kThetaQz = rep(1/(Q + 1), Q + 1)
  obs = seg.obj[["obs.scna"]]
  res = FindLocationModes(obs, mut, Q, kThetaQz, sigma.h, tau.dom, 
                          verbose = verbose)
  if (!is.null(res[["mode.flag"]])) {
    return(list(mode.flag = res[["mode.flag"]]))
  }
  else {
    mode.tab = res[["mode.tab"]]
  }
  if (is.na(mode.tab)) {
    return(list(mode.flag = "ERROR"))
  }
  n.modes = nrow(mode.tab)
  theta.qz.hat = matrix(NA, nrow = n.modes, ncol = Q + 1)
  theta.q.tab = array(NA, dim = c(n.modes, Q))
  if (verbose) {
    print(paste("Optimizing LL(data, theta.qz, sigma.h | comb) for ", 
                n.modes, " modes: ", sep = ""))
  }
  seg.z.tab = array(NA, dim = c(n.modes, length(obs[["W"]])))
  seg.qz.tab = array(NA, dim = c(n.modes, length(obs[["W"]]), 
                                 Q + 1))
  seg.q.tab = array(NA, dim = c(n.modes, length(obs[["W"]]), 
                                Q))
  if (obs[["data.type"]] == "ALLELIC") {
    ab.tab = array(NA, dim = c(n.modes, Q + 1))
    chr.arm.tab = array(NA, dim = c(n.modes, 2, nrow(chr.arms.dat), 
                                    Q))
  }
  if (obs[["data.type"]] == "TOTAL") {
    ab.tab = NULL
    chr.arm.tab = array(NA, dim = c(n.modes, 1, nrow(chr.arms.dat), 
                                    Q))
  }
  dimnames(chr.arm.tab)[[3]] = rownames(chr.arms.dat)
  mode.tab = cbind(mode.tab, NA, NA, NA, NA, NA, NA, NA, NA)
  colnames(mode.tab)[c((ncol(mode.tab) - 7):ncol(mode.tab))] = c("genome mass", 
                                                                 "sigma.h.hat", "theta.z.hat", "frac.het", "log.partition", 
                                                                 "entropy", "clust_LL", "post_LL")
  mode.tab = cbind(mode.tab, somatic.mut.ll = NA)
  for (i in 1:n.modes) {
    delta = mode.tab[i, "delta"]
    b = mode.tab[i, "b"]
    obs[["error.model"]][["fit.at"]] = mode.tab[i, "AT"]
    comb = GetCopyRatioComb(Q, delta, b, obs[["error.model"]])
    res = OptThetaQzSigmaH(obs, comb, sigma.h.dom, pi_theta_qz, 
                           verbose = verbose)
    mode.tab[i, "sigma.h.hat"] = res[["sigma.h.hat"]]
    mode.tab[i, "log.partition"] = res[["LL"]]
    mode.tab[i, "post_LL"] = res[["LL"]]
    lambda.qz.res = res[["lambda.qz.res"]]
    mode.tab[i, "theta.z.hat"] = res$theta_qz_hat[(Q + 1)]
    res = ABSOLUTE:::GetSegQzPost(obs[["d.tx"]], obs[["d.stderr"]], 
                       obs[["W"]], comb, NA, mode.tab[i, "sigma.h.hat"], 
                       theta.qz = res[["theta_qz_hat"]])
    seg.z.tab[i, ] = res[["QZ"]][, (Q + 1)]
    seg.qz.tab[i, , ] = res[["QZ"]]
    seg.q.tab[i, , ] = res[["Q"]]
    theta.q.tab[i, ] = colSums(res[["Q"]] * obs[["W"]])
    mode.tab[i, "theta.z.hat"] = sum(seg.qz.tab[i, , (Q + 
                                                        1)] * obs$W)
    mode.tab[i, "frac.het"] = sum(obs[["W"]] * seg.z.tab[i, 
                                                         ])
    if (obs[["data.type"]] == "ALLELIC") {
      ab.tab[i, ] = CalcAbDistr(obs, res[["QZ"]])
      mode.tab[i, "genome mass"] = 2 * sum(c((1:Q) - 1) * 
                                             colSums(res[["Q"]] * obs[["W"]]))
    }
    if (obs[["data.type"]] == "TOTAL") {
      mode.tab[i, "genome mass"] = 1 * sum(c((1:Q) - 1) * 
                                             colSums(res[["Q"]] * obs[["W"]]))
    }
    mode.tab[i, "entropy"] = CalcFitEntropy(obs, res[["QZ"]])
    chr.arm.tab[i, , , ] = CalcChrArmDistr(seg.obj, res[["Q"]], 
                                           chr.arms.dat)
  }
  return(list(mode.tab = mode.tab, mode.posts = NA, theta.q.tab = theta.q.tab, 
              b = theta.qz.hat, seg.z.tab = seg.z.tab, seg.qz.tab = seg.qz.tab, 
              seg.q.tab = seg.q.tab, ab.tab = ab.tab, chr.arm.tab = chr.arm.tab, 
              mode.flag = NA))
}
<environment: namespace:ABSOLUTE>
  
  >
  
  
  
  ABSOLUTE:::FindLocationModes
function (obs, mut, Q, theta.qz, sigma.h, tau.dom, verbose = FALSE) 
{
  kAlphaDom <- c(0, 1)
  kDom1 <- c(0, 1)
  kDom2 <- log(c(0.08, 1.05))
  lambda.qz.res <- GetLambda(theta.qz, obs[["W"]], Q + 1, verbose = verbose)
  pz <- theta.qz[Q + 1]
  mode.tab <- MargModeFinder(obs, mut, kDom1, kDom2, Q, lambda.qz.res, 
                             pz, sigma.h, tau.dom, verbose = verbose)
  if (nrow(mode.tab) == 0) {
    return(list(mode.flag = "DELTA_B_DOM"))
  }
  mode.tab <- cbind(mode.tab, NA, NA, NA)
  for (i in 1:nrow(mode.tab)) {
    par <- mode.tab[i, c(1, 2)]
    if (obs[["platform"]] == "ARRAY") {
      cur.par = par
      if (is.null(obs$error.model$at)) {
        stop("at/AT issue, contact Jeff Gentry jgentry@broadinstitute.org")
      }
      mode.params = list(b = par[1], delta = par[2], AT = obs$error.model$at)
      obs$error.model$fit.at = obs$error.model$at
    }
    else {
      cur.par <- par
      mode.params <- list(b = par[1], delta = par[2], AT = NA)
    }
    LL <- CombLL(cur.par, Q = Q, obs = obs, dom1 = kDom1, 
                 dom2 = kDom2, lambda.qz.res, pz, sigma.h)
    mode.hess <- hessian(CombLL, x = cur.par, method = "Richardson", 
                         Q = Q, obs = obs, dom1 = kDom1, dom2 = kDom2, lambda.qz.res = lambda.qz.res, 
                         pz = pz, sigma.h = sigma.h)
    if (!is.na(mode.hess[1])) {
      mode.curv <- CalcModeLogCurv(par, mode.hess, verbose = verbose)
    }
    else {
      LL <- NA
      mode.curv <- NA
    }
    mode.tab[i, ] <- c(mode.params[["b"]], mode.params[["delta"]], 
                       mode.params[["AT"]], LL, mode.curv)
  }
  b <- mode.tab[, 1]
  delta <- exp(mode.tab[, 2])
  at <- mode.tab[, 3]
  res <- GetAlphaAndTau(b, delta)
  alpha <- res[["alpha"]]
  tau <- res[["tau"]]
  LL <- mode.tab[, 4]
  mode.curv <- mode.tab[, 5]
  mode.tab <- cbind(alpha, tau, at, b, delta, LL, mode.curv)
  colnames(mode.tab) <- c("alpha", "tau", "AT", "b", "delta", 
                          "LL", "mode_curv")
  if (verbose) {
    print(mode.tab)
  }
  mode.ix <- (mode.tab[, "alpha"] >= kAlphaDom[1] & mode.tab[, 
                                                             "alpha"] <= kAlphaDom[2] & mode.tab[, "tau"] >= tau.dom[1] & 
                mode.tab[, "tau"] <= tau.dom[2])
  if (verbose) {
    print(paste("removing ", sum(!mode.ix), " / ", length(mode.ix), 
                " modes outside of alpha/tau range.", sep = ""))
  }
  mode.tab <- mode.tab[mode.ix, , drop = FALSE]
  if (nrow(mode.tab) == 0) {
    return(list(mode.flag = "ALPHA_TAU_DOM"))
  }
  return(list(mode.tab = mode.tab))
}
<environment: namespace:ABSOLUTE>
  
 
  ABSOLUTE:::MargModeFinder
function (obs, mut, dom1, dom2, Q, lambda.qz.res, pz, sigma.h, 
          tau_dom, b.res = 0.125, d.res = 0.125, verbose = FALSE) 
{
  b.grid <- seq(dom1[1], dom1[2], b.res)
  d.grid <- seq(dom2[1], dom2[2], d.res)
  n.b <- length(b.grid)
  n.d <- length(d.grid)
  mode.tab <- array(NA, dim = c(n.b * n.d, 3))
  for (i in seq_len(n.b)) {
    for (j in seq_len(n.d)) {
      cur.par <- c(b.grid[i], d.grid[j])
      res <- RunOpt(cur.par, obs, dom1, dom2, lambda.qz.res, 
                    pz, sigma.h, Q, verbose = verbose)
      if (!is.na(res)) {
        mode.tab[(i - 1) * n.d + j, ] <- c(res[[1]], 
                                           res[[2]], res[[3]])
      }
    }
    if (verbose) {
      cat("\n")
    }
  }
  delta_dom = log(c(1/tau_dom[2] - 0.05, 1))
  res_1d = run_1d_opt(obs, delta_dom, d.res, lambda.qz.res, 
                      pz, sigma.h, Q, verbose = verbose)
  if (verbose) {
    print("1d mode opt: ")
    print(res_1d)
  }
  mode.tab = rbind(mode.tab, res_1d)
  if (!is.na(mut)) {
    alpha_dom = c(0.1, 1)
    res_snv_only = run_diploid_snv_purity_opt(obs, mut$mut_CN_dat, 
                                              mut$Pi_som_theta_Q, mut$mut_class_W, alpha_dom, verbose = verbose)
    if (!is.na(res_snv_only)) {
      if (verbose) {
        print("SNV_only opt: ")
        print(res_snv_only)
      }
      mode.tab = rbind(mode.tab, res_snv_only)
    }
  }
  ix <- !is.na(mode.tab[, 1])
  mode.list <- mode.tab[ix, c(1, 2), drop = FALSE]
  umodes = unique(round(mode.list, 2))
  if (verbose) {
    print(paste(nrow(umodes), " unique modes found", sep = ""))
  }
  umodes <- umodes[umodes[, 1] >= dom1[1] & umodes[, 2] >= 
                     dom2[1], , drop = FALSE]
  if (verbose) {
    print(paste(nrow(umodes), " modes in b / delta range.", 
                sep = ""))
  }
  return(umodes)
}
<environment: namespace:ABSOLUTE>
  
  
  
  
  ExtractReviewedResults
function (reviewed.pp.calls.fn, analyst.id, modes.fn, out.dir.base, 
          obj.name, copy_num_type, verbose = FALSE) 
{
  if (copy_num_type == "total") {
    set_total_funcs()
  }
  else if (copy_num_type == "allelic") {
    ABSOLUTE:::set_allelic_funcs()
  }
  else {
    stop("Unsupported copy number type: ", copy_num_type)
  }
  load(modes.fn)
  dat = read.delim(reviewed.pp.calls.fn, row.names = NULL, 
                   stringsAsFactors = FALSE, header = 1, check.names = FALSE)
  if (colnames(dat)[3] == "sample") {
    call_override = dat[, 1]
    pp.calls = dat[, c(3:ncol(dat))]
    rownames(pp.calls) = dat[, "sample"]
    names(call_override) = dat[, "sample"]
    found = intersect(names(segobj.list), rownames(pp.calls))
    segobj.list = segobj.list[found]
    call_override = call_override[found]
    called.segobj.list = override_absolute_calls(segobj.list, 
                                                 call_override)
  }
  else {
    if (colnames(dat)[2] != "sample") {
      stop("Invalid reviewed.pp.calls.fn!")
    }
    pp.calls = dat[, c(2:ncol(dat))]
    rownames(pp.calls) = dat[, "sample"]
    found = intersect(names(segobj.list), rownames(pp.calls))
    segobj.list = segobj.list[found]
    mode.ix = ABSOLUTE:::MatchPpModes(pp.calls, segobj_list = segobj.list)[["matched"]]
    called.segobj.list = ABSOLUTE:::ReduceSegobjListToCalled(mode.ix, 
                                                  segobj_list = segobj.list)
  }
  process_extract_reviewed_results(out.dir.base, called.segobj.list, 
                                   pp.calls, obj.name, analyst.id)
}
<environment: namespace:ABSOLUTE>
  
  
  ABSOLUTE:::MatchPpModes 
  function (pp.calls, segobj_list) 
  {
    N = length(segobj_list)
    mode.ix = rep(NA, N)
    wrong.ix = c()
    names(mode.ix) = names(segobj_list)
    for (i in seq_len(N)) {
      res = MatchPpMode(pp.calls, sid = names(segobj_list)[i], seg.dat = segobj_list[[i]])
      mode.ix[i] = res[1]
      if (!is.na(res) && res[2]) {
        wrong.ix = c(wrong.ix, i)
      }
    }
    return(list(matched = mode.ix, wrong.ix = wrong.ix))
  }
<environment: namespace:ABSOLUTE>
  
  ABSOLUTE:::MatchPpMode
function (pp.calls, sid, seg.dat) 
{
  if ((any(is.na(pp.calls[sid, c("purity", "ploidy")]))) || 
        (!is.na(seg.dat[["mode.res"]][["mode.flag"]]))) {
    return(NA)
  }
  mode.tab = seg.dat[["mode.res"]][["mode.tab"]]
  vals = mode.tab[, c("alpha", "genome mass"), drop = FALSE]
  call.vals = pp.calls[sid, c("purity", "ploidy")]
  call.vals = matrix(as.numeric(call.vals), ncol = 2, nrow = nrow(mode.tab), 
                     byrow = TRUE)
  min.mode = which.min(rowSums((vals - call.vals)^2))
  if ((abs(pp.calls[sid, "ploidy"] - vals[min.mode, "genome mass"]) > 
         1) || (abs(pp.calls[sid, "purity"] - vals[min.mode, "alpha"]) > 
                  0.1)) {
    return(NA)
  }
  mode.ix = min.mode
  wrong = TRUE
  return(c(mode.ix, wrong))
}
<environment: namespace:ABSOLUTE>
  

  
  ABSOLUTE:::ReduceSegobjListToCalled
function (mode.ix, segobj_list) 
{
  segobj_list = segobj_list[!is.na(mode.ix)]
  mode.ix = mode.ix[!is.na(mode.ix)]
  for (i in seq_along(segobj_list)) {
    segobj_list[[i]][["mode.res"]] = ReorderModeRes(mode.res = segobj_list[[i]][["mode.res"]], 
                                                    ix = mode.ix[i])
  }
  return(segobj_list)
}
<environment: namespace:ABSOLUTE>
  
  >
  > ABSOLUTE:::ReorderModeRes
function (mode.res, ix, DROP = FALSE) 
{
  mode.res[["mode.tab"]] = mode.res[["mode.tab"]][ix, , drop = DROP]
  mode.res[["seg.z.tab"]] = mode.res[["seg.z.tab"]][ix, , drop = DROP]
  mode.res[["seg.qz.tab"]] = mode.res[["seg.qz.tab"]][ix, , 
                                                      , drop = DROP]
  mode.res[["seg.q.tab"]] = mode.res[["seg.q.tab"]][ix, , , 
                                                    drop = DROP]
  if (!is.null(mode.res[["ab.tab"]])) {
    mode.res[["ab.tab"]] = mode.res[["ab.tab"]][ix, , drop = DROP]
  }
  mode.res[["theta.q.tab"]] = mode.res[["theta.q.tab"]][ix, 
                                                        , drop = DROP]
  mode.res[["theta.qz.hat"]] = mode.res[["theta.qz.hat"]][ix, 
                                                          , drop = DROP]
  mode.res[["chr.arm.tab"]] = mode.res[["chr.arm.tab"]][ix, 
                                                        , , , drop = DROP]
  mode.res[["mode.clust.p"]] = mode.res[["mode.clust.p"]][ix, 
                                                          , drop = DROP]
  mode.res$subclonal_SCNA_res$subclonal_SCNA_tab = mode.res$subclonal_SCNA_res$subclonal_SCNA_tab[ix, 
                                                                                                  , , drop = DROP]
  mode.res$subclonal_SCNA_res$CCF_dens = mode.res$subclonal_SCNA_res$CCF_dens[ix, 
                                                                              , , drop = DROP]
  if (!is.null(mode.res[["muts.post.prs"]])) {
    mode.res[["muts.post.prs"]] = mode.res[["muts.post.prs"]][, 
                                                              , ix, drop = DROP]
  }
  mode.res[["mode.posts"]] = mode.res[["mode.posts"]][ix]
  return(mode.res)
}
<environment: namespace:ABSOLUTE>
  
  
  ABSOLUTE:::process_extract_reviewed_results
function (out.dir.base, called.segobj.list, pp.calls, obj.name, 
          analyst.id) 
{
  seg.maf.dir = file.path(out.dir.base, "reviewed", "SEG_MAF")
  dir.create(seg.maf.dir, recursive = TRUE)
  write_called_seg_maf(called_segobj_list = called.segobj.list, pp_calls = pp.calls, 
                       out_dir = seg.maf.dir)
  out.fn = file.path(out.dir.base, "reviewed", paste(obj.name, 
                                                     ".", analyst.id, ".ABSOLUTE.table.txt", sep = ""))
  PrintPpCallTable(called.segobj.list, out.fn)
  pdf.fn = file.path(out.dir.base, "reviewed", paste(obj.name, 
                                                     ".called.ABSOLUTE.plots.pdf", sep = ""))
  PlotModes(segobj.list = called.segobj.list, pdf.fn, n.print = 1)
  indv.called.dir = file.path(out.dir.base, "reviewed", "samples")
  dir.create(indv.called.dir, recursive = TRUE)
  file.base = file.path(paste(names(called.segobj.list), ".ABSOLUTE.", 
                              analyst.id, ".called", sep = ""))
  called.files = file.path(indv.called.dir, paste(file.base, 
                                                  "RData", sep = "."))
  for (i in seq_along(called.files)) {
    seg.obj = called.segobj.list[[i]]
    save(seg.obj, file = called.files[i])
  }
}
<environment: namespace:ABSOLUTE>
  
  >

  
  
  
ABSOLUTE:::write_called_seg_maf
function (called_segobj_list, pp_calls, out_dir, verbose = FALSE) 
{
  dir.create(out_dir)
  for (i in seq_along(called_segobj_list)) {
    called_segobj = called_segobj_list[[i]]
    s_name = called_segobj$array.name
    at = called_segobj$mode.res$mode.tab[1, "AT"]
    abs_seg = GetAbsSegDat(segobj = called_segobj, at)
    WriteAbsSegtab(list(abs_seg), s_name, file.path(out_dir, 
                                                    paste(s_name, "segtab.txt", sep = ".")))
    if (!is.null(called_segobj$mode.res$muts.post.prs)) {
      maf_out_fn = file.path(out_dir, paste(s_name, "_ABS_MAF.txt", 
                                            sep = ""))
      modeled = called_segobj$mode.res$muts.post.prs[, 
                                                     , 1, drop = TRUE]
      if (is.null(dim(modeled))) {
        modeled = matrix(modeled, nrow = 1)
        colnames(modeled) = dimnames(called_segobj$mode.res$muts.post.prs)[[2]]
      }
      mut_dat = cbind(called_segobj$mut.cn.dat, modeled)
      ccf_grid = seq(0, 1, by = 0.01)
      alpha = mut_dat[1, "purity"]
      ccf_dens = matrix(NA, nrow = nrow(mut_dat), ncol = length(ccf_grid))
      for (i in seq_along(nrow(mut_dat))) {
        if (mut_dat[i, "q_hat"] == 0) {
          next
        }
        ccf_dens[i, ] = calc_ccf_posterior_grid(mut_dat[i, 
                                                        "alt"], mut_dat[i, "ref"], alpha, mut_dat[i, 
                                                                                                  "q_hat"], ccf_grid)
      }
      colnames(ccf_dens) = paste("CCF_", ccf_grid, sep = "")
      out_mut_dat = cbind(mut_dat, ccf_dens)
      write.table(file = maf_out_fn, out_mut_dat, row.names = FALSE, 
                  sep = "\t", quote = FALSE)
    }
    if (verbose) {
      cat(".")
    }
  }
}
<environment: namespace:ABSOLUTE>
  
  >
  
  GetAbsSegDat
function (segobj, at) 
{
  Q <- 15
  qq <- Q
  seg.qz.tab <- segobj[["mode.res"]][["seg.qz.tab"]][1, , ]
  seg.q.tab <- segobj[["mode.res"]][["seg.q.tab"]][1, , ]
  max.mat <- apply(seg.qz.tab, 1, which.max)
  subclonal.ix <- max.mat == (Q + 1)
  max.mat <- apply(seg.q.tab, 1, which.max)
  exp.mat <- apply(seg.q.tab, 1, function(x) {
    x <- x[1:qq]/sum(x[1:qq])
    return(sum(x * c(1:qq)))
  })
  seg.list <- segobj[["as.seg.dat"]]
  u <- unique(seg.list[, "seg.ix"])
  modal.a1 <- vector(mode = "numeric", length = length(u))
  expected.a1 <- vector(mode = "numeric", length = length(u))
  modal.a2 <- vector(mode = "numeric", length = length(u))
  expected.a2 <- vector(mode = "numeric", length = length(u))
  LOH <- vector(mode = "numeric", length = length(u))
  HZ <- vector(mode = "numeric", length = length(u))
  subclonal.a1 <- vector(mode = "numeric", length = length(u))
  subclonal.a2 <- vector(mode = "numeric", length = length(u))
  copy.ratio <- vector(mode = "numeric", length = length(u))
  hscr.a1 <- vector(mode = "numeric", length = length(u))
  hscr.a2 <- vector(mode = "numeric", length = length(u))
  cancer.cell.frac.a1 <- vector(mode = "numeric", length = length(u))
  cancer.cell.frac.a2 <- vector(mode = "numeric", length = length(u))
  ccf.ci95.low.a1 <- vector(mode = "numeric", length = length(u))
  ccf.ci95.low.a2 <- vector(mode = "numeric", length = length(u))
  ccf.ci95.high.a1 <- vector(mode = "numeric", length = length(u))
  ccf.ci95.high.a2 <- vector(mode = "numeric", length = length(u))
  sc_tab = segobj$mode.res$subclonal_SCNA_res$subclonal_SCNA_tab[1, 
                                                                 , ]
  ccf_hat = round(sc_tab[, "CCF_hat"], 5)
  ccf_ci95 = round(sc_tab[, c("CI95_low", "CI95_high")], 5)
  for (i in seq_along(u)) {
    seg <- u[i]
    usegs <- which(seg.list[, "seg.ix"] == seg)
    hscr <- InvAtten(sort(seg.list[usegs, "copy_num"]), at)
    copy.ratio[i] <- round(sum(hscr)/2, 5)
    hscr.a1[i] <- round(hscr[1], 5)
    hscr.a2[i] <- round(hscr[2], 5)
    modal.a1[i] <- max.mat[usegs[1]] - 1
    modal.a2[i] <- max.mat[usegs[2]] - 1
    expected.a1[i] <- round(exp.mat[usegs[1]] - 1, 5)
    expected.a2[i] <- round(exp.mat[usegs[2]] - 1, 5)
    subclonal.a1[i] <- subclonal.ix[usegs[1]]
    subclonal.a2[i] <- subclonal.ix[usegs[2]]
    cancer.cell.frac.a1[i] = ccf_hat[usegs[1]]
    cancer.cell.frac.a2[i] = ccf_hat[usegs[2]]
    ccf.ci95.low.a1[i] = ccf_ci95[usegs[1], 1]
    ccf.ci95.high.a1[i] = ccf_ci95[usegs[1], 2]
    ccf.ci95.low.a2[i] = ccf_ci95[usegs[2], 1]
    ccf.ci95.high.a2[i] = ccf_ci95[usegs[2], 2]
    if ((modal.a1[i] == 0) && (modal.a2[i] == 0)) {
      HZ[i] <- 1
    }
    else if ((modal.a1[i] == 0) || (modal.a2[i] == 0)) {
      LOH[i] <- 1
    }
  }
  ix <- which(colnames(seg.list) %in% c("seg.ix", "copy.num"))
  tab <- round(seg.list[1:length(u), c(-ix)], 5)
  return(cbind(tab, copy.ratio, hscr.a1, hscr.a2, modal.a1, 
               modal.a2, expected.a1, expected.a2, subclonal.a1, subclonal.a2, 
               cancer.cell.frac.a1, ccf.ci95.low.a1, ccf.ci95.high.a1, 
               cancer.cell.frac.a2, ccf.ci95.low.a2, ccf.ci95.high.a2, 
               LOH, HZ))
} 
  

ABSOLUTE:::override_absolute_calls
function (segobj.list, call_override, verbose = FALSE) 
{
  mode_str = call_override
  mode.ix = rep(NA, length(mode_str))
  status.vec = rep(NA, length(mode_str))
  empty_ix = (call_override == "") | (is.na(call_override))
  call_override[empty_ix] = "1"
  called_ix = !is.na(as.integer(call_override))
  status.vec[called_ix] = "called"
  mode.ix[called_ix] = as.integer(call_override[called_ix])
  mode.ix[!called_ix] = 1
  status.vec[!called_ix] = call_override[!called_ix]
  segobj.list = SelectMatchedModes(segobj_list = segobj.list, mode.ix, status.vec, 
                                   verbose = verbose)
  return(segobj.list)
}
<environment: namespace:ABSOLUTE>
  
  ABSOLUTE:::SelectMatchedModes
function (segobj_list, mode.ix, status.vec, verbose = FALSE) 
{
  new.segobj.list = segobj_list[!is.na(mode.ix)]
  mode.ix = mode.ix[!is.na(mode.ix)]
  for (i in seq_along(mode.ix)) {
    segobj = new.segobj.list[[i]]
    segobj[["mode.res"]][["call.status"]] = status.vec[i]
    if (!(status.vec[i] %in% c("non-clonal", "non-aneuploid", 
                               "low purity", "FAILED"))) {
      mode.res = segobj[["mode.res"]]
      ix = mode.ix[i]
      segobj[["mode.res"]] = ReorderModeRes(mode.res, ix)
    }
    new.segobj.list[[i]] = segobj
    if (verbose) {
      cat(".")
    }
  }
  if (verbose) {
    cat("\n")
  }
  return(new.segobj.list)
} 
  


ABSOLUTE:::PostPlot
function (mode.tab, mode.colors, alpha.dom, tau.dom, sample.name, 
          obs, called.mode.ix, debug.info, call.status = NA, model.id = NA) 
{
  n.print <- nrow(mode.tab)
  mode.tab <- mode.tab[c(1:n.print), , drop = FALSE]
  mode.colors <- mode.colors[c(1:NROW(mode.tab))]
  a1 <- mode.tab[1, "alpha"]
  t1 <- mode.tab[1, "tau"]
  den <- 2 * (1 - a1) + a1 * t1
  delta <- a1/den
  b <- 2 * (1 - a1)/den
  delta.set <- c(delta/2, delta, 2 * delta)
  b.set <- c(b, b - 2 * delta, b + 2 * delta)
  b.set <- b.set[b.set >= 0 & b.set < 1]
  old.mar <- par("mar")
  if (debug.info) {
    par(mar = old.mar + c(0, 0, 1, 0))
  }
  ylab <- substitute(paste("Fraction cancer nuclei", (hat(alpha)), 
                           sep = ""))
  xlab <- substitute(paste("Ploidy", (hat(tau)), sep = ""))
  plot(0, type = "n", ylab = ylab, xlab = xlab, xlim = tau.dom, 
       ylim = alpha.dom, main = "", bty = "n", las = 1, cex = 1.5, 
       xaxt = "n")
  at <- seq(2, max(tau.dom), by = 2)
  axis(side = 1, at = at, labels = paste(at, "N", sep = ""))
  for (i in seq_along(b.set)) {
    curve(ABSOLUTE:::AlphaBFunc(x, b.set[i]), from = tau.dom[1], to = tau.dom[2], 
          col = "grey30", add = TRUE, lty = 3)
  }
  for (i in seq_along(delta.set)) {
    curve(ABSOLUTE:::AlphaFunc(x, delta.set[i]), from = tau.dom[1], 
          to = tau.dom[2], col = "grey30", add = TRUE, lty = 3, 
          n = 500)
  }
  points(mode.tab[, "tau"], mode.tab[, "alpha"], col = "black", 
         bg = mode.colors, pch = 21, cex = par("cex") * 2)
  if (!is.na(called.mode.ix)) {
    ix <- called.mode.ix
    points(mode.tab[ix, "tau"], mode.tab[ix, "alpha"], col = "red", 
           pch = "*", cex = par("cex") * 4)
  }
  if (debug.info) {
    title(sample.name, line = 3)
    mtext(call.status, line = 0, adj = 0)
    mtext(model.id, line = 0, adj = 1)
  }
  if ((obs[["platform"]] == "ARRAY") && debug.info) {
    hscn.params <- obs[["error.model"]]
    mtext(substitute(paste(sigma[eta] == x, ",", sigma[epsilon] == 
                             y, sep = ""), list(x = round(hscn.params[["sigma.eta"]], 
                                                          3), y = round(hscn.params[["sigma.nu"]], 3))), line = 1.5, 
          adj = 0)
  }
  par(mar = old.mar)
}



ABSOLUTE:::PpModeScorerBarplot = 
function (mode.tab, mode.colors, obs, n.print) 
{
  n.print <- min(n.print, nrow(mode.tab))
  mode.colors <- mode.colors[c(1:NROW(mode.tab))]
  if (("somatic.mut.ll" %in% colnames(mode.tab)) && (!any(is.na(mode.tab[1:n.print, 
                                                                         "somatic.mut.ll"])))) {
    mat <- mode.tab[1:n.print, c("log.partition", "clust_LL", 
                                 "somatic.mut.ll"), drop = FALSE]
    colnames(mat) <- c("SCNAs", "karyotype", "SSNVs")
  }
  else {
    mat <- mode.tab[c(1:n.print), c("log.partition", "clust_LL"), 
                    drop = FALSE]
    colnames(mat) <- c("SCNAs", "karyotype")
  }
  ix <- apply(!is.finite(mat), 1, sum) > 0
  mat <- mat[!ix, , drop = FALSE]
  mat <- t(t(mat) - apply(mat, 2, min))
  mat <- mat + max(mat) * 0.1
  mat <- cbind(mat, rowMeans(mat))
  colnames(mat)[ncol(mat)] <- "combined"
  barplot(mat, beside = TRUE, col = mode.colors[1:n.print[!ix]], 
          axes = FALSE, ylab = "", space = c(0, 2), cex.names = par("cex.axis"))
  mtext("Log-likelihood", side = 2, line = 1, las = 3, cex = par("cex") * 
          par("cex.axis"))
  mtext("Model-based evaluation", side = 3, line = 0, adj = 0, 
        cex = par("cex.axis"))
  axis(side = 2, labels = FALSE)
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

GetScnaStderrGridDensity
function (obs, grid, sigma.h, i = NA) 
{
  sigma.nu <- obs[["error.model"]][["sigma.nu"]]
  sigma.eta <- obs[["error.model"]][["sigma.eta"]]
  b <- (sqrt(exp(sigma.eta^2) - 1))/sigma.nu
  a <- 0
  N <- length(obs[["d.tx"]])
  grid.dens <- matrix(NA, nrow = N, ncol = length(grid))
  dx <- c(0, diff(grid))
  for (i in seq_len(N)) {
    grid.dens[i, ] <- dnorm(ABSOLUTE:::HTx(grid, a, b), obs[["d.tx"]][i], 
                            obs[["d.stderr"]][i] + sigma.h) * ABSOLUTE:::HTxVar(grid, a, 
                                                                     b) * dx
  }
  return(grid.dens)
}


ABSOLUTE:::GetLambda
function (theta, W, n.states, lambda.init = NA, w.thresh = 0, 
          verbose = FALSE) 
{
  if (all(theta == theta[1])) {
    lambda.res <- list()
    lambda.res[["Lambda"]] <- rep(0, (n.states - 1))
    lambda.res[["norm.term"]] <- rep(1/n.states, length(W))
    lambda.res[["resid"]] <- 0
    return(lambda.res)
  }
  if (!is.na(lambda.init)) {
    lambda.res <- CalcLambda(theta, W, n.states, first.lambda.guess = lambda.init, 
                             verbose = verbose)
    resid <- lambda.res[["resid"]]
    if (resid < 0.05) {
      return(lambda.res)
    }
  }
  wsix <- sort(W, index.return = TRUE)[["ix"]]
  orig.w <- W
  if (w.thresh > 0) {
    n.small <- max(which(cumsum(W[wsix]) < w.thresh))
    small.ix <- wsix[c(1:n.small)]
    W <- W[-small.ix]
    W <- W/sum(W)
  }
  theta[theta < .Machine$double.xmin] = .Machine$double.xmin
  resid <- Inf
  wix <- sort(W, decreasing = TRUE, index.return = TRUE)[["ix"]]
  iix <- 1
  init.ix <- wix[iix]
  while (resid > 0.05 & iix <= length(W)) {
    lambda.init <- -1/W[init.ix] * log(theta[c(2:length(theta))]/theta[1])/c(1:(length(theta) - 
                                                                                  1))
    lambda.init[lambda.init == Inf] = .Machine$double.xmax/2
    lambda.init[lambda.init == -Inf] = -.Machine$double.xmax/2
    lambda.res <- CalcLambda(theta, W, n.states, first.lambda.guess = lambda.init, 
                             verbose = verbose)
    resid <- lambda.res[["resid"]]
    iix <- iix + 1
    init.ix <- wix[iix]
  }
  if (resid > 0.05) {
    if (verbose) {
      print(paste("Warning: Lambda resid = ", resid, sep = ""))
    }
    stop()
  }
  if (w.thresh > 0) {
    norm.term <- rep(NA, length(orig.w))
    norm.term[small.ix] <- 1/n.states
    norm.term[-small.ix] <- lambda.res[["norm.term"]]
    lambda.res[["norm.term"]] <- norm.term
  }
  return(lambda.res)
}
<environment: namespace:ABSOLUTE>
  
  >
  
  
 
  ABSOLUTE:::CalcNormLoglik
function (d, sigma, w, comb, lambda.qz.res, pz, sigma.h, comb.s = 1) 
{
  if (is.list(lambda.qz.res)) {
    log.pqz <- t((ABSOLUTE:::LambdaExpr(lambda.qz.res[["Lambda"]], w))) + 
      log(lambda.qz.res[["norm.term"]])
    q <- length(lambda.qz.res[["Lambda"]])
    comb <- comb[1:q]
  } else {
    q <- length(comb)
    log.pqz <- matrix(log(pz), nrow = length(w), ncol = q + 
                        1, byrow = TRUE)
  }
  log.ss.p.1 <- sapply(comb, dnorm, d, sqrt(sigma^2 + comb.s^2 * 
                                              sigma.h^2), log = TRUE)
  log.comb.seg.mat <- cbind(log.pqz[, c(1:q)] + log.ss.p.1, 
                            log(1/7) + log.pqz[, q + 1])
  LL <- sum(LogAdd(log.comb.seg.mat))
  return(LL)
}
<environment: namespace:ABSOLUTE>
  
  
  ABSOLUTE:::AbsoluteResultPlot
function (sample.pdf.fn, seg.dat, called.mode.ix = NA, verbose = FALSE) 
{
  pdf(sample.pdf.fn, 14, 15)
  d.mar <- par("mar")
  layout(mat = matrix(1:4, nrow = 2, ncol = 2, byrow = TRUE), 
         widths = c(6, 2), heights = c(6, 6))
  par(las = 1)
  if (!is.null(seg.dat[["allele.segs"]])) {
    PlotHscrAndSeghist(seg.dat, called.mode.ix)
    frame()
    frame()
  }
  par(mar = d.mar)
  par(mfrow = c(4, 4))
  PlotModes(list(seg.dat), sample.pdf.fn, n.print = 19, debug.info = TRUE, 
            add = TRUE, fig.mode = FALSE, sideways = FALSE, called.mode.ix = called.mode.ix, 
            verbose = verbose)
  dev.off()
}
<environment: namespace:ABSOLUTE>

