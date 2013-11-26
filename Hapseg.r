  
myRunHapSeg = function (plate.name, array.name, seg.fn, snp.fn, genome.build, 
          results.dir, platform, use.pop, impute.gt, plot.segfit, merge.small, 
          merge.close, min.seg.size, normal, out.p, seg.merge.thresh, 
          use.normal, adj.atten, phased.bgl.dir, force.diploid = normal, 
          drop.x = FALSE, drop.y = TRUE, calls.fn = NULL, mn.sample = NULL, 
          out.file = NULL, calibrate.data = NULL, clusters.fn = NULL, 
          prev.theta.fn = NULL, snp.file.parser = DefaultSnpFileParser, 
          clusters.file.parser = DefaultClustersFileParser, verbose = FALSE) 
{
  if (!HAPSEG:::CheckGenomeBuild(genome.build)) {
    stop("Unsupported genome build: ", genome.build)
  }
  if ((!impute.gt) && (verbose)) {
    print("Note: impute.gt is set to FALSE. It is recommended that this be set to TRUE")
  }
  if (!identical(adj.atten, FALSE)) {
    print("Note: The adj.atten argument has been disabled and will be removed in a future release")
  }
  tmp.dir = HAPSEG:::CreateTmpDir(results.dir)
  on.exit(try(unlink(tmp.dir, recursive = TRUE), silent = TRUE), 
          add = TRUE)
  iams.res = HAPSEG:::InitAndMergeSmall(array.name, genome.build, use.pop, 
                               use.normal, normal, impute.gt, platform, seg.fn, snp.fn, 
                               drop.x, drop.y, calls.fn, mn.sample, min.seg.size, merge.close, 
                               out.p, merge.small, force.diploid, calibrate.data, clusters.fn, 
                               prev.theta.fn, snp.file.parser = snp.file.parser, clusters.file.parser = clusters.file.parser, 
                               verbose = verbose)
  if (merge.close == TRUE) {
    mcaf.res = HAPSEG:::MergeCloseAndFit(iams.res, out.p, seg.merge.thresh, 
                                impute.gt, platform, force.diploid = force.diploid, 
                                verbose = verbose)
  }
  else {
    if (verbose) {
      print("Merged Loci: nulled ")
    }
    seg.dat[["merged.loci"]] = NA
    seg.dat[["final.merge.prob"]] = matrix(NA, ncol = 2, 
                                           nrow = length(iams.res[["h.d"]]) - 1)
    seg.dat[["seg.expected.phase"]] = iams.res[["em.fit"]][["seg.expected.phase"]]
    mcaf.res = iams.res
    mcaf.res[["seg.dat"]] = seg.dat
  }
  if (impute.gt == TRUE) {
    gt.res = HAPSEG:::ImputeGt(mcaf.res, platform, tmp.dir, plate.name, 
                      out.p, phased.bgl.dir, verbose = verbose)
    gt.res = HAPSEG:::PostImputeSegFit(gt.res, out.p, platform, force.diploid = force.diploid, 
                              verbose = verbose)
  }
  else {
    gt.res = mcaf.res
  }
  seg.dat = HAPSEG:::FinishHapsegAndSave(gt.res, plate.name, normal, 
                                array.name, results.dir, platform, out.file)
  if (plot.segfit) {
    if (nrow(seg.dat$merged.loci) == 0) seg.dat$merged.loci = NA
    HAPSEG:::DoPlots(seg.dat, results.dir, platform, verbose = verbose)
  }
  return(TRUE)
}

HAPSEG:::FinishHapsegAndSave
function (seg.obj, plate.name, normal, array.name, results.dir, 
          platform, out.file) 
{
  SetPlatformSpecificFuncs(platform)
  h.d <- seg.obj[["h.d"]]
  h.snp.annot <- seg.obj[["h.snp.annot"]]
  em.fit <- seg.obj[["em.fit"]]
  h.snp.gt.p <- seg.obj[["h.snp.gt.p"]]
  seg.dat <- seg.obj[["seg.dat"]]
  as.seg.mat <- ConstructAsSegMat(h.d, h.snp.annot, em.fit)
  em.fit[["h.het.prob"]] <- h.snp.gt.p
  em.fit[["h.snp.annot"]] <- h.snp.annot
  seg.dat <- HAPSEG:::BuildFinalSegDat(seg.dat, as.seg.mat, em.fit, 
                              h.d, plate.name, normal, array.name)
  if (is.null(out.file)) {
    out.file <- paste(plate.name, array.name, "segdat.RData", 
                      sep = "_")
  }
  save(seg.dat, file = file.path(results.dir, out.file))
  return(invisible(seg.dat))
}

RunHapSeg
function (plate.name, array.name, seg.fn, snp.fn, genome.build, 
          results.dir, platform, use.pop, impute.gt, plot.segfit, merge.small, 
          merge.close, min.seg.size, normal, out.p, seg.merge.thresh, 
          use.normal, adj.atten, phased.bgl.dir, force.diploid = normal, 
          drop.x = FALSE, drop.y = TRUE, calls.fn = NULL, mn.sample = NULL, 
          out.file = NULL, calibrate.data = NULL, clusters.fn = NULL, 
          prev.theta.fn = NULL, snp.file.parser = DefaultSnpFileParser, 
          clusters.file.parser = DefaultClustersFileParser, verbose = FALSE) 
{
  if (!CheckGenomeBuild(genome.build)) {
    stop("Unsupported genome build: ", genome.build)
  }
  if ((!impute.gt) && (verbose)) {
    print("Note: impute.gt is set to FALSE. It is recommended that this be set to TRUE")
  }
  if (!identical(adj.atten, FALSE)) {
    print("Note: The adj.atten argument has been disabled and will be removed in a future release")
  }
  tmp.dir = CreateTmpDir(results.dir)
  on.exit(try(unlink(tmp.dir, recursive = TRUE), silent = TRUE), 
          add = TRUE)
  iams.res = InitAndMergeSmall(array.name, genome.build, use.pop, 
                               use.normal, normal, impute.gt, platform, seg.fn, snp.fn, 
                               drop.x, drop.y, calls.fn, mn.sample, min.seg.size, merge.close, 
                               out.p, merge.small, force.diploid, calibrate.data, clusters.fn, 
                               prev.theta.fn, snp.file.parser = snp.file.parser, clusters.file.parser = clusters.file.parser, 
                               verbose = verbose)
  if (merge.close == TRUE) {
    mcaf.res = MergeCloseAndFit(iams.res, out.p, seg.merge.thresh, 
                                impute.gt, platform, force.diploid = force.diploid, 
                                verbose = verbose)
  }
  else {
    if (verbose) {
      print("Merged Loci: nulled ")
    }
    seg.dat[["merged.loci"]] = NA
    seg.dat[["final.merge.prob"]] = matrix(NA, ncol = 2, 
                                           nrow = length(iams.res[["h.d"]]) - 1)
    seg.dat[["seg.expected.phase"]] = iams.res[["em.fit"]][["seg.expected.phase"]]
    mcaf.res = iams.res
    mcaf.res[["seg.dat"]] = seg.dat
  }
  if (impute.gt == TRUE) {
    gt.res = ImputeGt(mcaf.res, platform, tmp.dir, plate.name, 
                      out.p, phased.bgl.dir, verbose = verbose)
    gt.res = PostImputeSegFit(gt.res, out.p, platform, force.diploid = force.diploid, 
                              verbose = verbose)
  }
  else {
    gt.res = mcaf.res
  }
  seg.dat = FinishHapsegAndSave(gt.res, plate.name, normal, 
                                array.name, results.dir, platform, out.file)
  if (plot.segfit) {
    DoPlots(seg.dat, results.dir, platform, verbose = verbose)
  }
  return(TRUE)
}
<environment: namespace:HAPSEG>
  
  >
  
  

HAPSEG:::DoPlots
function (seg.dat, results.dir, platform, verbose = FALSE) 
{
  if (!file.exists(results.dir)) {
    dir.create(results.dir, recursive = TRUE)
  }
  HAPSEG:::SetPlatformSpecificFuncs(platform)
  h.d <- seg.dat[["as.seg.dat"]]
  n.segs <- length(h.d)
  seg.plot.ix <- seq_along(h.d)
  theta <- seg.dat[["em.res"]][["theta"]]
  theta[["p.snp.cond.out"]] <- 1/25
  segment.chroms <- seg.dat[["allele.segs"]][, "Chromosome"]
  unique.seg.chroms <- unique(segment.chroms)
  for (i in seq_along(unique.seg.chroms)) {
    dir.create(HAPSEG:::GetChromPlotDir(unique.seg.chroms[i], results.dir))
  }
  merge.prob <- rbind(seg.dat[["final.merge.prob"]], c(NA, 
                                                       NA))
  em.fit <- seg.dat[["em.res"]]
  h.snp.annot <- em.fit[["h.snp.annot"]]
  for (i in seq_len(n.segs)) {
    segfit.fn <- file.path(HAPSEG:::GetChromPlotDir(segment.chroms[i], 
                                           results.dir), paste("HAPSEG_SEG_", i, ".jpg", sep = ""))
    jpeg(segfit.fn, 7, 5, units = "in", type = "cairo", res = 200, 
         quality = 100)
    par(cex = 0.5)
    par(cex.axis = 0.75)
    par(cex.lab = 0.75)
    par(las = 1)
    par(mfcol = c(2, 3))
    if (ncol(h.d[[i]]) > 5) {
      HAPSEG:::PlotSegFit(h.d[[i]] - theta[["bg"]], em.fit[["h.het.prob"]][[i]], 
                 em.fit[["snp.clust.p"]][[i]], h.snp.annot[[i]], 
                 em.fit[["e.mu"]][i, ], theta, seg.dat[["merged.loci"]], 
                 merge.prob[i, ], em.fit[["seg.log.ev"]][i], verbose = verbose)
    }
    dev.off()
    if (verbose) {
      cat(".")
    }
  }
}
<environment: namespace:HAPSEG>
  
  HAPSEG:::PlotSegFit
function (d, snp.gt.p, snp.clust.p, snp.annot, e.mu, theta, merged.loci, 
          merge.prob, seg.log.ev, marker.gt.pal = NA, verbose = FALSE) 
{
  data.name <- "Calibrated Intensities"
  par(bty = "n")
  if (verbose) {
    print(summary(snp.gt.p))
  }
  sigma.nu <- theta[["sigma.nu"]]
  sigma.eta <- theta[["sigma.eta"]]
  nu <- theta[["nu"]]
  het.cov <- theta[["het.cov"]]
  het.cov <- min(het.cov, (het.cov * e.mu[1] * e.mu[2]))
  mu0 <- 0
  sigma.h <- HAPSEG:::GetSigmaH(sigma.nu, sigma.eta)
  ct <- colSums(snp.gt.p)
  ct <- ct/sum(ct)
  ngt.t <- c((ct[1] + ct[3]), (ct[2] + ct[4]))
  cex <- 1.2
  min <- -0.5
  max <- 5.25
  t.min <- -0.5
  t.max <- 2.75
  het.col <- "green"
  hom.col <- "black"
  df <- d[((d >= min) & (d < max))]
  hist(df, breaks = 100, freq = FALSE, main = "", xlab = "Allelic copy-ratio", 
       xlim = c(min, max))
  xgl <- 1001
  xg <- array(NA, dim = c(4, xgl))
  x <- seq(min(df), max(df), length.out = xgl)
  x.tx <- HAPSEG:::HTx(x, sigma.nu, sigma.eta)
  e.mu <- HAPSEG:::HTx(e.mu, sigma.nu, sigma.eta)
  mu0.tx <- HAPSEG:::HTx(mu0, sigma.nu, sigma.eta)
  xg[1, ] <- ngt.t[1] * HAPSEG:::DFunc(x.tx, mu0.tx, sigma.h, nu) * 
    HAPSEG:::HTxVar(x, sigma.nu, sigma.eta)
  xg[2, ] <- ngt.t[2] * HAPSEG:::DFunc(x.tx, e.mu[1], sigma.h, nu) * 
    HAPSEG:::HTxVar(x, sigma.nu, sigma.eta)
  xg[3, ] <- ngt.t[2] * HAPSEG:::DFunc(x.tx, e.mu[2], sigma.h, nu) * 
    HAPSEG:::HTxVar(x, sigma.nu, sigma.eta)
  xg[4, ] <- ngt.t[1] * HAPSEG:::DFunc(x.tx, e.mu[3], sigma.h, nu) * 
    HAPSEG:::HTxVar(x, sigma.nu, sigma.eta)
  csum <- apply(xg, 2, sum)
  crds <- par("usr")
  scale <- 1/2
  for (i in 1:4) {
    lines(x, scale * xg[i, ], lwd = 2, col = hom.col)
  }
  lines(x, scale * csum, lwd = 2, col = "coral")
  xg[1, ] <- ngt.t[1] * HAPSEG:::DFunc(x.tx, mu0.tx, sigma.h, nu)
  xg[2, ] <- ngt.t[2] * HAPSEG:::DFunc(x.tx, e.mu[1], sigma.h, nu)
  xg[3, ] <- ngt.t[2] * HAPSEG:::DFunc(x.tx, e.mu[2], sigma.h, nu)
  xg[4, ] <- ngt.t[1] * HAPSEG:::DFunc(x.tx, e.mu[3], sigma.h, nu)
  csum <- apply(xg, 2, sum)
  d.tx <- HAPSEG:::HTx(d, sigma.nu, sigma.eta)
  d.tx.f <- d.tx[(d.tx > t.min) & (d.tx < t.max)]
  hist(d.tx.f, breaks = 100, freq = FALSE, main = "", xlab = "Allelic copy-ratio")
  mtext("Variance-stabilizing transformation", side = 3, line = 0, 
        adj = 0, cex = par("cex") * par("cex.axis"))
  for (i in 1:4) {
    lines(x.tx, scale * xg[i, ], lwd = 2, col = hom.col)
  }
  lines(x.tx, scale * csum, lwd = 2, col = "coral")
  gl <- 250
  x2d <- seq(t.min, t.max, length.out = gl)
  x1 <- rep(x2d, gl)
  x2 <- as.vector(sapply(x2d, rep, gl))
  x <- cbind(x1, x2)
  hom.cov <- 0
  sigma <- matrix(c(sigma.h^2, het.cov, het.cov, sigma.h^2), 
                  nrow = 2, ncol = 2)
  invert.sigma <- solve(sigma)
  hom.sigma <- matrix(c(sigma.h^2, 0, 0, sigma.h^2), nrow = 2, 
                      ncol = 2)
  invert.hom.sigma <- solve(hom.sigma)
  het.1 <- array(scale * ngt.t[2] * DmvFunc(x, c(e.mu[1], e.mu[2]), 
                                            sigma, invert.sigma, nu), dim = c(gl, gl))
  het.2 <- array(scale * ngt.t[2] * DmvFunc(x, c(e.mu[2], e.mu[1]), 
                                            sigma, invert.sigma, nu), dim = c(gl, gl))
  hom.1 <- array(scale * ngt.t[1] * DmvFunc(x, c(mu0.tx, e.mu[3]), 
                                            hom.sigma, invert.hom.sigma, nu), dim = c(gl, gl))
  hom.2 <- array(scale * ngt.t[1] * DmvFunc(x, c(e.mu[3], mu0.tx), 
                                            hom.sigma, invert.hom.sigma, nu), dim = c(gl, gl))
  csum2d <- het.1 + het.2 + hom.1 + hom.2
  het.md <- DmvFunc(c(e.mu[1], e.mu[2]), c(e.mu[1], e.mu[2]), 
                    sigma, invert.sigma, nu) * scale * ngt.t[2]
  hom.md <- DmvFunc(c(mu0.tx, e.mu[3]), c(mu0.tx, e.mu[3]), 
                    hom.sigma, invert.hom.sigma, nu) * scale * ngt.t[1]
  level_scales <- c(0.99, 0.95, 0.8, 0.5, 0.2)
  out.pch <- rep(".", ncol(d))
  if (any(is.na(marker.gt.pal))) {
    pal <- HAPSEG:::GetHapsegMarkerPal()
  }
  else {
    pal <- marker.gt.pal
  }
  cols <- rev(pal(1000))
  state.mat <- matrix(c(0, 1, 2, 3), ncol = 4, nrow = nrow(snp.clust.p), 
                      byrow = TRUE)
  col.ix <- round(rowSums(state.mat * snp.clust.p[, c(3, 2, 
                                                      4, 1)])/3 * 999, 0) + 1
  pcols <- cols[col.ix]
  x <- d[1, ]
  y <- d[2, ]
  ix <- (x < max) & (x > min) & (y < max) & (y > min)
  if (sum(ix) > 5) {
    plot(0, type = "n", main = "", xlim = c(min, max), ylim = c(min, 
                                                                max), xlab = "A allele copy-ratio", ylab = "B allele copy-ratio", 
         bty = "n")
  }
  else {
    smoothScatter(x, y, main = data.name, xlab = "A allele copy-ratio", 
                  ylab = "B allele copy-ratio", bty = "n")
  }
  points(d[1, ], d[2, ], col = pcols, pch = out.pch, cex = cex)
  cl.het = contourLines(x = x2d, y = x2d, z = het.1, levels = level_scales * 
                          het.md)
  cl.hom <- contourLines(x = x2d, y = x2d, z = hom.1, levels = level_scales * 
                           hom.md)
  for (i in seq_along(cl.het)) {
    txi.x <- HAPSEG:::HTxInv(cl.het[[i]][["x"]], sigma.nu, sigma.eta)
    txi.y <- HAPSEG:::HTxInv(cl.het[[i]][["y"]], sigma.nu, sigma.eta)
    lines(txi.x, txi.y, col = het.col)
    txi.x <- HAPSEG:::HTxInv(cl.het[[i]][["y"]], sigma.nu, sigma.eta)
    txi.y <- HAPSEG:::HTxInv(cl.het[[i]][["x"]], sigma.nu, sigma.eta)
    lines(txi.x, txi.y, col = het.col)
  }
  for (i in seq_along(cl.hom)) {
    txi.x <- HAPSEG:::HTxInv(cl.hom[[i]][["x"]], sigma.nu, sigma.eta)
    txi.y <- HAPSEG:::HTxInv(cl.hom[[i]][["y"]], sigma.nu, sigma.eta)
    lines(txi.x, txi.y, col = hom.col)
    txi.x <- HAPSEG:::HTxInv(cl.hom[[i]][["y"]], sigma.nu, sigma.eta)
    txi.y <- HAPSEG:::HTxInv(cl.hom[[i]][["x"]], sigma.nu, sigma.eta)
    lines(txi.x, txi.y, col = hom.col)
  }
  abline(v = mu0, lty = 2)
  abline(h = mu0, lty = 2)
  d.tx <- HAPSEG:::HTx(d, sigma.nu, sigma.eta)
  plot(0, type = "n", main = "", xlim = c(t.min, t.max), ylim = c(t.min, 
       t.max), xlab = "A allele copy-ratio", ylab = "B allele copy-ratio", 
       bty = "n")
  mtext("Variance-stabilizing transformation", side = 3, line = 0, 
        adj = 0, cex = par("cex") * par("cex.axis"))
  points(d.tx[1, ], d.tx[2, ], col = pcols, pch = out.pch, 
         cex = cex)
  contour(x = x2d, y = x2d, z = het.1, add = TRUE, drawlabels = FALSE, 
          levels = level_scales * het.md, col = het.col)
  contour(x = x2d, y = x2d, z = het.2, add = TRUE, drawlabels = FALSE, 
          levels = level_scales * het.md, col = het.col)
  contour(x = x2d, y = x2d, z = hom.1, add = TRUE, drawlabels = FALSE, 
          levels = level_scales * hom.md, col = hom.col)
  contour(x = x2d, y = x2d, z = hom.2, add = TRUE, drawlabels = FALSE, 
          levels = level_scales * hom.md, col = hom.col)
  abline(v = 0, lty = 2)
  abline(h = 0, lty = 2)
  HAPSEG:::PlotSnpAscn(d, snp.clust.p, HAPSEG:::HTxInv(e.mu, sigma.nu, sigma.eta), 
              "Allelic copy-ratio", xcrds = snp.annot[["pos"]]/1e+06, 
              chr = snp.annot[["chr"]], seg.wd = 1, loci.annot = merged.loci, 
              merge.prob = merge.prob, seg.log.ev = seg.log.ev, min = min, 
              max = max)
  PlotSnpAscn(HTx(d, sigma.nu, sigma.eta), snp.clust.p, e.mu, 
              y.lab = "Variance stabilized ", xcrds = snp.annot[["pos"]]/1e+06, 
              chr = snp.annot[["chr"]], seg.wd = 1, loci.annot = merged.loci, 
              merge.prob = merge.prob, seg.log.ev = seg.log.ev, min = t.min, 
              max = t.max)
}
  
HAPSEG:::PlotSnpAscn
function (d, snp.clust.p, e.mu, y.lab, xcrds = NA, chr = NA, 
          seg.wd = 4, add = FALSE, vertical = FALSE, plot.markers = TRUE, 
          col.markers = TRUE, plot.segs = TRUE, het.seg.col.1 = NA, 
          het.seg.col.2 = NA, loci.annot = NA, merge.prob = NA, seg.log.ev = NA, 
          plot.tcn.seg = TRUE, min = NA, max = NA, marker.gt.pal = NA, 
          pch = ".", cex = 1.2, xaxt = "s") 
{
  max.c <- apply(snp.clust.p, 1, which.max)
  out.color <- "white"
  het.col <- "green"
  hom.col <- "black"
  if (any(is.na(marker.gt.pal))) {
    pal <- GetHapsegMarkerPal()
  }
  else {
    pal <- marker.gt.pal
  }
  cols <- pal(1000)
  state.mat <- matrix(c(0, 1, 2, 3), ncol = 4, nrow = nrow(snp.clust.p), 
                      byrow = TRUE)
  col.ix.1 <- round(rowSums(state.mat * snp.clust.p[, c(3, 
                                                        2, 4, 1)])/3 * 999, 0) + 1
  col.ix.2 <- round(rowSums(state.mat * snp.clust.p[, c(1, 
                                                        4, 2, 3)])/3 * 999, 0) + 1
  if (col.markers) {
    pcols.1 <- cols[col.ix.1]
    pcols.2 <- cols[col.ix.2]
  }
  else {
    pcols.1 <- rep(1, length(col.ix.1))
    pcols.2 <- rep(1, length(col.ix.2))
  }
  if (is.na(het.seg.col.1)) {
    het.seg.col.1 <- het.col
  }
  if (is.na(het.seg.col.2)) {
    het.seg.col.2 <- het.col
  }
  hom.seg.col <- hom.col
  if (!add) {
    x.lab <- ifelse(xaxt == "n", "", paste("Chr", chr, "position (Mb)"))
    plot(0, type = "n", xlab = x.lab, ylab = y.lab, ylim = c(min, 
                                                             max), xlim = range(xcrds), main = "", xaxt = xaxt)
    if (any(is.na(xcrds))) {
      xcrds <- c(1:ncol(d))
      plot(0, type = "n", xlab = "Marker order", ylab = y.lab, 
           ylim = c(min, max), xlim = range(xcrds), main = "")
    }
  }
  if (plot.markers) {
    hom.ix <- which((max.c == 1) | (max.c == 3))
    points(xcrds[hom.ix], d[1, hom.ix], col = pcols.1[hom.ix], 
           pch = pch, cex = cex)
    points(xcrds[hom.ix], d[2, hom.ix], col = pcols.2[hom.ix], 
           pch = pch, cex = cex)
    if (sum(max.c == 2) > 0) {
      het.ix <- which((max.c == 2) | (max.c == 4))
      points(xcrds[het.ix], d[1, het.ix], col = pcols.1[het.ix], 
             pch = pch)
      points(xcrds[het.ix], d[2, het.ix], col = pcols.2[het.ix], 
             pch = pch)
    }
  }
  wd <- 0.04
  if (plot.segs) {
    if (!vertical) {
      if (add) {
        rect(xleft = min(xcrds), ybottom = e.mu[1] - 
               wd/2, xright = max(xcrds), ytop = e.mu[1] + 
               wd/2, border = NA, col = het.seg.col.1)
        rect(xleft = min(xcrds), ybottom = e.mu[2] - 
               wd/2, xright = max(xcrds), ytop = e.mu[2] + 
               wd/2, border = NA, col = het.seg.col.2)
      }
      else {
        lines(x = range(xcrds), y = c(e.mu[1], e.mu[1]), 
              col = het.seg.col.1, lwd = seg.wd)
        lines(x = range(xcrds), y = c(e.mu[2], e.mu[2]), 
              col = het.seg.col.2, lwd = seg.wd)
      }
      if (plot.tcn.seg) {
        lines(x = range(xcrds), y = c(e.mu[3], e.mu[3]), 
              col = hom.seg.col, lwd = seg.wd)
        lines(x = range(xcrds), y = c(0, 0), col = hom.seg.col, 
              lwd = seg.wd)
      }
    }
  }
  if (!is.na(loci.annot)) {
    n <- length(xcrds)
    msg <- paste("N = ", n, ", log_ev/N = ", round(seg.log.ev/n, 
                                                   3), ", MPr = ", round(merge.prob[1], 5), ",\nlog10 MPr = ", 
                 round(merge.prob[2], 2), sep = "")
    mtext(msg, line = 0, side = 3, adj = 0, cex = par("cex") * 
            par("cex.axis"))
    ix <- (loci.annot[, "chr"] == chr)
    if (sum(ix) == 0) {
      return()
    }
    loci.annot <- loci.annot[ix, ]
    gcrd <- loci.annot[, "pos"]/1e+06
    if (!vertical) {
      segments(x0 = gcrd, y0 = -1, x1 = gcrd, y1 = 4, lty = 2)
      text(gcrd, y = rep(4.5, length(gcrd)), labels = paste("Pr=", 
                                                            round(loci.annot[, "log10_prob"], 5), sep = ""), 
           cex = 0.5, srt = 90)
    }
    else {
      lines(y = loci.annot[, "pos"]/1e+06, x = c(-1, 4), 
            lty = 2)
      text(y = gcrd, x = rep(5, length(gcrd)), labels = paste("Pr=", 
                                                              round(loci.annot[, "prob"], 5), sep = ""), cex = 0.5)
    }
  }
}
<environment: namespace:HAPSEG>
  
  >                                     
  
myInitAndMergeSmall = function (array.name, genome.build, use.pop, use.normal, normal, 
          impute.gt, platform, seg.fn, snp.fn, drop.x, drop.y, calls.fn, 
          mn.sample, min.seg.size, merge.close, out.p, merge.small, 
          force.diploid, calibrate.data, clusters.fn, prev.theta.fn, 
          snp.file.parser = DefaultSnpFileParser, clusters.file.parser = DefaultClustersFileParser, 
          verbose = FALSE) 
{
  SetPlatformSpecificFuncs(platform)
  res <- myExtractSnpDat(array.name, genome.build, use.pop, use.normal, 
                       normal, impute.gt, platform, seg.fn, snp.fn, drop.x, 
                       drop.y, calls.fn, mn.sample, calibrate.data = calibrate.data, 
                       clusters.fn = clusters.fn, verbose = verbose, snp.file.parser = snp.file.parser, 
                       clusters.file.parser = clusters.file.parser)
  if (verbose) {
    print(paste("Setting platform specific functions for platform", 
                platform))
    PrintHapSegStartMessage()
  }
  h.seg.dat <- res[["as.res"]][["h.seg.dat"]]
  h.d <- h.seg.dat[["h.d"]]
  h.snp.gt.p <- h.seg.dat[["h.snp.gt.p"]]
  h.snp.annot <- h.seg.dat[["h.snp.annot"]]
  seg.info <- h.seg.dat[["segmentation.information"]]
  seg.dat <- res[["seg.dat"]]
  seg.dat[["matched.name"]] <- res[["as.res"]][["matched.name"]]
  seg.dat[["found.matched.normal"]] <- res[["as.res"]][["found.matched.normal"]]
  if (merge.small) {
    merge.res <- JoinSmallSegs(h.d, min.seg.size, h.snp.gt.p, 
                               h.snp.annot, verbose = verbose)
    h.d <- merge.res[["h.d"]]
    h.snp.gt.p <- merge.res[["h.snp.gt.p"]]
    h.snp.annot <- merge.res[["h.snp.annot"]]
    if (verbose) {
      print(paste(length(h.d), " new segments.", sep = ""))
    }
  }
  theta <- InitTheta(array.name, prev.theta.fn, verbose)
  use.eps <- ifelse((merge.close || impute.gt), 0.005, 1e-04)
  em.fit <- HscrSegFit(h.d, h.snp.gt.p, h.snp.annot, theta, 
                       seg.info, eps = use.eps, out.p = out.p, force.diploid = force.diploid, 
                       verbose = verbose)
  theta <- em.fit[["theta"]]
  return(list(h.seg.dat = h.seg.dat, h.d = h.d, h.snp.gt.p = h.snp.gt.p, 
              h.snp.annot = h.snp.annot, seg.info = seg.info, seg.dat = seg.dat, 
              theta = theta, use.eps = use.eps, em.fit = em.fit, theta = theta))
}
  
  
  
myExtractSnpDat = function (sample.name, genome.build, use.pop, use.normal, normal, 
          impute.gt, platform, glad.fn, snp.fn, drop.x, drop.y, calls.fn, 
          mn.sample, calibrate.data = FALSE, clusters.fn = NULL, verbose = FALSE, 
          snp.file.parser = DefaultSnpFileParser, clusters.file.parser = DefaultClusterFileParser) 
{
  data(list = GetPlatformDataName(platform), package = "HAPSEG")
  if (!genome.build %in% names(platform.annots)) {
    stop("Unsupported genome build, ", genome.build, ", for this platform: ", 
         platform)
  }
  platform.vals <- mget(c("snp.freqs", "dbSNP.annot", "snp.annot", 
                          "post.birdseed.calibration"), env = platform.annots[[genome.build]])
  snp.freqs <- platform.vals[["snp.freqs"]]
  dbSNP.annot <- platform.vals[["dbSNP.annot"]]
  snp.annot <- platform.vals[["snp.annot"]]
  post.birdseed.calibration <- platform.vals[["post.birdseed.calibration"]]
  allele.data <- snp.file.parser(snp.fn, sample.name, verbose = verbose)
  calibrate.data = infer_calibration_status(calibrate.data, 
                                            snp.file.parser)
  if (calibrate.data) {
    if (is.null(clusters.fn)) {
      stop("Calibrating data, but there is no clusters file!")
    }
    if (verbose) {
      print("Calibrating data")
    }
    as.d <- CalibrateAsDat(allele.data, clusters.fn, clusters.file.parser = clusters.file.parser, 
                           verbose = verbose)
  }
  else {
    if (!verbose) {
      print("Not calibrating data")
    }
    as.d <- allele.data
  }
  as.d <- PostCalibrateAsDat(as.d, post.birdseed.calibration[["snp.tx"]], 
                             verbose = verbose)
  nas <- apply(is.na(as.d), 1, sum)
  as.d <- as.d[nas == 0, ]
  if (is.null(glad.fn)) {
    glad.mat = run_cbs(as.d, genome.build, sample.name, verbose = verbose)
  }
  else {
    if (verbose) {
      print(paste("loading glad mat", glad.fn))
    }
    glad.mat <- read.delim(glad.fn, row.names = NULL, as.is = TRUE)
  }
  seg.dat <- ReadGladMat(glad.mat, sample.name, glad.log = TRUE, 
                         drop.x = drop.x, drop.y = drop.y, verbose = verbose)
  as.res <- GetAlleleSegData(as.d, snp.annot, seg.dat[["seg.info"]], 
                             snp.freqs, use.pop, use.normal, normal, NA, dbSNP.annot, 
                             impute.gt, calls.fn, mn.sample, verbose = verbose)
  return(list(seg.dat = seg.dat, as.res = as.res))
}
  


HAPSEG:::MergeCloseAndFit
function (in.list, out.p, seg.merge.thresh, impute.gt, platform, 
          force.diploid = FALSE, verbose = FALSE) 
{
  SetPlatformSpecificFuncs(platform)
  h.seg.dat <- in.list[["h.seg.dat"]]
  h.d <- in.list[["h.d"]]
  h.snp.gt.p <- in.list[["h.snp.gt.p"]]
  h.snp.annot <- in.list[["h.snp.annot"]]
  theta <- in.list[["theta"]]
  seg.dat <- in.list[["seg.dat"]]
  seg.info <- in.list[["seg.info"]]
  use.eps <- in.list[["use.eps"]]
  mrg.res <- JoinCloseSegs(h.d, h.snp.gt.p, h.snp.annot, theta, 
                           force.diploid, out.p, merge.thresh = seg.merge.thresh, 
                           verbose = verbose)
  h.d <- mrg.res[["h.d"]]
  h.snp.gt.p <- mrg.res[["h.snp.gt.p"]]
  h.snp.annot <- mrg.res[["h.snp.annot"]]
  if (verbose) {
    print("Merged Loci:")
  }
  seg.dat[["merged.loci"]] <- mrg.res[["merged.loci"]]
  seg.dat[["final.merge.prob"]] <- mrg.res[["final.merge.prob"]]
  if (verbose) {
    print(paste(length(h.d), " merged segs.", sep = ""))
  }
  use.eps <- ifelse(impute.gt, 0.01, 1e-04)
  em.fit <- HscrSegFit(h.d, h.snp.gt.p, h.snp.annot, theta, 
                       seg.info, eps = use.eps, out.p = out.p, force.diploid = force.diploid, 
                       verbose = verbose)
  theta <- em.fit[["theta"]]
  seg.dat[["seg.expected.phase"]] <- em.fit[["seg.expected.phase"]]
  return(list(h.d = h.d, h.snp.gt.p = h.snp.gt.p, h.seg.dat = h.seg.dat, 
              h.snp.annot = h.snp.annot, seg.dat = seg.dat, use.eps = use.eps, 
              em.fit = em.fit, theta = theta, seg.info = seg.info))
}
<environment: namespace:HAPSEG>
  > 
 
  
  
HAPSEG:::InitAndMergeSmall
function (array.name, genome.build, use.pop, use.normal, normal, 
          impute.gt, platform, seg.fn, snp.fn, drop.x, drop.y, calls.fn, 
          mn.sample, min.seg.size, merge.close, out.p, merge.small, 
          force.diploid, calibrate.data, clusters.fn, prev.theta.fn, 
          snp.file.parser = DefaultSnpFileParser, clusters.file.parser = DefaultClustersFileParser, 
          verbose = FALSE) 
{
  SetPlatformSpecificFuncs(platform)
  res <- ExtractSnpDat(array.name, genome.build, use.pop, use.normal, 
                       normal, impute.gt, platform, seg.fn, snp.fn, drop.x, 
                       drop.y, calls.fn, mn.sample, calibrate.data = calibrate.data, 
                       clusters.fn = clusters.fn, verbose = verbose, snp.file.parser = snp.file.parser, 
                       clusters.file.parser = clusters.file.parser)
  if (verbose) {
    print(paste("Setting platform specific functions for platform", 
                platform))
    PrintHapSegStartMessage()
  }
  h.seg.dat <- res[["as.res"]][["h.seg.dat"]]
  h.d <- h.seg.dat[["h.d"]]
  h.snp.gt.p <- h.seg.dat[["h.snp.gt.p"]]
  h.snp.annot <- h.seg.dat[["h.snp.annot"]]
  seg.info <- h.seg.dat[["segmentation.information"]]
  seg.dat <- res[["seg.dat"]]
  seg.dat[["matched.name"]] <- res[["as.res"]][["matched.name"]]
  seg.dat[["found.matched.normal"]] <- res[["as.res"]][["found.matched.normal"]]
  if (merge.small) {
    merge.res <- JoinSmallSegs(h.d, min.seg.size, h.snp.gt.p, 
                               h.snp.annot, verbose = verbose)
    h.d <- merge.res[["h.d"]]
    h.snp.gt.p <- merge.res[["h.snp.gt.p"]]
    h.snp.annot <- merge.res[["h.snp.annot"]]
    if (verbose) {
      print(paste(length(h.d), " new segments.", sep = ""))
    }
  }
  theta <- InitTheta(array.name, prev.theta.fn, verbose)
  use.eps <- ifelse((merge.close || impute.gt), 0.005, 1e-04)
  em.fit <- HscrSegFit(h.d, h.snp.gt.p, h.snp.annot, theta, 
                       seg.info, eps = use.eps, out.p = out.p, force.diploid = force.diploid, 
                       verbose = verbose)
  theta <- em.fit[["theta"]]
  return(list(h.seg.dat = h.seg.dat, h.d = h.d, h.snp.gt.p = h.snp.gt.p, 
              h.snp.annot = h.snp.annot, seg.info = seg.info, seg.dat = seg.dat, 
              theta = theta, use.eps = use.eps, em.fit = em.fit, theta = theta))
}


HAPSEG:::ExtractSnpDat
function (sample.name, genome.build, use.pop, use.normal, normal, 
          impute.gt, platform, glad.fn, snp.fn, drop.x, drop.y, calls.fn, 
          mn.sample, calibrate.data = FALSE, clusters.fn = NULL, verbose = FALSE, 
          snp.file.parser = DefaultSnpFileParser, clusters.file.parser = DefaultClusterFileParser) 
{
  data(list = GetPlatformDataName(platform), package = "HAPSEG")
  if (!genome.build %in% names(platform.annots)) {
    stop("Unsupported genome build, ", genome.build, ", for this platform: ", 
         platform)
  }
  platform.vals <- mget(c("snp.freqs", "dbSNP.annot", "snp.annot", 
                          "post.birdseed.calibration"), env = platform.annots[[genome.build]])
  snp.freqs <- platform.vals[["snp.freqs"]]
  dbSNP.annot <- platform.vals[["dbSNP.annot"]]
  snp.annot <- platform.vals[["snp.annot"]]
  post.birdseed.calibration <- platform.vals[["post.birdseed.calibration"]]
  allele.data <- snp.file.parser(snp.fn, sample.name, verbose = verbose)
  calibrate.data = infer_calibration_status(calibrate.data, 
                                            snp.file.parser)
  if (calibrate.data) {
    if (is.null(clusters.fn)) {
      stop("Calibrating data, but there is no clusters file!")
    }
    if (verbose) {
      print("Calibrating data")
    }
    as.d <- CalibrateAsDat(allele.data, clusters.fn, clusters.file.parser = clusters.file.parser, 
                           verbose = verbose)
  }
  else {
    if (!verbose) {
      print("Not calibrating data")
    }
    as.d <- allele.data
  }
  as.d <- PostCalibrateAsDat(as.d, post.birdseed.calibration[["snp.tx"]], 
                             verbose = verbose)
  nas <- apply(is.na(as.d), 1, sum)
  as.d <- as.d[nas == 0, ]
  if (is.null(glad.fn)) {
    glad.mat = run_cbs(as.d, genome.build, sample.name, verbose = verbose)
  }
  else {
    if (verbose) {
      print(paste("loading glad mat", glad.fn))
    }
    glad.mat <- read.delim(glad.fn, row.names = NULL, as.is = TRUE)
  }
  seg.dat <- ReadGladMat(glad.mat, sample.name, glad.log = TRUE, 
                         drop.x = drop.x, drop.y = drop.y, verbose = verbose)
  as.res <- GetAlleleSegData(as.d, snp.annot, seg.dat[["seg.info"]], 
                             snp.freqs, use.pop, use.normal, normal, NA, dbSNP.annot, 
                             impute.gt, calls.fn, mn.sample, verbose = verbose)
  return(list(seg.dat = seg.dat, as.res = as.res))
}
<environment: namespace:HAPSEG>
  
  >
  
  
  
  HAPSEG:::GetAlleleSegData
function (as.d, snp.annot, glad.mat, snp.freqs, use.pop, use.normal, 
          normal, bad.snps, dbSNP.annot, impute.gt, mn.calls.fn, mn.sample, 
          verbose = FALSE) 
{
  sna <- apply(is.na(as.d), 1, sum)
  as.d <- as.d[sna == 0, ]
  as.d <- as.d[intersect(rownames(snp.annot), rownames(as.d)), 
               ]
  as.d <- cbind(snp.annot[rownames(as.d), ], as.d)
  if (!is.na(bad.snps)) {
    bad <- which(rownames(as.d) %in% bad.snps)
    as.d <- as.d[-bad, ]
    if (verbose) {
      print(paste("Removing ", sum(bad), " 'bad' SNPs", 
                  sep = ""))
    }
  }
  found.matched.normal <- FALSE
  mn.calls <- NULL
  if (!normal && use.normal && (!is.null(mn.sample)) && (!is.na(mn.sample)) && 
        (mn.sample != "NA")) {
    if (is.null(mn.calls.fn)) {
      stop("Matched normal sample specified, but no calls file provided")
    }
    if (verbose) {
      print("reading column data")
    }
    mn.calls <- ReadCol(mn.calls.fn, mn.sample, save.rownames = TRUE)
  }
  if (verbose) {
    print(summary(mn.calls))
  }
  if (use.normal && (!is.null(mn.calls))) {
    found.matched.normal <- TRUE
    mn.calls <- mn.calls[rownames(as.d), ]
    snp.gt.p <- matrix(0, nrow = 4, ncol = length(mn.calls), 
                       byrow = TRUE)
    snp.gt.p[, mn.calls == 1] <- c(0.015, 0.985, 0.015, 0.985)/2
    snp.gt.p[, mn.calls == 0] <- c(0.005, 0.005, 0.985, 0.005)
    snp.gt.p[, mn.calls == 2] <- c(0.985, 0.005, 0.005, 0.005)
    snp.gt.p[, mn.calls == -1] <- c(0.39, 0.11, 0.39, 0.11)
    snp.gt.p <- t(snp.gt.p)
    if (verbose) {
      n.het <- sum(mn.calls == 1)
      print(paste("Normal ", round(n.het/length(mn.calls) * 
                                     100, 2), "% Het", sep = ""))
    }
  }
  else {
    pop.cols <- c(paste(use.pop, "_A", sep = ""), paste(use.pop, 
                                                        "_B", sep = ""))
    populations = gsub("_A", "", colnames(snp.freqs)[grep("_A$", 
                                                          colnames(snp.freqs))])
    if (!use.pop %in% populations) {
      stop("Supplied population ", use.pop, " is not supported. Supported populations are ", 
           paste(populations, collapse = ", "))
    }
    if (is.na(use.pop) || use.pop == "NA") {
      if (verbose) {
        print("Using Null pop.")
      }
      snp.gt.p <- matrix(c(0.39, 0.11, 0.39, 0.11), nrow = nrow(as.d), 
                         ncol = 4, byrow = TRUE)
    }
    else {
      if (verbose) {
        print(paste("Using population allele-frequency data for: ", 
                    use.pop, sep = ""))
      }
      snp.freqs <- snp.freqs[intersect(rownames(as.d), 
                                       rownames(snp.freqs)), ]
      freqs <- snp.freqs[, pop.cols] + 0.01
      freqs <- freqs/rowSums(freqs)
      aa.probs <- freqs[, 1]^2
      bb.probs <- freqs[, 2]^2
      het.probs <- 1 - (aa.probs + bb.probs)
      snp.gt.p <- rbind(bb.probs, het.probs/2, aa.probs, 
                        het.probs/2)
      colnames(snp.gt.p) <- rownames(as.d)
      rownames(snp.gt.p)[c(1, 3)] <- c("BB_probs", "AA_probs")
      na.ix <- apply(is.na(snp.gt.p), 2, sum) > 0
      snp.gt.p[, na.ix] <- apply(snp.gt.p, 1, mean, na.rm = TRUE)
      snp.gt.p <- t(snp.gt.p)
    }
  }
  as.d[, 1] <- as.character(as.d[, 1])
  h.seg.dat <- GetAsSegs(glad.mat, as.d, snp.gt.p, dbSNP.annot, 
                         impute.gt, verbose = verbose)
  return(list(h.seg.dat = h.seg.dat, mn.sample = mn.sample, 
              found.matched.normal = found.matched.normal))
}
<environment: namespace:HAPSEG>
  
  >  
  
  HAPSEG:::ReadGladMat
function (glad.mat, sample.name = NULL, drop.x = TRUE, drop.y = TRUE, 
          glad.log = FALSE, verbose = FALSE) 
{
  glad.mat[, 2] <- gsub("X", "23", glad.mat[, 2])
  glad.mat[, 2] <- gsub("Y", "24", glad.mat[, 2])
  glad.mat[, 2] <- as.numeric(glad.mat[, 2])
  if (verbose) {
    print("summary pre-processing")
  }
  glad.mat[, 2] <- sub("chr", "", glad.mat[, 2])
  colnames(glad.mat) <- c("Sample", "Chromosome", "Start.bp", 
                          "End.bp", "NUM.SNPs", "Seg.CN")
  if (!is.null(sample.name)) {
    if (verbose) {
      print(paste("Extracting out sample", sample.name, 
                  "from the seg file"))
    }
    idx <- glad.mat[, 1] == sample.name
    if (sum(idx) == 0) {
      stop("Unable to find any data in the sample column matching the sample name, ", 
           "even though it was specified.  Check your segmentation file")
    }
    glad.mat <- glad.mat[idx, ]
    if (nrow(glad.mat) == 0) {
      stop()
    }
    if (verbose) {
      print("read GLAD mat: ")
      print(dim(glad.mat))
    }
  }
  dropChroms <- character()
  if (drop.x) {
    dropChroms <- c(dropChroms, c("23", "X"))
  }
  if (drop.y) {
    dropChroms <- c(dropChroms, c("24", "Y"))
  }
  if (length(dropChroms) > 0) {
    glad.mat <- DropChromosomes(glad.mat, dropChroms)
  }
  seg.lens <- (glad.mat[, "End.bp"] - glad.mat[, "Start.bp"])/1e+07
  glad.mat <- cbind(glad.mat, seg.lens)
  colnames(glad.mat)[ncol(glad.mat)] <- "length"
  data <- list(glad = glad.mat)
  seg.info <- list()
  NSEG <- nrow(glad.mat)
  cnames <- c("Chromosome", "Start.bp", "End.bp", "n_probes", 
              "length", "copy_num")
  seg.info <- matrix(NA, nrow = NSEG, ncol = length(cnames))
  colnames(seg.info) <- cnames
  cols <- c("Chromosome", "Start.bp", "End.bp")
  seg.info[, cols[1]] <- glad.mat[, cols[1]]
  seg.info[, cols[2]] <- glad.mat[, cols[2]]
  seg.info[, cols[3]] <- glad.mat[, cols[3]]
  if (glad.log == TRUE) {
    seg.info[, "copy_num"] <- 2^data[["glad"]][, "Seg.CN"]
  }
  else {
    seg.info[, "copy_num"] <- data[["glad"]][, "Seg.CN"]
  }
  seg.info[, "length"] <- data[["glad"]][, "length"]
  seg.info[, "n_probes"] <- data[["glad"]][, "NUM.SNPs"]
  seg.info <- matrix(as.numeric(seg.info), nrow = nrow(seg.info), 
                     ncol = ncol(seg.info))
  colnames(seg.info) <- cnames
  if (verbose) {
    print(seg.info[1:10, 1])
  }
  return(list(seg.info = seg.info))
}
<environment: namespace:HAPSEG>
  
  >
  
  
  
  
  HAPSEG:::run_cbs
function (as_d, genome_build, sample_name, verbose = FALSE) 
{
  set.seed(0)
  if (verbose) {
    print("No segmentation data provided, creating ....")
  }
  require("DNAcopy") || stop("DNAcopy required to perform segmentation")
  data(markers)
  markers = markers[[genome_build]]
  shared_markers = intersect(rownames(as_d), markers[, 1])
  markers = markers[which(markers[, 1] %in% shared_markers), 
                    ]
  chroms = markers[, 2]
  chroms = gsub("X", "23", chroms)
  chroms = gsub("Y", "24", chroms)
  chroms = as.numeric(chroms)
  pos = as.numeric(markers[, 3])
  as_d = as_d[markers[, 1], ]
  dat = as_d[, 1] + as_d[, 2]
  names(dat) = rownames(as_d)
  dat = log2(dat) - 1
  cna = CNA(dat, chroms, pos, data.type = "logratio", sampleid = sample_name)
  smoothed_cna = smooth.CNA(cna)
  segment_cna = segment(smoothed_cna, min.width = 2, nperm = 10000, 
                        alpha = 0.01, undo.splits = "sdundo", undo.prune = 0.05, 
                        undo.SD = 1, verbose = 2)$output
  colnames(segment_cna) = c("Sample", "Chromosome", "Start", 
                            "End", "Num_Probes", "Segment_Mean")
  return(segment_cna)
}
<environment: namespace:HAPSEG>
  
  >
  
  
  
  getAnywhere('CNA')
A single object matching ‘CNA’ was found
It was found in the following places
namespace:DNAcopy
with value

function (genomdat, chrom, maploc, data.type = c("logratio", 
                                                 "binary"), sampleid = NULL, presorted = FALSE) 
{
  if (is.data.frame(genomdat)) 
    genomdat <- as.matrix(genomdat)
  if (!is.numeric(genomdat)) 
    stop("genomdat must be numeric")
  if (!is.numeric(maploc)) 
    stop("maploc must be numeric")
  data.type <- match.arg(data.type)
  ina <- (!is.na(chrom) & is.finite(maploc))
  if (sum(!ina) > 0) 
    warning("markers with missing chrom and/or maploc removed\n")
  if (!presorted) {
    sortindex <- which(ina)[order(chrom[ina], maploc[ina])]
  }
  else {
    sortindex <- which(ina)
  }
  if (is.factor(chrom)) 
    chrom <- as.character(chrom)
  if (is.vector(genomdat)) 
    genomdat <- as.matrix(genomdat)
  if (!missing(sampleid)) {
    if (length(sampleid) != ncol(genomdat)) {
      warning("length(sampleid) and ncol(genomdat) differ, names ignored\n")
      sampleid <- paste("Sample", 1:ncol(genomdat))
    }
  }
  else {
    sampleid <- paste("Sample", 1:ncol(genomdat))
  }
  colnames(genomdat) <- sampleid
  zzz <- data.frame(chrom = I(chrom), maploc = maploc, genomdat)
  zzz <- zzz[sortindex, ]
  if (length(ii <- which(diff(maploc) == 0)) > 0) {
    if (any(chrom[ii] == chrom[ii + 1])) 
      warning("array has repeated maploc positions\n")
  }
  attr(zzz, "data.type") <- data.type
  class(zzz) <- c("CNA", "data.frame")
  zzz
}
<environment: namespace:DNAcopy>
  
  >
  
  
HAPSEG:::PostCalibrateAsDat
function (dat, snp.tx, verbose = FALSE) 
{
  InvAtten <- function(r, at) {
    r/(1 + at - (at * r))
  }
  dat <- dat[intersect(rownames(dat), rownames(snp.tx)), ]
  snp.tx <- snp.tx[rownames(dat), ]
  n.snp <- nrow(snp.tx)
  tx.snp.dat <- array(NA, dim = c(n.snp, 2))
  tx.snp.dat[, 1] <- dat[, 1] + snp.tx[, 2] * dat[, 2]
  tx.snp.dat[, 2] <- dat[, 1] * snp.tx[, 1] + dat[, 2] * (1 + 
                                                            apply(snp.tx[, c(1:2)], 1, prod))
  tx.snp.dat <- (tx.snp.dat + snp.tx[, c(4, 3)]) * snp.tx[, 
                                                          c(6, 5)]
  snp.tx[is.na(snp.tx)] <- 0
  return(tx.snp.dat)
}
<environment: namespace:HAPSEG>
  
  >
  
  
HAPSEG:::DoPlots
function (seg.dat, results.dir, platform, verbose = FALSE) 
{
  if (!file.exists(results.dir)) {
    dir.create(results.dir, recursive = TRUE)
  }
  HAPSEG:::SetPlatformSpecificFuncs(platform)
  h.d <- seg.dat[["as.seg.dat"]]
  n.segs <- length(h.d)
  seg.plot.ix <- seq_along(h.d)
  theta <- seg.dat[["em.res"]][["theta"]]
  theta[["p.snp.cond.out"]] <- 1/25
  segment.chroms <- seg.dat[["allele.segs"]][, "Chromosome"]
  unique.seg.chroms <- unique(segment.chroms)
  for (i in seq_along(unique.seg.chroms)) {
    dir.create(HAPSEG:::GetChromPlotDir(unique.seg.chroms[i], results.dir))
  }
  merge.prob <- rbind(seg.dat[["final.merge.prob"]], c(NA, 
                                                       NA))
  em.fit <- seg.dat[["em.res"]]
  h.snp.annot <- em.fit[["h.snp.annot"]]
  for (i in seq_len(n.segs)) {
    segfit.fn <- file.path(HAPSEG:::GetChromPlotDir(segment.chroms[i], 
                                           results.dir), paste("HAPSEG_SEG_", i, ".jpg", sep = ""))
    jpeg(segfit.fn, 7, 5, units = "in", type = "cairo", res = 200, 
         quality = 100)
    par(cex = 0.5)
    par(cex.axis = 0.75)
    par(cex.lab = 0.75)
    par(las = 1)
    par(mfcol = c(2, 3))
    if (ncol(h.d[[i]]) > 5) {
      myPlotSegFit(d = h.d[[i]] - theta[["bg"]], snp.gt.p = em.fit[["h.het.prob"]][[i]], 
                snp.clust.p = em.fit[["snp.clust.p"]][[i]], snp.annot = h.snp.annot[[i]], 
                          e.mu = em.fit[["e.mu"]][i, ], theta, merged.loci = seg.dat[["merged.loci"]], 
                          merge.prob = merge.prob[i, ], seg.log.ev = em.fit[["seg.log.ev"]][i], verbose = verbose)
    }
    dev.off()
    if (verbose) {
      cat(".")
    }
  }
} 



myPlotSegFit = function (d, snp.gt.p, snp.clust.p, snp.annot, e.mu, theta, merged.loci, 
          merge.prob, seg.log.ev, marker.gt.pal = NA, verbose = FALSE) 
{
  data.name <- "Calibrated Intensities"
  par(bty = "n")
  if (verbose) {
    print(summary(snp.gt.p))
  }
  sigma.nu <- theta[["sigma.nu"]]
  sigma.eta <- theta[["sigma.eta"]]
  nu <- theta[["nu"]]
  het.cov <- theta[["het.cov"]]
  het.cov <- min(het.cov, (het.cov * e.mu[1] * e.mu[2]))
  mu0 <- 0
  sigma.h <- HAPSEG:::GetSigmaH(sigma.nu, sigma.eta)
  ct <- colSums(snp.gt.p)
  ct <- ct/sum(ct)
  ngt.t <- c((ct[1] + ct[3]), (ct[2] + ct[4]))
  cex <- 1.2
  min <- -0.5
  max <- 5.25
  t.min <- -0.5
  t.max <- 2.75
  het.col <- "green"
  hom.col <- "black"
  df <- d[((d >= min) & (d < max))]
  if (length(df)<1){
    min <- min(d)  #I have changed to this because of high CN on chr 17
    max <- max(d)    
    df <- d[((d >= min) & (d < max))]
  }
  hist(df, breaks = 100, freq = FALSE, main = "", xlab = "Allelic copy-ratio", 
       xlim = c(min, max))
  xgl <- 1001
  xg <- array(NA, dim = c(4, xgl))
  x <- seq(min(df), max(df), length.out = xgl)
  x.tx <- HAPSEG:::HTx(x, sigma.nu, sigma.eta)
  e.mu <- HAPSEG:::HTx(e.mu, sigma.nu, sigma.eta)
  mu0.tx <- HAPSEG:::HTx(mu0, sigma.nu, sigma.eta)
  xg[1, ] <- ngt.t[1] * HAPSEG:::DFunc(x.tx, mu0.tx, sigma.h, nu) * 
    HAPSEG:::HTxVar(x, sigma.nu, sigma.eta)
  xg[2, ] <- ngt.t[2] * HAPSEG:::DFunc(x.tx, e.mu[1], sigma.h, nu) * 
    HAPSEG:::HTxVar(x, sigma.nu, sigma.eta)
  xg[3, ] <- ngt.t[2] * HAPSEG:::DFunc(x.tx, e.mu[2], sigma.h, nu) * 
    HAPSEG:::HTxVar(x, sigma.nu, sigma.eta)
  xg[4, ] <- ngt.t[1] * HAPSEG:::DFunc(x.tx, e.mu[3], sigma.h, nu) * 
    HAPSEG:::HTxVar(x, sigma.nu, sigma.eta)
  csum <- apply(xg, 2, sum)
  crds <- par("usr")
  scale <- 1/2
  for (i in 1:4) {
    lines(x, scale * xg[i, ], lwd = 2, col = hom.col)
  }
  lines(x, scale * csum, lwd = 2, col = "coral")
  xg[1, ] <- ngt.t[1] * HAPSEG:::DFunc(x.tx, mu0.tx, sigma.h, nu)
  xg[2, ] <- ngt.t[2] * HAPSEG:::DFunc(x.tx, e.mu[1], sigma.h, nu)
  xg[3, ] <- ngt.t[2] * HAPSEG:::DFunc(x.tx, e.mu[2], sigma.h, nu)
  xg[4, ] <- ngt.t[1] * HAPSEG:::DFunc(x.tx, e.mu[3], sigma.h, nu)
  csum <- apply(xg, 2, sum)
  d.tx <- HAPSEG:::HTx(d, sigma.nu, sigma.eta)
  d.tx.f <- d.tx[(d.tx > t.min) & (d.tx < t.max)]
  hist(d.tx.f, breaks = 100, freq = FALSE, main = "", xlab = "Allelic copy-ratio")
  mtext("Variance-stabilizing transformation", side = 3, line = 0, 
        adj = 0, cex = par("cex") * par("cex.axis"))
  for (i in 1:4) {
    lines(x.tx, scale * xg[i, ], lwd = 2, col = hom.col)
  }
  lines(x.tx, scale * csum, lwd = 2, col = "coral")
  gl <- 250
  x2d <- seq(t.min, t.max, length.out = gl)
  x1 <- rep(x2d, gl)
  x2 <- as.vector(sapply(x2d, rep, gl))
  x <- cbind(x1, x2)
  hom.cov <- 0
  sigma <- matrix(c(sigma.h^2, het.cov, het.cov, sigma.h^2), 
                  nrow = 2, ncol = 2)
  invert.sigma <- solve(sigma)
  hom.sigma <- matrix(c(sigma.h^2, 0, 0, sigma.h^2), nrow = 2, 
                      ncol = 2)
  invert.hom.sigma <- solve(hom.sigma)
  het.1 <- array(scale * ngt.t[2] * DmvFunc(x, c(e.mu[1], e.mu[2]), 
                                            sigma, invert.sigma, nu), dim = c(gl, gl))
  het.2 <- array(scale * ngt.t[2] * DmvFunc(x, c(e.mu[2], e.mu[1]), 
                                            sigma, invert.sigma, nu), dim = c(gl, gl))
  hom.1 <- array(scale * ngt.t[1] * DmvFunc(x, c(mu0.tx, e.mu[3]), 
                                            hom.sigma, invert.hom.sigma, nu), dim = c(gl, gl))
  hom.2 <- array(scale * ngt.t[1] * DmvFunc(x, c(e.mu[3], mu0.tx), 
                                            hom.sigma, invert.hom.sigma, nu), dim = c(gl, gl))
  csum2d <- het.1 + het.2 + hom.1 + hom.2
  het.md <- DmvFunc(c(e.mu[1], e.mu[2]), c(e.mu[1], e.mu[2]), 
                    sigma, invert.sigma, nu) * scale * ngt.t[2]
  hom.md <- DmvFunc(c(mu0.tx, e.mu[3]), c(mu0.tx, e.mu[3]), 
                    hom.sigma, invert.hom.sigma, nu) * scale * ngt.t[1]
  level_scales <- c(0.99, 0.95, 0.8, 0.5, 0.2)
  out.pch <- rep(".", ncol(d))
  if (any(is.na(marker.gt.pal))) {
    pal <- HAPSEG:::GetHapsegMarkerPal()
  } else {
    pal <- marker.gt.pal
  }
  cols <- rev(pal(1000))
  state.mat <- matrix(c(0, 1, 2, 3), ncol = 4, nrow = nrow(snp.clust.p), 
                      byrow = TRUE)
  col.ix <- round(rowSums(state.mat * snp.clust.p[, c(3, 2, 
                                                      4, 1)])/3 * 999, 0) + 1
  pcols <- cols[col.ix]
  x <- d[1, ]
  y <- d[2, ]
  ix <- (x < max) & (x > min) & (y < max) & (y > min)
  if (sum(ix) > 5) {
    plot(0, type = "n", main = "", xlim = c(min, max), ylim = c(min, 
                                                                max), xlab = "A allele copy-ratio", ylab = "B allele copy-ratio", 
         bty = "n")
  } else {
    smoothScatter(x, y, main = data.name, xlab = "A allele copy-ratio", 
                  ylab = "B allele copy-ratio", bty = "n")
  }
  points(d[1, ], d[2, ], col = pcols, pch = out.pch, cex = cex)
  cl.het = contourLines(x = x2d, y = x2d, z = het.1, levels = level_scales * 
                          het.md)
  cl.hom <- contourLines(x = x2d, y = x2d, z = hom.1, levels = level_scales * 
                           hom.md)
  for (i in seq_along(cl.het)) {
    txi.x <- HAPSEG:::HTxInv(cl.het[[i]][["x"]], sigma.nu, sigma.eta)
    txi.y <- HAPSEG:::HTxInv(cl.het[[i]][["y"]], sigma.nu, sigma.eta)
    lines(txi.x, txi.y, col = het.col)
    txi.x <- HAPSEG:::HTxInv(cl.het[[i]][["y"]], sigma.nu, sigma.eta)
    txi.y <- HAPSEG:::HTxInv(cl.het[[i]][["x"]], sigma.nu, sigma.eta)
    lines(txi.x, txi.y, col = het.col)
  }
  for (i in seq_along(cl.hom)) {
    txi.x <- HAPSEG:::HTxInv(cl.hom[[i]][["x"]], sigma.nu, sigma.eta)
    txi.y <- HAPSEG:::HTxInv(cl.hom[[i]][["y"]], sigma.nu, sigma.eta)
    lines(txi.x, txi.y, col = hom.col)
    txi.x <- HAPSEG:::HTxInv(cl.hom[[i]][["y"]], sigma.nu, sigma.eta)
    txi.y <- HAPSEG:::HTxInv(cl.hom[[i]][["x"]], sigma.nu, sigma.eta)
    lines(txi.x, txi.y, col = hom.col)
  }
  abline(v = mu0, lty = 2)
  abline(h = mu0, lty = 2)
  d.tx <- HAPSEG:::HTx(d, sigma.nu, sigma.eta)
  plot(0, type = "n", main = "", xlim = c(t.min, t.max), ylim = c(t.min, 
                                                                  t.max), xlab = "A allele copy-ratio", ylab = "B allele copy-ratio", 
       bty = "n")
  mtext("Variance-stabilizing transformation", side = 3, line = 0, 
        adj = 0, cex = par("cex") * par("cex.axis"))
  points(d.tx[1, ], d.tx[2, ], col = pcols, pch = out.pch, 
         cex = cex)
  contour(x = x2d, y = x2d, z = het.1, add = TRUE, drawlabels = FALSE, 
          levels = level_scales * het.md, col = het.col)
  contour(x = x2d, y = x2d, z = het.2, add = TRUE, drawlabels = FALSE, 
          levels = level_scales * het.md, col = het.col)
  contour(x = x2d, y = x2d, z = hom.1, add = TRUE, drawlabels = FALSE, 
          levels = level_scales * hom.md, col = hom.col)
  contour(x = x2d, y = x2d, z = hom.2, add = TRUE, drawlabels = FALSE, 
          levels = level_scales * hom.md, col = hom.col)
  abline(v = 0, lty = 2)
  abline(h = 0, lty = 2)
  HAPSEG:::PlotSnpAscn(d, snp.clust.p, HAPSEG:::HTxInv(e.mu, sigma.nu, sigma.eta), 
              "Allelic copy-ratio", xcrds = snp.annot[["pos"]]/1e+06, 
              chr = snp.annot[["chr"]], seg.wd = 1, loci.annot = merged.loci, 
              merge.prob = merge.prob, seg.log.ev = seg.log.ev, min = min, 
              max = max)
  HAPSEG:::PlotSnpAscn(HAPSEG:::HTx(d, sigma.nu, sigma.eta), snp.clust.p, e.mu, 
              y.lab = "Variance stabilized ", xcrds = snp.annot[["pos"]]/1e+06, 
              chr = snp.annot[["chr"]], seg.wd = 1, loci.annot = merged.loci, 
              merge.prob = merge.prob, seg.log.ev = seg.log.ev, min = t.min, 
              max = t.max)
}
