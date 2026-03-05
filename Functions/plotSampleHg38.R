plotSampleHg38 <- function (vcf, cn, sv = NULL, title = "", regions = NULL, ylim.cn = c(0, 5), 
                            layout.height = c(4, 1.2, 3.5), y.sv = ylim.cn[2] - 1, UCSC = F) {
  if (UCSC) {
    chrOffset = c(0, cumsum(c(1.025*seqlengths(BSgenome.Hsapiens.UCSC.hg38)[1], as.numeric(seqlengths(BSgenome.Hsapiens.UCSC.hg38))[1:22])))
    names(chrOffset) =  c(1:22, "X", "Y")
    seqlevels(vcf) = gsub("chr", '', seqlevels(vcf))
    seqlevels(cn) = gsub("chr", '', seqlevels(cn))
    if (is.null(regions)) 
      regions <- GRanges(seqnames = c(1:22, "X", "Y"), IRanges(rep(1, 24), seqlengths(BSgenome.Hsapiens.UCSC.hg38)[1:24]))
  } else {
    chrOffset = c(0, cumsum(c(1.025*seqlengths(BSgenome.Hsapiens.NCBI.GRCh38)[1], as.numeric(seqlengths(BSgenome.Hsapiens.NCBI.GRCh38))[1:22])))
    names(chrOffset) = c(1:22, "X", "Y")
    if (is.null(regions)) 
      regions <- GRanges(seqnames = c(1:22, "X", "Y"), IRanges(rep(1, 24), seqlengths(BSgenome.Hsapiens.NCBI.GRCh38)[1:24]))
  }
  p <- par()
  layout(matrix(1:3, ncol = 1), height = layout.height)
  par(mar = c(0.5, 3, 0.5, 0.5), mgp = c(2, 0.25, 0), bty = "L", 
      las = 2, tcl = -0.25, cex = 1)
  xlim = c(min(chrOffset[as.character(seqnames(regions))] + 
                 start(regions)), max(chrOffset[as.character(seqnames(regions))] + 
                                        end(regions)))
  bbb <- cn[cn %over% regions]
  cat('plot VCF', fill = T)
  MutationTimeR:::.plotVcf(vcf[vcf %over% regions], bbb, legend = FALSE, col.grid = "white", 
           xaxt = FALSE, cex = 0.33, xlim = xlim)
  mtext(line = -1, side = 3, title, las = 1)
  cat('plot BB', fill = T)
  MutationTimeR:::.plotBB(bbb, ylim = ylim.cn, legend = FALSE, type = "bar", 
          col.grid = "white", col = c("lightgrey", "darkgrey"), 
          xaxt = FALSE, xlim = xlim)
  tryCatch({
    par(xpd = NA)
    if (!is.null(sv)) 
      MutationTimeR:::.plotSv(sv, y1 = y.sv, regions = regions, add = TRUE)
    par(xpd = FALSE)
  }, error = function(x) warning(x))
  par(mar = c(3, 3, 0.5, 0.5))
  cat('plot time', fill = T)
  MutationTimeR:::.plotTiming(bbb, xlim = xlim, legend = FALSE, col.grid = NA)
  if (length(regions) == 1) 
    axis(side = 1, at = pretty(c(start(regions), end(regions))) + 
           chrOffset[as.character(seqnames(regions))], labels = sitools::f2si(pretty(c(start(regions), 
                                                                                       end(regions)))))
  if (any(!is.na(cn$time))) {
    y0 <- seq(0.005, 0.995, 0.01)
    s <- MutationTimeR:::.histBeta(cn)
    g <- colorRampPalette(RColorBrewer::brewer.pal(4, "Set1")[c(3, 
                                                                2, 4)])(100)
    segments(x0 = chrOffset["MT"], y0 = y0, x1 = chrOffset["MT"] + 
               s/max(s) * 1e+08, col = g, lend = 3)
    getMode <- function(s) {
      if (all(is.na(s))) 
        return(NA)
      w <- which.max(s)
      if (w %in% c(1, length(s))) {
        m <- which(c(0, diff(s)) > 0 & c(diff(s), 0) < 
                     0)
        if (length(m) == 0) 
          return(w)
        m <- m[which.max(s[m])]
        return(if (s[w] > 2 * s[m]) w else m)
      }
      else return(w)
    }
    abline(h = y0[getMode(s)], lty = 5)
    if ("time.2nd" %in% colnames(mcols(cn))) 
      if (any(!is.na(cn$time.2nd))) {
        s2 <- MutationTimeR:::.histBeta(cn, time = "time.2nd")
        segments(x0 = chrOffset["MT"] + s/max(s) * 1e+08, 
                 y0 = y0, x1 = chrOffset["MT"] + s/max(s) * 
                   1e+08 + s2/max(s) * 1e+08, col = paste0(g, 
                                                           "44"), lend = 3)
        abline(h = y0[getMode(s2)], lty = 3)
      }
  }
  par(p[setdiff(names(p), c("cin", "cra", "csi", "cxy", "din", 
                            "page"))])
}
