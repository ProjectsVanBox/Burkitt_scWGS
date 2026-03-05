detach_package <- function(pkg) {
  pkg_name <- paste0("package:", pkg)
  if (pkg_name %in% search()) {
    detach(pkg_name, unload = TRUE, character.only = TRUE)
    message(paste("Package", pkg, "detached successfully."))
  } else {
    message(paste("Package", pkg, "is not attached."))
  }
}

data_summary <- function(data, varname, groupnames){
  require(plyr)
  summary_func <- function(x, col){
    c(mean = mean(x[[col]], na.rm=TRUE),
      sd = sd(x[[col]], na.rm=TRUE))
  }
  data_sum<-ddply(data, groupnames, .fun=summary_func,
                  varname)
  data_sum <- rename(data_sum, c("mean" = varname))
  return(data_sum)
}

## Grep genes
grep_genes <- function(patterns_vector, dataset_genes = rownames(srat)) {
  unique_patterns_vec <- paste0("--", patterns_vector, "$")

  output <- vector(mode = "character", length = 0)

  for (i in 1:length(unique_patterns_vec)){
    output_i <- grep(unique_patterns_vec[i], dataset_genes)
    output <- c(output, dataset_genes[output_i])
    output_g <- gsub("ENS...............--", "", output)
  }
  output
}
merge_sparse = function(mat_list) {
  A <- mat_list[[1]]
  for (i in 2:length(mat_list)){ #i indexes of matrices
    B <- mat_list[[i]]
    # finding what's missing
    misA <- rownames(B)[!rownames(B) %in% rownames(A)]
    misB <- rownames(A)[!rownames(A) %in% rownames(B)]

    # adding missing rows to initial matrices
    misAlm <- matrix(numeric(length(misA) * ncol(A)), ncol = ncol(A))
    rownames(misAlm) <- misA
    misBlm <- matrix(numeric(length(misB) * ncol(B)), ncol = ncol(B))
    rownames(misBlm) <- misB

    An <- rbind(A, misAlm)
    Bn <- rbind(B, misBlm)

    # order them the same way
    Bn <- Bn[rownames(An), ]

    # final bind
    A = cbind(An, Bn)
    print(c(length(mat_list), i))
  }
  A
}

steps1 <- c("#FFFFFF64", "grey90", "grey80", "grey70", "grey60",'#4575b4','#74add1','#abd9e9','#e0f3f8','#ffffbf','#fee090','#fdae61','#f46d43', '#d73027') #first is transparent white, so that you can overlay it
steps2 <- c("white", "white", "grey90", "grey80", "grey70",'#4575b4','#74add1','#abd9e9','#e0f3f8','#ffffbf','#fee090','#fdae61','#f46d43', '#d73027')
steps3 <- c('#a50026','#d73027','#f46d43','#fdae61','#fee090','#ffffbf','#e0f3f8','#abd9e9','#74add1','#4575b4','#313695')
steps4 <- c("#607ab6", "#52a82f", "#ad51e6", "#728b2f", "#6236b1", "#518549", "#d147b6",
            "#379982", "#d43e56", "#6776dd", "#d45927", "#9958a0", "#ab7b31","#b54e75","#a25a40") #http://tools.medialab.sciences-po.fr/iwanthue/
steps5 <- c("#e48a9c","#66d1bf","#e48a3f","#71afe3","#d0af78","#a0b7cf","#bdc0a1")
steps6 <- c("grey60", "grey20", '#4575b4','#74add1','#abd9e9','#e0f3f8','#ffffbf','#fee090','#fdae61','#f46d43', '#d73027')
steps.linear <- c("white", "grey90", "grey60", '#74add1','#abd9e9','#ffffbf','#fee090','#fdae61','#f46d43', '#d73027')
steps.logscale2 <- c("white", "grey90", "grey60", '#4575b4','#74add1','#abd9e9','#fee090','#fdae61','#f46d43', '#d73027')
steps.logscale2 <- c("grey90", "grey80", "grey60", '#4575b4','#74add1','#abd9e9','#fee090','#fdae61','#f46d43', '#d73027')
steps.logscale3 <- c("grey25", '#4575b4','#74add1','#abd9e9','#fee090','#fdae61','#f46d43', '#d73027')
steps.logscale4 <- c("grey30", "grey40", '#4575b4','#74add1','#abd9e9','#fee090','#fdae61','#f46d43', '#d73027')
removeDoublets <- function(frame){
  channels <- colnames(frame)
  chan <- channels[which(grepl("FS.*W|SS.*W", channels, ignore.case = TRUE))][1]
  if (!is.na(chan)){
    doublet.gate <- deGate(frame, chan)
    doublet.gate <- max(doublet.gate, quantile(exprs(frame)[, chan],
                                               0.95, na.rm=TRUE))
    frame <- frame[which(exprs(frame)[, chan] < doublet.gate)]
  } else {
    chans <- channels[which(grepl("FS.*H|FS.*A", channels, ignore.case = TRUE))]
    if (length(chans) == 2){
      frame <- getflowFrame(flowDensity(frame, channels = chans, scale = 0.9,
                                        position = c(TRUE, TRUE), gates = c(-Inf, -Inf), ellip.gate = TRUE))
    } else {
      warning("No scatter width channel found, only one of height and area available. Doublets not removed.")
    }
  }
  return (frame)
}

## Distinct color pallet
distcols <- c('#e6194b', '#3cb44b', '#ffe119', '#4363d8', '#f58231', '#911eb4',
              '#46f0f0', '#f032e6', '#bcf60c', '#fabebe', '#008080', '#e6beff',
              '#9a6324', '#fffac8', '#800000', '#aaffc3', '#808000', '#ffd8b1',
              '#000075', '#808080', '#b8ccd1', '#000000')
distc50 <- c('#808080', '#c0c0c0', '#2f4f4f', '#556b2f', '#8b4513', '#2e8b57', '#7f0000', '#191970',
       '#006400', '#808000', '#5f9ea0', '#b8860b', '#4682b4', '#d2691e', '#9acd32', '#cd5c5c',
       '#00008b', '#32cd32', '#8fbc8f', '#8b008b', '#d2b48c', '#ff4500', '#ffa500', '#ffd700',
       '#0000cd', '#00ff00', '#9400d3', '#00ff7f', '#4169e1', '#dc143c', '#00ffff', '#00bfff',
       '#adff2f', '#ff6347', '#ff00ff', '#db7093', '#f0e68c', '#ffff54', '#6495ed', '#dda0dd',
       '#87ceeb', '#ff1493', '#ffa07a', '#ee82ee', '#98fb98', '#7fffd4', '#ff69b4', '#fffacd', '#e0ffff', '#ffb6c1')
distc31 <- c('#ff0000', '#7f2200', '#ffe1bf', '#474d00', '#ccff00', '#00593c', '#00ccff', '#002e73', '#220033',
             '#ff00aa', '#ff4059', '#330000', '#ffa280', '#8c5e00', '#fbffbf', '#71b359', '#00ffcc', '#005580',
             '#0044ff', '#cc00ff', '#7f0044', '#806460', '#ff8800', '#ffcc00', '#64664d', '#00ff22', '#7ca69d',
             '#0088ff', '#d0bfff', '#9553a6', '#ffbfd9')
gradcol <- c(gplots::colorpanel(n = 50, low = '#2f2bad', mid = '#93cfe2', high = '#fffaaa'),
             gplots::colorpanel(n = 50, low = '#fffaaa', mid = "#ffad66", high = '#d60000'))

## highlight empty and UNK cells
VlnPlots_overview <- function(srat, ngene_l, mito_cut, txpt_l, raw = F) {
  plots <- list()
  plots[['pxpt']] <- VlnPlot(srat, "n_txpt", group.by = "factor")
  plots[['txpt2']] <- VlnPlot(srat, "log2_n_txpt", group.by = "factor") +
    geom_hline(yintercept = log2(txpt_l), col = "gray50", lty = 2)
  plots[['txpt_rel']] <- VlnPlot(srat, "M_mintotal", group.by = "factor") +
    geom_hline(yintercept = 0, col="gray50", lty = 2)

  plots[['ngene']] <- VlnPlot(srat, "n_gene", group.by = "factor") +
    geom_hline(yintercept = c(ngene_l), col="gray50", lty = 2)
  plots[['ngene2']] <- VlnPlot(srat, "log2_n_gene", group.by = "factor") +
    geom_hline(yintercept = log2(ngene_l), col="gray50", lty = 2)

  plots[['mito']] <- VlnPlot(srat, "percent_mito", group.by = "factor") +
    geom_hline(yintercept = mito_cut, col="gray50", lty = 2)
  plots <- lapply(plots, function(plt) plt + theme(axis.title.x = element_blank()))
  if (raw) plots else plot_grid(plotlist = plots)
}

scatter_overview <- function(srat, mito_cut, ngene_l, txpt_l) {
  highl <- rep("normal", ncol(srat))
  names(highl) <- colnames(srat)
  highl[grep("UNK", srat@meta.data$well_id)] <- "UNK"
  highl[names(highl) %in% rownames(meta)[meta$factor == 'empty']] <- "empty"
  hl <- data.frame(highl)
  hl <- cbind.data.frame(highl, srat@meta.data[names(highl), c("n_txpt", "log2_n_txpt", "M_mintotal", "n_gene", "log2_n_gene", "percent_mito")])

  plot1 <- FeatureScatter(srat, feature1 = "n_txpt", feature2 = "percent_mito") +
    geom_hline(yintercept = mito_cut, color = "gray50", lty = 2) +
    ggtitle("percent_mito/n_txpt") +
    theme(legend.position = "none") +
    geom_point(data = hl[hl$highl == "empty", ], col = "orange", shape = 17, size = 2) +
    geom_point(data = hl[hl$highl == "UNK", ], col = "purple", shape = 17, size = 2)
  plot2 <- FeatureScatter(srat, feature1 = "n_txpt", feature2 = "n_gene") +
    geom_hline(yintercept = c(ngene_l), color = "gray50", lty = 2) +
    geom_vline(xintercept = txpt_l, color = "gray50", lty = 2) +
    ggtitle("ngene/n_txpt") +
    geom_point(data = hl[hl$highl == "empty", ], col = "orange", shape = 17, size = 2) +
    geom_point(data = hl[hl$highl == "UNK", ], col = "purple", shape = 17, size = 2)
  plot_grid(plot1, plot2)
}
reToTags = function (retok, ...) {
  require(GO.db)
  allt = dbGetQuery(GO_dbconn(), "select  go_id, term from go_term")
  inds = grep(retok, as.character(allt[,"term"]), value=FALSE, ...)
  # will error if value set at call
  if (length(inds)>0) return(as.character(allt[inds,"go_id"]))
  stop("retok did not grep to any term")
}

dotPlot <- function(srat, cluster_by, genes, ct_order = NULL, makeItBars = FALSE, max_dotsize = 6) {
  if (is.null(ct_order)) ct_order = sort(unique(as.character(srat@meta.data[ ,cluster_by])))
  gene_order = gsub(".*--", "", genes)
  to_plot = srat@assays$RNA@data[genes, ] %>%
    as.matrix() %>%
    t() %>%
    as.data.frame %>%
    rownames_to_column('sample') %>%
    merge(srat@meta.data[ ,cluster_by, drop = FALSE] %>% `colnames<-`('cluster_by') %>% rownames_to_column('sample'), by = 'sample') %>%
    reshape2::melt() %>%
    group_by(variable, cluster_by) %>%
    dplyr::summarise(perc = sum(value > 0)/length(value),
                     mean = mean(value)) %>% #[value > 0])) %>%
    mutate(gene = gsub(".*--", "", variable),
           gene = factor(gene, levels = gene_order),
           cluster_by = factor(cluster_by, levels = ct_order),
           mean =  scale(mean))
    to_plot$mean[is.na(to_plot$mean)] = 0
    expr_range = range(to_plot$mean)
    ratio = abs(expr_range[1])/abs(expr_range[2])
    clrs = c(gplots::colorpanel(n = 100*ratio, low = '#b0b0b0', high = '#f2f2f2'),
             gplots::colorpanel(n = 100, low = '#f2f2f2', mid = '#5990f0', high = '#140033'))
    if (makeItBars) {
      ggplot(to_plot, aes(x = cluster_by, y = perc, fill = mean)) +
        geom_bar(stat = 'identity', position = 'dodge') +
        scale_fill_gradientn(colors = clrs) +
        labs(x = '', y = '') +
        theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5)) +
        facet_wrap(~gene)
    } else {
      ggplot(to_plot, aes(x = gene, y = cluster_by, size = perc, color = mean)) +
        geom_point() +
        scale_color_gradientn(colors = clrs) +
        labs(x = '', y = '') +
        theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5)) +
        scale_size(range = c(0.1, max_dotsize))
    }
}
