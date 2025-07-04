library(ggplot2)
library(roll)
library(mclust)
library(tidyverse)
library(LaplacesDemon)
library(GetoptLong)
library(ggdendro)
library(ggh4x)

myCols <- c("darkgray","steelblue3","red3","lightgray")
names(myCols) <- c("neutral1","gain","loss","neutral2")

rolling_windows <- function( mydf, n=5, o=3) {
  df_roll <- data.frame()
  
  print( "ROLLING WINDOW" )
  for (i in c(1:length(chroms))) {
    chrom <- chroms[i]
    # Select chromosome
    tmp <- subset( mydf, Chromosome==chrom)
    
    tmp_roll <- data.frame( Chromosome = tmp$Chromosome,
                            Start = roll_min(tmp$start.pos,n,min_obs=o),
                            End = roll_max(tmp$end.pos,n,min_obs=o),
                            CN = roll_median(tmp$CopyNumberNorm,n,min_obs=o))
    
    tmp_roll$Pos <- round((tmp_roll$Start+tmp_roll$End)/2)
    
    #    tmp_roll <- na.omit(tmp_roll)[-c(1:2),]
    df_roll <- rbind(df_roll,tmp_roll)
  }
  
  return( df_roll )
}

estimate_mode <- function(x) {
  d <- density(x)
  d$x[which.max(d$y)]
}

maxleveltoplot <- 4
ploidy <- 2
prefix <- ""
cytoband_file <- "/hpc/pmc_vanboxtel/resources/PTATO/resources/hg38/svs/cytoBand_hg38.txt.gz"

GetoptLong(
  "input=s", "Input folder.",
  "maxleveltoplot=i", "Maximum y value of the plot.",
  "ploidy=i", "Ploidy",
  "cytoband_file=s", "Path to cytoband file",
  "prefix=s",  "Prefix "
)


if ( prefix == "") {
  prefix <- "CNAtool"
}

CYTOBAND <- read.delim(cytoband_file, header = FALSE,  col.names = c("chrom","chromStart","chromEnd","name","gieStain"))

readcounts_fname <- list.files(path=input, pattern="readcounts.filtered.1mb.txt", recursive = T, full.names = T)

heatmap_df <- NA

for ( fname in readcounts_fname ) {
  readcounts_df <- read.table(fname, header = T)
  # Normalize Read counts
  readcounts_df$medianReadCountNorm <- readcounts_df$medianReadCount/mean(readcounts_df$medianReadCount)
  readcounts_df$CopyNumberNorm <- ploidy*(readcounts_df$CopyNumber/mean(readcounts_df$CopyNumber))
  #readcounts_df$Sample <- gsub(".+/(.+).readcounts.filtered.1mb.txt","\\1",fname)
  
  chroms <- as.character(unique((readcounts_df$Chromosome)))

  # Create rolling window dataframe
  df_roll <- rolling_windows( readcounts_df )
  
  # Add chromsome arm information from the cytoband file to the rolling window dataframe
  CYTOBAND$arm <- gsub("^(.).+","\\1",CYTOBAND$name)
  CYTOBAND$chrom <- gsub("^chr","",CYTOBAND$chrom)
  chrom_arms <- CYTOBAND %>% group_by( chrom, arm) %>% summarise( armStart = min(chromStart), armEnd = max(chromEnd)) %>% rename( Chromosome = chrom)
  df_roll <- df_roll %>% full_join( chrom_arms, by="Chromosome") %>% filter( Pos >= armStart, Pos <= armEnd)
  
  df_roll_merged <- data.frame()

  print( "GMM" )
  for (i in c(1:length(chroms))) {
    for ( a in c("p","q")) {
      chrom <- chroms[i]
      df_roll_tmp <- df_roll[df_roll$Chromosome == chrom & df_roll$arm == a,]
      if ( nrow(df_roll_tmp) == 0 ) {
        next
      }
      g <- length(Modes(df_roll_tmp$CN)$modes)
      
      x.gmm = Mclust(df_roll_tmp$CN, G=g)
      if ( is_null(x.gmm) ) {
        next
      }
      summary(x.gmm)
      x.gmm$G
      cn_raw <- x.gmm$parameters$mean
      cn_raw[cn_raw > 1.25 & cn_raw < 2.75] <- 2
      cn_states <- round(cn_raw)
      
      probs <- x.gmm$z
      colnames(probs) <- cn_states
      
      probs_corr <- t(rowsum(t(probs), group = colnames(probs), na.rm = T))
      
      df_roll_tmp$GMMclass <- x.gmm$classification
      df_roll_tmp$CopyNumberCorrected <- cn_states[df_roll_tmp$GMMclass]
      df_roll_tmp$RatioCorrection <- 1-apply(probs_corr, 1, max)
      # df_roll_tmp$RatioCorrection <- 1-apply(probs, 1, max)
      df_roll_tmp[df_roll_tmp$CN < df_roll_tmp$CopyNumberCorrected,]$RatioCorrection <- df_roll_tmp[df_roll_tmp$CN < df_roll_tmp$CopyNumberCorrected,]$RatioCorrection*-1
      df_roll_tmp$RatioCorrected <- df_roll_tmp$CopyNumberCorrected+df_roll_tmp$RatioCorrection
      if ((i%%2) == 0) {
        df_roll_tmp$Cols <- 'neutral1'
      } else {
        df_roll_tmp$Cols <- 'neutral2'
      }
      df_roll_merged <- rbind(df_roll_merged, df_roll_tmp)
    }
  }
  
  df_roll_merged$Chromosome <- factor(df_roll_merged$Chromosome, levels=chroms)
  df_roll_merged$Sample <- gsub(".+/(.+).readcounts.filtered.1mb.txt","\\1",fname)
  if ( length( df_roll_merged[df_roll_merged$RatioCorrected >=ploidy+0.8,]$Cols) > 0 ) {
    df_roll_merged[df_roll_merged$RatioCorrected >=ploidy+0.8,]$Cols <- "gain"
  }
  if ( length( df_roll_merged[df_roll_merged$RatioCorrected <=ploidy-0.8,]$Cols) > 0 ) {
    df_roll_merged[df_roll_merged$RatioCorrected <=ploidy-0.8,]$Cols <- "loss"
  }
  
  if ( is.null(nrow(heatmap_df) ) ) {
    heatmap_df <- df_roll_merged
  } else {
    heatmap_df <- rbind(heatmap_df, df_roll_merged)
  }
}

heatmap_df[heatmap_df$RatioCorrected >4,"RatioCorrected"] <- 4

heatmap_data_wide <- heatmap_df %>% select( Chromosome, arm, Pos, RatioCorrected, Sample ) %>% pivot_wider(names_from = c(Chromosome, arm, Pos), values_from = RatioCorrected)
heatmap_matrix <- as.matrix(heatmap_data_wide[,-1])
rownames(heatmap_matrix) <- heatmap_data_wide$Sample
hc <- hclust(dist(heatmap_matrix))
dendro_data <- dendro_data(hc)

p1 <- ggplot(heatmap_df, aes(x=Pos,y=Sample)) +
  geom_tile(aes(fill=RatioCorrected)) +
  scale_fill_gradientn(limits = c(0,4), values=c(1, .6, .6, .5, .4, 0), colours=c("darkred", "red", "white", "white","white", "blue")) +
  facet_nested( . ~ factor(Chromosome, levels = c(1:22,"X","Y")) + arm, scales = "free_x", space = "free_x") +
  scale_y_discrete(limits = dendro_data$labels$label) +
  theme_classic() +
  theme( panel.spacing = unit(0,'lines')) 

  
heatmap_df2 <- heatmap_df %>%
  mutate(group_id = cumsum(RatioCorrected != lag(RatioCorrected, default = first(RatioCorrected)) |
                           Sample != lag(Sample, default = first(Sample)) |
                           Chromosome != lag(Chromosome, default = first(Chromosome)) |
                           arm != lag(arm, default = first(arm))
                                     ))

heatmap_df2_sub <- heatmap_df2 %>% group_by( group_id) %>% summarise( Sample = first(Sample), Pos = mean(Pos), start = min(Start), end = max(End), Chromosome = first(Chromosome), arm = first(arm), RatioCorrected = first(RatioCorrected), .groups = 'drop'  )
heatmap_df2_sub$w1 <- heatmap_df2_sub$Pos-heatmap_df2_sub$start
heatmap_df2_sub$w2 <- heatmap_df2_sub$end-heatmap_df2_sub$Pos
heatmap_df2_sub$w <- (heatmap_df2_sub$w1+heatmap_df2_sub$w2)/2

heatmap_data_wide2 <- heatmap_df2_sub %>% select( Chromosome, arm, Pos, RatioCorrected, Sample ) %>% pivot_wider(names_from = c(Chromosome, arm, Pos), values_from = RatioCorrected)
heatmap_matrix2 <- as.matrix(heatmap_data_wide2[,-1])
rownames(heatmap_matrix2) <- heatmap_data_wide2$Sample
hc2 <- hclust(dist(heatmap_matrix2))
dendro_data2 <- dendro_data(hc2)

P3G6_samples <- c("P3G6GPDABC31", "PB11197-BLASC-BCELLP2F4", "PB11197-BLASC-BCELLP2B4", "PB11197-BLASC-BCELLP2E4", "PB11197-BLASC-BCELLP2D4", "PB11197-BLASC-BCELLP2C4", "P3G6GPDABC28", "P3G6GPDABC27", "PB11197-BLASC-BCELLP1P3", "PB11197-BLASC-BCELLP1C4", "PB11197-BLASC-BCELLP1L3", "PB11197-BLASC-BCELLP1J3", "P3G6GPDABC26", "PB11197-BLASC-BCELLP1K4" , "PB11197-BLASC-BCELLP1O3", "PB11197-BLASC-BCELLP1I4", "P3G6GDDABC71", "PB11197-BLASC-BCELLP1B4")


p2 <- ggplot(heatmap_df2_sub) +
  geom_tile(aes(Pos,Sample, width = w, fill=RatioCorrected)) +
  scale_fill_gradientn(limits = c(0,4), values=c(1, .6, .6, .5, .4, 0), colours=c("darkred", "red", "white", "white","white", "blue")) +
  facet_nested( . ~ factor(Chromosome, levels = c(1:22,"X","Y")) + arm, scales = "free_x", space = "free_x") +
  scale_y_discrete(limits = P3G6_samples) +
  theme_classic() +
  theme( panel.spacing = unit(0,'lines')) 


p3 <- ggplot(heatmap_df2_sub) +
  geom_tile(aes(Pos,Sample, width = w, fill=RatioCorrected)) +
  scale_fill_gradientn(limits = c(0,4), values=c(1, .6, .6, .5, .4, 0), colours=c("darkred", "red", "white", "white","white", "blue")) +
  facet_nested( . ~ factor(Chromosome, levels = c(1:22,"X","Y")) + arm, scales = "free_x", space = "free_x") +
  scale_y_discrete(limits = P3G6_samples) +
  theme_classic() +
  theme( panel.spacing = unit(0,'lines')) 

pdf(file=paste0(prefix, "_heatmap_p1.pdf"), width=10, height=4, pointsize=6, useDingbats=FALSE)
print(p1)
dev.off()

pdf(file=paste0(prefix, "_heatmap_p2.pdf"), width=10, height=4, pointsize=6, useDingbats=FALSE)
print(p2)
dev.off()

pdf(file=paste0(prefix, "_heatmap_p3.pdf"), width=10, height=4, pointsize=6, useDingbats=FALSE)
print(p3)
dev.off()

# library(cellPhyWrapperPlotting)
# library(str)
# library(ggtree)
# tree <- readRDS("3_Output/TreeBuilding_Markus/P3G6/CPW_04/TreeObject0.4.RDS")
# tree <- prepare_tree(tree)
# 
# p1 <- plot_gg_tree(tree)
# #tree <- readRDS("3_Output/TreeBuilding_Markus/P3G6/CPW_04/TreeObject0.4.RDS")
