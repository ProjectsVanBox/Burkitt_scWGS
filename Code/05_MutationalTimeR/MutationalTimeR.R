################################################################################
# Manuscript: Single-cell whole-genome sequencing reveals convergent evolution in Burkitt lymphoma
# Description: Script to run MutationalTimeR
# Author: Alexander Steemers
################################################################################

source("Functions/plotSampleHg38.R")
#devtools::install_github("mg14/mg14", force = TRUE)
#devtools::install_github("gerstung-lab/MutationTimeR", force = TRUE)

library(MutationTimeR)
library(ggplot2)
library(stringr)
library(dplyr)
library(MutationalPatterns)
library(reshape2)
library(ggbeeswarm)
library(cowplot)

library(BSgenome.Hsapiens.UCSC.hg38)
library(BSgenome.Hsapiens.NCBI.GRCh38)
ref_genome <- "BSgenome.Hsapiens.NCBI.GRCh38"

# For vcf_files take all "*.vep.SMuRF.filtered.sorted.vcf.gz" files found from Mendeley doi: 10.17632/phk5jhhm7d.1 in the Bulk samples and Diagnostic samples folder
# For cnv_files take all "*.purple.cnv.somatic.tsv" files produced after running PURPLE (https://github.com/hartwigmedical/hmftools/blob/master/purple/README.md)
# For purity_f take all "*..purple.purity.tsv" files produced after running PURPLE (https://github.com/hartwigmedical/hmftools/blob/master/purple/README.md)

samples_names <- c() # create a samples_name vector that matches the sample names

# load data
vcfs = lapply(vcf_files, function(v) { cat(v, fill = T) VariantAnnotation::readVcf(v)})
cnvs = lapply(cnv_files, read.table, sep = "\t", header = T)
purities = lapply(purity_f, read.table, sep = "\t", header = T)

names(vcfs) = names(cnvs) = names(purities) = samples_names

## try to reduce noise in purple cna data
for (n in 1:length(cnvs)){
  prevcn <- 2.0
  currentcn <- 2.0
  for (m in 1:nrow(cnvs[[n]])){
    cn <- cnvs[[n]][m,]$copyNumber
    if (cn > (0.9*prevcn) && cn < (1.1*prevcn)){ # prevent small changes from previous to switch the copy number (for example from 2.3 (2) to 2.5 (3))
      cnvs[[n]][m,]$copyNumber <- currentcn
    }
    else {
      cnvs[[n]][m,]$copyNumber <- round(cn)
      currentcn <- round(cn)
    }
    prevcn <- cn
  }
}

cnvs_list = list()
# run MT
mts = list()
for (samp in cnv_samps) {
  cat(samp, " ")
  cnv = cnvs[[samp]]
  vcf = vcfs[[samp]]
  purity = purities[[samp]][1,1]
  # samples = samples(header(vcf))
  # CNV simplify & calculate major/minor CN
  alleles = data.frame(
    copyNumber = cnv$copyNumber,
    allele = round(cnv$copyNumber * cnv$baf)
  ) %>%
    mutate(allele2 = round(cnv$copyNumber - allele),
           alls = paste0(allele, "_", allele2))
  cnv_gr = with(cnv, { GRanges(chromosome, IRanges(start, end)) })
  mcols(cnv_gr) = DataFrame(alleles)
  # simplify: merge connecting regions that have the same CN profile
  cnv_gr_split = split(cnv_gr, cnv_gr$alls)
  cnv_gr_simp = lapply(names(cnv_gr_split), function(cn1cn2) {
    cn = GenomicRanges::reduce(cnv_gr_split[[cn1cn2]])
    mcols(cn) = DataFrame(allele = as.numeric(gsub("_.*", "", cn1cn2)),
                          allele2 = as.numeric(gsub(".*_", "", cn1cn2)),
                          copyNumber = strsplit(cn1cn2, "_") %>% sapply(., function(x) sum(as.numeric(x))))
    cn
  }) %>% 
    do.call(c, .) %>%
    sort
  for (n in 1:length(cnv_gr_simp)){
    if (width(cnv_gr_simp[n])<1000000){ #1Mb
      if (n==1){
        cnv_gr_simp[n]$allele <- cnv_gr_simp[n+1]$allele
        cnv_gr_simp[n]$allele2 <- cnv_gr_simp[n+1]$allele2
        cnv_gr_simp[n]$copyNumber <- cnv_gr_simp[n]$allele + cnv_gr_simp[n]$allele2
      }
      else if (n==length(cnv_gr_simp)){
        cnv_gr_simp[n]$allele <- cnv_gr_simp[n-1]$allele
        cnv_gr_simp[n]$allele2 <- cnv_gr_simp[n-1]$allele2
        cnv_gr_simp[n]$copyNumber <- cnv_gr_simp[n]$allele + cnv_gr_simp[n]$allele2
      }
      else if (runValue(seqnames(cnv_gr_simp[n-1]) == seqnames(cnv_gr_simp[n+1]))){ # chromosomes around equal
        if (cnv_gr_simp[n-1]$allele == cnv_gr_simp[n+1]$allele && cnv_gr_simp[n-1]$allele2 == cnv_gr_simp[n+1]$allele2){ # both same
          cnv_gr_simp[n]$allele <- cnv_gr_simp[n-1]$allele
          cnv_gr_simp[n]$allele2 <- cnv_gr_simp[n-1]$allele2
          cnv_gr_simp[n]$copyNumber <- cnv_gr_simp[n]$allele + cnv_gr_simp[n]$allele2
        }
        else {
          if (width(cnv_gr_simp[n-1]) > width(cnv_gr_simp[n+1])){
            cnv_gr_simp[n]$allele <- cnv_gr_simp[n-1]$allele
            cnv_gr_simp[n]$allele2 <- cnv_gr_simp[n-1]$allele2
            cnv_gr_simp[n]$copyNumber <- cnv_gr_simp[n]$allele + cnv_gr_simp[n]$allele2
          }
          else {
            cnv_gr_simp[n]$allele <- cnv_gr_simp[n+1]$allele
            cnv_gr_simp[n]$allele2 <- cnv_gr_simp[n+1]$allele2
            cnv_gr_simp[n]$copyNumber <- cnv_gr_simp[n]$allele + cnv_gr_simp[n]$allele2
          }
        }
      }
      else if (runValue(seqnames(cnv_gr_simp[n-1]) == seqnames(cnv_gr_simp[n]))){
        cnv_gr_simp[n]$allele <- cnv_gr_simp[n-1]$allele
        cnv_gr_simp[n]$allele2 <- cnv_gr_simp[n-1]$allele2
        cnv_gr_simp[n]$copyNumber <- cnv_gr_simp[n]$allele + cnv_gr_simp[n]$allele2
      }
      else if (runValue(seqnames(cnv_gr_simp[n+1]) == seqnames(cnv_gr_simp[n]))){
        cnv_gr_simp[n]$allele <- cnv_gr_simp[n+1]$allele
        cnv_gr_simp[n]$allele2 <- cnv_gr_simp[n+1]$allele2
        cnv_gr_simp[n]$copyNumber <- cnv_gr_simp[n]$allele + cnv_gr_simp[n]$allele2
      }
      else {
        print(paste(samp, n)) # verify that no small events are unmerged
      }
    }
  }
  cnv_gr_simp$alls <- paste0(cnv_gr_simp$allele, '_', cnv_gr_simp$allele2)
  # simplify again: merge connecting regions that have the same CN profile
  cnv_gr_split = split(cnv_gr_simp, cnv_gr_simp$alls)
  cnv_gr_simp = lapply(names(cnv_gr_split), function(cn1cn2) {
    cn = GenomicRanges::reduce(cnv_gr_split[[cn1cn2]])
    mcols(cn) = DataFrame(allele = as.numeric(gsub("_.*", "", cn1cn2)),
                          allele2 = as.numeric(gsub(".*_", "", cn1cn2)),
                          copyNumber = strsplit(cn1cn2, "_") %>% sapply(., function(x) sum(as.numeric(x))))
    cn
  }) %>% 
    do.call(c, .) %>%
    sort
  # major/minor allele
  CN = data.frame(major_cn = rowMax(as.matrix(mcols(cnv_gr_simp)[ ,1:2])),
                  minor_cn = rowMins(as.matrix(mcols(cnv_gr_simp)[ ,1:2])),
                  tot = mcols(cnv_gr_simp)[ ,'copyNumber']) %>%
    mutate(minor_cn = replace(minor_cn, minor_cn<0, 0),
           diff = tot - major_cn - minor_cn,
           tot_int = major_cn + minor_cn)
  # process objects for MutationTimeR & run...
  main_vaf = purity
  cnv_gr = granges(cnv_gr_simp)
  cnvs_list[[samp]] = cnv_gr_simp
  #cnvs_list = append( cnvs_list, list(cnv_gr_simp))
  
  mcols(cnv_gr) = DataFrame(major_cn = CN$major_cn, minor_cn = CN$minor_cn, clonal_frequency = main_vaf)
  genome(cnv_gr) = 'hg38'
  vcf_adj = vcf
  genome(vcf_adj) = 'hg38'
  samples = samples(header(vcf)) 
  full_name <- full_names[match(samp, cnv_samps)]
  if (nchar(full_name) <= 7) {
    full_name <- full_name
  } else if ( samp == "PMCID211AAO" ) {
    full_name <- "PMABM000IBZ"
  } else {
    full_name <- gsub(".+/(PM.+).bwamem2.+","\\1",cnv_files[grepl(samp,cnv_files)])
  }
  if (nchar(full_name) <= 7) {
    which_tum = !grepl("MS", samples)
  } else {
    which_tum = grepl(full_name, samples)
  }
  vcf_adj@info$t_ref_count = vcf_adj@assays@data$AD[ ,which_tum] %>% sapply('[[', 1)
  vcf_adj@info$t_alt_count = vcf_adj@assays@data$AD[ ,which_tum] %>% sapply('[[', 2)
  vcf_adj@metadata$header@header$INFO = rbind(vcf_adj@metadata$header@header$INFO,
                                              DataFrame(Number = 1, Type = 'Float',
                                                        Description = c('number of reads supporting ref',
                                                                        'number of reads supporting alt'),
                                                        row.names = c('t_ref_count', 't_alt_count')))
  gender = ifelse(cnv[nrow(cnv),'chromosome'] == 'Y', 'male', 'female') # works, no mixups due to loss-of-Y or anything
  # add chr
  #seqlevels(vcf_adj, pruning.mode = 'tidy') =  paste0('chr', c(1:22, "X", "Y"))
  seqlevels(vcf_adj, pruning.mode = 'tidy') = c(1:22, "X", "Y")
  seqlevels(cnv_gr, pruning.mode = 'tidy') = c(1:22, "X", "Y")
  #seqlevels(vcf_adj) = paste0('chr', seqlevels(vcf_adj))
  seqlengths(cnv_gr) = seqlengths(get(ref_genome))[seqlevels(cnv_gr)]
  #seqlevels(cnv_gr) = paste0('chr', seqlevels(cnv_gr))
  # Run MutationTimeR
  cat('running MT...\n')
  mt = mutationTime(vcf = vcf_adj, cn = cnv_gr, purity = purity, 
                    n.boot = 20, gender = gender)
  mcols(cnv_gr) = cbind(mcols(cnv_gr), mt$T)
  vcf_adj = addMutTime(vcf_adj, mt$V)
  # return
  mts = c(mts, list(cnv_gr, vcf_adj))
#  cnvs_list <- append(cnvs_list, cnv_gr)
  png(file=paste0("MT_", samp, "_", purity, ".png"),
      width=1000, height=1250)
  plotSampleHg38(vcf_adj, cnv_gr, UCSC  = T)
  dev.off()
}

saveRDS(mts, file = "mts_noise_small_merged.rds")
saveRDS(cnvs_list, file = "cnvs_list.rds")


mt_vcfs <- mts[seq(from = 2, to = length(mts), by = 2)]
mt_granges <- mts[seq(from = 1, to = length(mts), by = 2)]
names(mt_vcfs) <- samples_names
names(mt_granges) <- samples_names

## Determine timing of events, and also the width of event
df_mt_events <- data.frame(matrix(ncol=3,nrow=0))
colnames(df_mt_events) <- c('sample', 'timing', 'width')
#times <- c()
for (n in 1:length(mt_granges)){
  # sum(widths(mt_granges[[n]]))
  for (m in 1:length(mt_granges[[n]])){
    CI <- (mt_granges[[n]][m]$time.up - mt_granges[[n]][m]$time.lo)
    if (!is.na(CI) && CI < 0.8){ # only take events with a confidence interval smaller than 0.5 so that unsure (narrow) events are filtered out
      if (width(mt_granges[[n]][m]) > 7000000){ #7000000
        df_mt_events[nrow(df_mt_events)+1,] <- c(names(mt_granges)[n], mt_granges[[n]][m]$time, width(mt_granges[[n]][m]))
        #times <- c(times, time)
      }
    }
  }
}

df_mt_events$width <- as.numeric(df_mt_events$width)
df_mt_events$timing <- as.numeric(df_mt_events$timing)
df_mt_events[, ncol(df_mt_events)+1] <- ''
colnames(df_mt_events)[ncol(df_mt_events)] <- 'tumor_types'
for (n in 1:nrow(df_mt_events)){
  cont_name <- df_mt_events[n,'sample']
  for (metarow in 1:nrow(meta)){
    full_names <- meta[metarow,'Biomaterial_Ids']
    if ((substring(full_names, 9, 11) == cont_name) || (substring(full_names, 21, 23) == cont_name)){
      df_mt_events[n,'tumor_types'] <- meta[metarow,'X01_Tumor_type']
    }
  }
}
for (n in seq(length(newnames))) { df_mt_events[,'tumor_types'][df_mt_events[,'tumor_types'] == unname(newnames[n])] = names(newnames)[n] }

ggplot(df_mt_events, aes(x=timing, y=width)) + geom_point(aes(color = tumor_types)) + theme_light() + xlab('Timing') + ylab('Width of CNA') + geom_hline(yintercept=7000000,linetype=2) + theme(legend.title=element_blank())

## Average timing of recurrent events
## Recurrent events: 1q, 7p, 7q, 1p
## same thresholds are used as for the plotting of distribution of timings above (CI<0.5 & width<7Mb)

hg38 <- rCGH::hg38
chromarms <- c('1p', '1q', '7p', '7q')
recurrent_events <- data.frame(matrix(ncol=length(chromarms), nrow=length(mt_granges)), row.names = samples_names[-c(31,34)])
colnames(recurrent_events) <- chromarms 

for (n in 1:length(mt_granges)){
  for (chromarm in colnames(recurrent_events)){
    times <- c()
    chrom <- substr(chromarm, 1, nchar(chromarm)-1)
    chr_gr <- mt_granges[[n]][mt_granges[[n]]@seqnames == paste0('', chrom)]
    for (m in 1:length(chr_gr)){
      if (substr(chromarm, nchar(chromarm), nchar(chromarm)) == 'p' && start(ranges(chr_gr[m])) < hg38$centromerStart[as.integer(chrom)]){
        CI <- (chr_gr[m]$time.up - chr_gr[m]$time.lo)
        if (!is.na(CI) && CI < 0.8 && width(chr_gr[m])>7000000){ # only take events with a confidence interval smaller than 0.7 so that unsure (narrow) events are filtered out
          times <- c(times, chr_gr[m]$time)
        }
      }
      else if (substr(chromarm, nchar(chromarm), nchar(chromarm)) == 'q' && end(ranges(chr_gr[m])) > hg38$centromerEnd[as.integer(chrom)]){
        CI <- (chr_gr[m]$time.up - chr_gr[m]$time.lo)
        if (!is.na(CI) && CI < 0.8 && width(chr_gr[m])>7000000){ # only take events with a confidence interval smaller than 0.7 so that unsure (narrow) events are filtered out
          times <- c(times, chr_gr[m]$time)
        } 
      }
    }
    if (length(times) > 0 && mean(as.numeric(times)) > 0.0001){
      recurrent_events[n,chromarm] <- mean(as.numeric(times))
    }
  }
}
recurrent_events <- recurrent_events[!(rownames(recurrent_events) %in% c('EQD')),]

melt_events <- melt(recurrent_events, variable.name = "chromarm", value.name = "time")
ggsave("MT_recurrent_CNA_timing.pdf", width = 3.5, height = 4)