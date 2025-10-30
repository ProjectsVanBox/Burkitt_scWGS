
.libPaths("/hpc/local/Rocky8/pmc_vanboxtel/bin/R-4.4.2/lib64/R/library/")
source("~/hpc/pmc_vanboxtel/projects/Burkitt/2_Code/MutationalTimeR/plotSampleHg38.R")
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
ref_genome <- "BSgenome.Hsapiens.UCSC.hg38"
ref_genome <- "BSgenome.Hsapiens.NCBI.GRCh38"

#patientID <- "PMCID104AAO"

#### Get filenames of the vcf files that are processed
vcf_files <- list.files(paste0("~/hpc/pmc_vanboxtel/projects/Burkitt/1_Input/Diagnostic_samples/"),
                        pattern = "*filtered.sorted.vcf.gz$", full.names = TRUE, recursive = T
)

#vcf_files <- "/hpc/pmc_vanboxtel/projects/Burkitt/1_Input/Diagnostic_samples/PMCID104AAO/SMuRF/PMCID104AAO.vep.SMuRF.filtered.sorted.vcf.gz"
vcf_files <- vcf_files[-6] #211AAO initial file but need 211AAO_2


vcf_files <- c(vcf_files, list.files(paste0("~/hpc/pmc_vanboxtel/projects/Burkitt/1_Input/ASAP_FROM_CLOUD/"),
                                     pattern = "*filtered.sorted.vcf.gz$", full.names = TRUE, recursive = T) )
vcf_files <- vcf_files[-c(13, 20, 23)] # 13 - no myc translocation ; 20 - no specific P856 sample ; 23 - duplicate of 458AAQ


#remove AUH and HND in future:
# vcf_files <- vcf_files[-c(16,19)]
samples_names <- c()
full_names <- c()
for (vcf_file in vcf_files) {
  if (grepl("P856_BM", vcf_file)) {
    sample_name <- "P856_BM"
  } else if (grepl("P856_PL", vcf_file)) {
    sample_name <- "P856_PL"
  } else if (grepl("ASAP_FROM_CLOUD", vcf_file)) {
    sample_name <- tail(strsplit(vcf_file, "/")[[1]], n=8)[1]
  } else {
    sample_name <- tail(strsplit(vcf_file, "/")[[1]], n=3)[1]
  }
  full_names <- c(full_names, sample_name)
  threeletters <- sample_name
  samples_names <- c(samples_names, threeletters)
}
cnv_samps <- samples_names

cnv_files <- c()
purity_f <- c()
for (full_name in full_names){
  if (grepl("P856_BM", full_name)) {
    full_name <- "P856/vcf_batches/batch_bulk_BM"
    cnv_files <- c(cnv_files, list.files(paste0("~/hpc/pmc_vanboxtel/projects/Burkitt/1_Input/ASAP_FROM_CLOUD/", full_name), "*.purple.cnv.somatic.tsv", recursive = T, full.names = T))
    purity_f <- c(purity_f, list.files(paste0("~/hpc/pmc_vanboxtel/projects/Burkitt/1_Input/ASAP_FROM_CLOUD/", full_name), "*.purple.purity.tsv", recursive = T, full.names = T))
  } else if (grepl("P856_PL", full_name)) {
    full_name <- "P856/vcf_batches/batch_bulk_PL"
    cnv_files <- c(cnv_files, list.files(paste0("~/hpc/pmc_vanboxtel/projects/Burkitt/1_Input/ASAP_FROM_CLOUD/", full_name), "*.purple.cnv.somatic.tsv", recursive = T, full.names = T))
    purity_f <- c(purity_f, list.files(paste0("~/hpc/pmc_vanboxtel/projects/Burkitt/1_Input/ASAP_FROM_CLOUD/", full_name), "*.purple.purity.tsv", recursive = T, full.names = T))
  } else  if (nchar(full_name) == 4) {
    cnv_files <- c(cnv_files, list.files(paste0("~/hpc/pmc_vanboxtel/projects/Burkitt/1_Input/ASAP_FROM_CLOUD/", full_name), "*.purple.cnv.somatic.tsv", recursive = T, full.names = T))
    purity_f <- c(purity_f, list.files(paste0("~/hpc/pmc_vanboxtel/projects/Burkitt/1_Input/ASAP_FROM_CLOUD/", full_name), "*.purple.purity.tsv", recursive = T, full.names = T))
  } else {
    cnv_files <- c(cnv_files, list.files(paste0("~/hpc/pmc_vanboxtel/projects/Burkitt/1_Input/Diagnostic_samples/", full_name), "*.purple.cnv.somatic.tsv", recursive = T, full.names = T))
    purity_f <- c(purity_f, list.files(paste0("~/hpc/pmc_vanboxtel/projects/Burkitt/1_Input/Diagnostic_samples/", full_name), "*.purple.purity.tsv", recursive = T, full.names = T))
  }
}

cnv_files <- cnv_files[-which(grepl("sv_old",cnv_files))]
purity_f <- purity_f[-which(grepl("sv_old",purity_f))]
cnv_files <- cnv_files[-which(grepl("sv_wrong",cnv_files))]
purity_f <- purity_f[-which(grepl("sv_wrong",purity_f))]

#cnv_files <- "/hpc/pmc_vanboxtel/projects/Burkitt/1_Input/Diagnostic_samples/PMCID104AAO/sv/somatic/gridss/PMABM000GTE.bwamem2.samtools.gatk4spark.dedup_PMABM000GTJ.bwamem2.samtools.gatk4spark.dedup/purple/PMABM000GTJ.bwamem2.samtools.gatk4spark.dedup.purple.cnv.somatic.tsv"
#purity_f <- "/hpc/pmc_vanboxtel/projects/Burkitt/1_Input/Diagnostic_samples/PMCID104AAO/sv/somatic/gridss/PMABM000GTE.bwamem2.samtools.gatk4spark.dedup_PMABM000GTJ.bwamem2.samtools.gatk4spark.dedup/purple/PMABM000GTJ.bwamem2.samtools.gatk4spark.dedup.purple.purity.tsv"


# load data
vcfs = lapply(vcf_files, function(v) {
  cat(v, fill = T)
  VariantAnnotation::readVcf(v)
})
lengths(vcfs)
cnvs = lapply(cnv_files, read.table, sep = "\t", header = T)
purities = lapply(purity_f, read.table, sep = "\t", header = T)
names(vcfs) = names(cnvs) = names(purities) = cnv_samps


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
  # Merge small events to the surrounding events
  ##
  # for (n in 1:length(cnv_gr_simp)){
  #   if (width(cnv_gr_simp[n])<1000000){ #1Mb
  #     cnv_gr_simp[n]$allele <- floor((cnv_gr_simp[n-1]$allele + cnv_gr_simp[n+1]$allele)/2)
  #     cnv_gr_simp[n]$allele2 <- floor((cnv_gr_simp[n-1]$allele2 + cnv_gr_simp[n+1]$allele2)/2)
  #     cnv_gr_simp[n]$copyNumber <- cnv_gr_simp[n]$allele + cnv_gr_simp[n]$allele2
  #   }
  # }
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
  png(file=paste0("/hpc/pmc_vanboxtel/projects/Burkitt/3_Output/MutationalTimeR/Figures/cn_noise_small_merged/MT_", samp, "_", purity, ".png"),
      width=1000, height=1250)
  plotSampleHg38(vcf_adj, cnv_gr, UCSC  = T)
  dev.off()
}
# ??
# mts = lapply(seq(1, length(mts), 2), function(n) { mts[c(n, n + 1)] })
# names(mts) = cnv_samps
saveRDS(mts, file = "/hpc/pmc_vanboxtel/projects/Burkitt/3_Output/MutationalTimeR/data/mts_noise_small_merged.rds")


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
    if (!is.na(CI) && CI < 0.5){ # only take events with a confidence interval smaller than 0.5 so that unsure (narrow) events are filtered out
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

plot(density(df_mt_events$timing))
hist(df_mt_events$timing, breaks=40) 

hist(df_mt_events$width, breaks=200, xlim = range(0, 50000000))

ggplot(df_mt_events, aes(x=timing, y=width)) + geom_point(aes(color = tumor_types)) + theme_light() + xlab('Timing') + ylab('Width of CNA') + geom_hline(yintercept=7000000,linetype=2) + theme(legend.title=element_blank())
ggsave("/hpc/pmc_vanboxtel/projects/Burkitt/3_Output/MutationalTimeR/Figures/MT_timing_width_CNA_filtering.pdf", width = 5, height = 4) 

ggplot(df_mt_events, aes(x=width)) + geom_histogram() + xlim(0, 50000000)
ggplot(df_mt_events, aes(x=timing)) + geom_histogram(bins=20, color="black", fill="grey") + theme_light() + xlab('Timing') + ylab('Timed CNAs (Count)')
ggplot(df_mt_events[df_mt_events[,'sample']!='EQD',], aes(x=timing, fill=tumor_types)) + geom_histogram(bins=20, color="black", position="stack") + theme_light() + xlab('Timing') + ylab('Timed CNAs (Count)') + theme(legend.title=element_blank())
p5A <- ggplot(df_mt_events, aes(x=timing, fill=tumor_types)) + geom_histogram(bins=20, color="black", position="stack") + geom_histogram(bins=20, color="black", position="stack") + theme_light() + xlab('Timing') + ylab('Timed CNAs (Count)') + theme(legend.title=element_blank()) + facet_wrap(~tumor_types, ncol = 4)
ggsave("/hpc/pmc_vanboxtel/projects/Burkitt/3_Output/MutationalTimeR/Figures/MT_CNAs_timing.pdf", width = 8, height = 3.5)