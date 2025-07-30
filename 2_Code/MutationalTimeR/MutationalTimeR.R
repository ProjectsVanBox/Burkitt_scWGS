source("~/Downloads/plotSampleHg38.R")
devtools::install_github("mg14/mg14")
devtools::install_github("gerstung-lab/MutationTimeR")

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
ref_genome <- "BSgenome.Hsapiens.UCSC.hg38"
ref_genome <- "BSgenome.Hsapiens.NCBI.GRCh38"

#### Get filenames of the vcf files that are processed
vcf_files <- list.files("~/hpc/pmc_vanboxtel/projects/Lymphoma/3_Output/biobank/vcffilter_2023_04/vcffiltervafn0",
                        pattern = "*.vcf", full.names = TRUE
)
#remove AUH and HND in future:
vcf_files <- vcf_files[-c(16,19)]
samples_names <- c()
full_names <- c()
for (vcf_file in vcf_files) {
  sample_name <- tail(strsplit(vcf_file, "/")[[1]], n=1)
  full_names <- c(full_names, substring(sample_name, 1, 11))
  threeletters <- substring(sample_name, 9, 11)
  samples_names <- c(samples_names, threeletters)
}
cnv_samps <- samples_names

cnv_files <- c()
purity_f <- c()
for (full_name in full_names){
  cnv_files <- c(cnv_files, paste0("~/hpc/pmc_vanboxtel/projects/Lymphoma/3_Output/biobank/hmftools/", full_name, "/purple/", full_name, "T.purple.cnv.somatic.tsv"))
  purity_f <- c(purity_f, paste0("~/hpc/pmc_vanboxtel/projects/Lymphoma/3_Output/biobank/hmftools/", full_name, "/purple/", full_name, "T.purple.purity.tsv"))
}

# vcf_files <- "~/hpc/pmc_vanboxtel/projects/Lymphoma/3_Output/biobank/vcffilter_2023_04/vcffiltervafn0/PMABM000GOI_PMABM000GNV_PMCRZ615MLB_WGSVAF.vcf"
# cnv_files <- "~/hpc/pmc_vanboxtel/projects/Lymphoma/3_Output/biobank/hmftools/PMABM000GOI/purple/PMABM000GOIT.purple.cnv.somatic.tsv"
# purity_f <- "~/hpc/pmc_vanboxtel/projects/Lymphoma/3_Output/biobank/hmftools/PMABM000GOI/purple/PMABM000GOIT.purple.purity.tsv"


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
  which_tum = grepl(full_name, samples)
  vcf_adj@info$t_ref_count = vcf_adj@assays@data$AD[ ,which_tum] %>% sapply('[[', 1)
  vcf_adj@info$t_alt_count = vcf_adj@assays@data$AD[ ,which_tum] %>% sapply('[[', 2)
  vcf_adj@metadata$header@header$INFO = rbind(vcf_adj@metadata$header@header$INFO,
                                              DataFrame(Number = 1, Type = 'Float',
                                                        Description = c('number of reads supporting ref',
                                                                        'number of reads supporting alt'),
                                                        row.names = c('t_ref_count', 't_alt_count')))
  gender = ifelse(cnv[nrow(cnv),'chromosome'] == 'Y', 'male', 'female') # works, no mixups due to loss-of-Y or anything
  # add chr
  seqlevels(vcf_adj, pruning.mode = 'tidy') =  paste0('chr', c(1:22, "X", "Y"))
  #seqlevels(vcf_adj, pruning.mode = 'tidy') = c(1:22, "X", "Y")
  seqlevels(cnv_gr, pruning.mode = 'tidy') = c(1:22, "X", "Y")
  #seqlevels(vcf_adj) = paste0('chr', seqlevels(vcf_adj))
  seqlengths(cnv_gr) = seqlengths(get(ref_genome))[seqlevels(cnv_gr)]
  seqlevels(cnv_gr) = paste0('chr', seqlevels(cnv_gr))
  # Run MutationTimeR
  cat('running MT...\n')
  mt = mutationTime(vcf = vcf_adj, cn = cnv_gr, purity = purity, 
                    n.boot = 20, gender = gender)
  mcols(cnv_gr) = cbind(mcols(cnv_gr), mt$T)
  vcf_adj = addMutTime(vcf_adj, mt$V)
  # return
  mts = c(mts, list(cnv_gr, vcf_adj))
  png(file=paste0("~/surfdrive/Shared/pmc_vanboxtel/projects/Burkitt_github/3_Output/MutationalTimeR/Figures/cn_noise_small_merged/MT_", samp, "_", purity, ".png"),
      width=1000, height=1250)
  plotSampleHg38(vcf_adj, cnv_gr, UCSC = T)
  dev.off()
}
# ??
# mts = lapply(seq(1, length(mts), 2), function(n) { mts[c(n, n + 1)] })
# names(mts) = cnv_samps
saveRDS(mts, file = "~/surfdrive/Shared/pmc_vanboxtel/projects/Burkitt_github/3_Output/MutationalTimeR/data/mts_noise_small_merged.rds")


## doesn't work because the seqnames are not matching (see issue on github)
plotSample(vcf_adj,cnv_gr)



## Take first the vcfs and granges
mt_vcfs <- mts[seq(from = 2, to = length(mts), by = 2)]
mt_granges <- mts[seq(from = 1, to = length(mts), by = 2)]
names(mt_vcfs) <- samples_names
names(mt_granges) <- samples_names

## Take first the vcfs and granges, exclude GQG because of ganciclovir and EQD because of low tumor purity which gave noise
mts <- mts[-c(61,62, 67, 68)]

mt_vcfs <- mts[seq(from = 2, to = length(mts), by = 2)]
mt_granges <- mts[seq(from = 1, to = length(mts), by = 2)]
names(mt_vcfs) <- samples_names[-c(31,34)]
names(mt_granges) <- samples_names[-c(31,34)]

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
ggsave("~/surfdrive/Shared/pmc_vanboxtel/projects/Lymphoma/3_Output/Biobank/Report_figures/MT_timing_width_CNA_filtering.pdf", width = 5, height = 4) 

ggplot(df_mt_events, aes(x=width)) + geom_histogram() + xlim(0, 50000000)
ggplot(df_mt_events, aes(x=timing)) + geom_histogram(bins=20, color="black", fill="grey") + theme_light() + xlab('Timing') + ylab('Timed CNAs (Count)')
ggplot(df_mt_events[df_mt_events[,'sample']!='EQD',], aes(x=timing, fill=tumor_types)) + geom_histogram(bins=20, color="black", position="stack") + theme_light() + xlab('Timing') + ylab('Timed CNAs (Count)') + theme(legend.title=element_blank())
p5A <- ggplot(df_mt_events, aes(x=timing, fill=tumor_types)) + geom_histogram(bins=20, color="black", position="stack") + geom_histogram(bins=20, color="black", position="stack") + theme_light() + xlab('Timing') + ylab('Timed CNAs (Count)') + theme(legend.title=element_blank()) + facet_wrap(~tumor_types, ncol = 4)
ggsave("~/surfdrive/Shared/pmc_vanboxtel/projects/Lymphoma/3_Output/Biobank/Report_figures/MT_CNAs_timing.pdf", width = 8, height = 3.5)

ks.test(df_mt_events[df_mt_events[,'tumor_types']=='DLBCL',]$timing, punif)
ks.test(df_mt_events$timing, punif)

t.test(df_mt_events$timing)
## Average timing of recurrent events
## Recurrent events: 1q, 7p, 7q, 17q, 22p
## same thresholds are used as for the plotting of distribution of timings above (CI<0.5 & width<7Mb)
hg38 <- rCGH::hg38
chromarms <- c('1q', '7p', '7q', '17q', '22p')
recurrent_events <- data.frame(matrix(ncol=length(chromarms), nrow=length(mt_granges)), row.names = samples_names[-c(31,34)])
colnames(recurrent_events) <- chromarms 
for (n in 1:length(mt_granges)){
  for (chromarm in colnames(recurrent_events)){
    times <- c()
    chrom <- substr(chromarm, 1, nchar(chromarm)-1)
    chr_gr <- mt_granges[[n]][mt_granges[[n]]@seqnames == paste0('chr', chrom)]
    for (m in 1:length(chr_gr)){
      if (substr(chromarm, nchar(chromarm), nchar(chromarm)) == 'p' && start(ranges(chr_gr[m])) < hg38$centromerStart[as.integer(chrom)]){
        CI <- (chr_gr[m]$time.up - chr_gr[m]$time.lo)
        if (!is.na(CI) && CI < 0.5 && width(chr_gr[m])>7000000){ # only take events with a confidence interval smaller than 0.7 so that unsure (narrow) events are filtered out
          times <- c(times, chr_gr[m]$time)
        }
      }
      else if (substr(chromarm, nchar(chromarm), nchar(chromarm)) == 'q' && end(ranges(chr_gr[m])) > hg38$centromerEnd[as.integer(chrom)]){
        CI <- (chr_gr[m]$time.up - chr_gr[m]$time.lo)
        if (!is.na(CI) && CI < 0.5 && width(chr_gr[m])>7000000){ # only take events with a confidence interval smaller than 0.7 so that unsure (narrow) events are filtered out
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
p5B <- ggplot(melt_events, aes(x=chromarm, y=time)) + geom_boxplot() + geom_point(shape=1) + ylim(0,1) + theme_light() + xlab('Recurrent CNA event') + ylab('Timing')
ggsave("~/surfdrive/Shared/pmc_vanboxtel/projects/Lymphoma/3_Output/Biobank/Report_figures/MT_recurrent_CNA_timing.pdf", width = 3.5, height = 4)


## Load tumor types for sample names to split on burkitt and dlbcl later
meta <- read.delim("~/surfdrive/Shared/pmc_vanboxtel/projects/Lymphoma/1_Input/wgs_vcfs/meta_data_wgs_vcfs_edited.tsv", header=TRUE, sep='\t')

tumor_types <- vector( "character" , length(samples_names))
tt_df <- data.frame(tumor_types) 
rownames(tt_df) <- samples_names

for (row in 1:nrow(tt_df)){
  cont_name <- rownames(tt_df)[row]
  for (metarow in 1:nrow(meta)){
    full_names <- meta[metarow,'Biomaterial_Ids']
    if ((paste(substring(full_names, 9, 11)) == cont_name) || (paste(substring(full_names, 21, 23)) == cont_name)){
      tt_df[row,'tumor_types'] <- meta[metarow,'X01_Tumor_type']
    }
  }
}

# tt_df['GYD GXT',] <- "Burkitt lymphoma, NOS (Includes all variants, see also M-9826/3)"

newnames = setNames(
  c("Hodgkin lymphoma, nodular sclerosis, NOS",
    "Burkitt lymphoma, NOS (Includes all variants, see also M-9826/3)",
    "Malignant lymphoma, large B-cell, diffuse, NOS",
    "Hodgkin lymphoma, NOS",
    "Verdenking Maligniteit",
    "Hodgkin lymphoma, mixed cellularity, NOS",
    "Malignant lymphoma, non-Hodgkin, NOS",
    "Hodgkin lymphoma, nodular lymphocyte predominance",
    "Hodgkin lymphoma, lymphocyte depletion, NOS",
    "Precursor B-cell lymphoblastic lymphoma (see also M-9836/3)",
    "Mediastinal (thymic) large B-cell lymphoma (C38.3)",
    "B lymphoblastic leukemia/lymphoma, NOS",
    "Primary mediastinal large B-cell lymphoma"),
  c("other", "Burkitt", "DLBCL", 'other', 'other', 'other', 'other', 'other', 'other', 'other', 'other', 'other','other'))

for (n in seq(length(newnames))) { tt_df[,'tumor_types'][tt_df[,'tumor_types'] == unname(newnames[n])] = names(newnames)[n] }
tt_df$extra <- 1
burkitt_samps <- rownames(tt_df[tt_df$tumor_types=="Burkitt", ])
dlbcl_samps <- rownames(tt_df[tt_df$tumor_types=="DLBCL", ])



## Load driver mutations per sample in a GRangesList -> to determine if driver has a timing
driver_mutations <- read.csv("~/surfdrive/Shared/pmc_vanboxtel/projects/Lymphoma/3_Output/Biobank/2023_04/Driver_mutations_2023_04.csv")
driver_mutations[ncol(driver_mutations)+1] <- driver_mutations[,4]
#driver_mutations[,'chrom'] <- gsub("chr", "", driver_mutations[,'chrom'])
colnames(driver_mutations)[c(4,11)] <- c('start', 'end')
grl_drivers <- GRangesList()
for (sample in unique(driver_mutations[,'sample'])){
  gr_driversample <- GRanges(driver_mutations[driver_mutations[,'sample'] == sample, ])
  grl_drivers <- c(grl_drivers, GRangesList(gr_driversample))
  names(grl_drivers)[length(grl_drivers)] <- substr(sample, 1, 3)
}

## Get signatures that were NMF extracted only
# Create signatures dataframe of extracted signatures but values from cosmic (not -like to COSMIC ones)
signatures_subset <- signatures[,c(1,2,5,7,12,19,24,25)]
## Extract mutations based on timing estimation relative to cna event as well as timing of cna event, also immediately do an overlaps with drivers to determine timed drivers
grl_timing <- c()
n_cna_time_splits <- 3 # make n groups based on timing of cna event
driver_mutations$timing <- NA
for (n in 1:(2*3)){ # prepare list of granges for various timing groups
  grl_timing <- c(grl_timing, GRanges())
}
for (n in 1:length(mt_granges)){ #loop through cnas per sample
  print(n)
  samp <- names(mt_granges)[n]
  for (m in 1:length(mt_granges[[n]])){
    cna_grange <- mt_granges[[n]][m]
    cna_time <- mt_granges[[n]][m]$time
    if (!is.na(cna_time)){ # only take events with a timing, as these potentially have timed mutations
      vcf_grange <- rowRanges(mt_vcfs[n][[1]])
      vcf_hits <- queryHits(findOverlaps(vcf_grange, cna_grange)) # find the index of the mutations in this event
      for (hit in vcf_hits){
        if (info(mt_vcfs[n][[1]])$CLS[hit] == 'clonal [early]' && !is.na(info(mt_vcfs[n][[1]])$CLS[hit])){ # check timing of mutation
          for (split in 1:n_cna_time_splits){
            if (((split-1) / n_cna_time_splits) <= cna_time && cna_time <= (split / n_cna_time_splits)){
              grl_timing[[split]] <- c(grl_timing[[split]], vcf_grange[hit])
              timed_driver <- grl_drivers[[samp]][queryHits(findOverlaps(grl_drivers[[samp]], vcf_grange[hit]))] ## vcf_grange[hit] overlaps met drivers van die sample
              if (length(timed_driver)>0){
                #print(timed_driver)
                driver_mutations[driver_mutations$X == timed_driver$X, ]$timing <- split ## add timing to driver_mutations df
              }
            }
          }
        } 
        else if (info(mt_vcfs[n][[1]])$CLS[hit] == 'clonal [late]' && !is.na(info(mt_vcfs[n][[1]])$CLS[hit])){
          for (split in 1:n_cna_time_splits){
            if (((split-1) / n_cna_time_splits) <= cna_time && cna_time <= (split / n_cna_time_splits)){
              grl_timing[[n_cna_time_splits+split]] <- c(grl_timing[[n_cna_time_splits+split]], vcf_grange[hit])
              timed_driver <- grl_drivers[[samp]][queryHits(findOverlaps(grl_drivers[[samp]], vcf_grange[hit]))]
              if (length(timed_driver)>0){
                #print(timed_driver)
                driver_mutations[driver_mutations$X == timed_driver$X, ]$timing <- n_cna_time_splits+split ## vcf_grange[hit] overlaps met drivers van die sample, add timing to driver_mutations df
              }
            }
          }
        }
      }
    }
  }  
}

timed_drivers <- driver_mutations[!is.na(driver_mutations$timing), ]
write.csv(timed_drivers, "~/surfdrive/Shared/pmc_vanboxtel/projects/Lymphoma/3_Output/Biobank/2023_04/MTimeR_timed_driver_mutations.csv", row.names=F)

counts <- timed_drivers %>% 
  group_by(gene) %>% 
  summarise(count = n())
timed_drivers <- timed_drivers %>% 
  filter(gene %in% counts$gene[count > 1])

mean <- timed_drivers %>% 
  group_by( gene ) %>% 
  summarise( mean_val = mean( timing ))

meta <- read.delim("~/surfdrive/Shared/pmc_vanboxtel/projects/Lymphoma/1_Input/wgs_vcfs/meta_data_wgs_vcfs_edited.tsv", header=TRUE, sep='\t')
newnames = setNames(
  c("Hodgkin lymphoma, nodular sclerosis, NOS",
    "Burkitt lymphoma, NOS (Includes all variants, see also M-9826/3)",
    "Malignant lymphoma, large B-cell, diffuse, NOS",
    "Hodgkin lymphoma, NOS",
    "Verdenking Maligniteit",
    "Hodgkin lymphoma, mixed cellularity, NOS",
    "Malignant lymphoma, non-Hodgkin, NOS",
    "Hodgkin lymphoma, nodular lymphocyte predominance",
    "Hodgkin lymphoma, lymphocyte depletion, NOS",
    "Precursor B-cell lymphoblastic lymphoma (see also M-9836/3)",
    "Mediastinal (thymic) large B-cell lymphoma (C38.3)",
    "B lymphoblastic leukemia/lymphoma, NOS",
    "Primary mediastinal large B-cell lymphoma"),
  c("HL", "Burkitt", "DLBCL", 'HL', 'other', 'HL', 'other', 'HL', 'HL', 'PBCLL', 'DLBCL', 'B-ALL','PMLBL'))

timed_drivers$sample <- substr(timed_drivers$sample, 1, 3)
timed_drivers[, ncol(timed_drivers)+1] <- ''
colnames(timed_drivers)[ncol(timed_drivers)] <- 'tumor_types'
for (n in 1:nrow(timed_drivers)){
  cont_name <- timed_drivers[n,'sample']
  for (metarow in 1:nrow(meta)){
    full_names <- meta[metarow,'Biomaterial_Ids']
    if ((substring(full_names, 9, 11) == cont_name) || (substring(full_names, 21, 23) == cont_name)){
      timed_drivers[n,'tumor_types'] <- meta[metarow,'X01_Tumor_type']
    }
  }
}
for (n in seq(length(newnames))) { timed_drivers[,'tumor_types'][timed_drivers[,'tumor_types'] == unname(newnames[n])] = names(newnames)[n] }

p5E <- ggplot(timed_drivers[!(timed_drivers$sample) %in% c('EQD', 'GQG'), ], aes(x=timing, y=gene, colour=tumor_types)) + geom_beeswarm() + scale_x_continuous(breaks=seq(1,6,1), limits = c(1, 6)) + theme_light() + xlab('Timing') + ylab('Driver Mutation (Gene)') + labs(colour = "Type")
ggsave("~/surfdrive/Shared/pmc_vanboxtel/projects/Lymphoma/3_Output/Biobank/Report_figures/MT_driver_mutations_timing.pdf", width = 4, height = 4)

ggplot(timed_drivers, aes(x=timing, y=gene)) + geom_boxplot()
#+ scale_x_continuous("ID", labels = as.character(ID), breaks = ID)


names(grl_timing) <- c('1', '2', '3', '4', '5', '6')
saveRDS(grl_timing, file = "~/surfdrive/Shared/pmc_vanboxtel/projects/Lymphoma/3_Output/Biobank/2023_04/MutationTimeR/data/mts_grl_timing.rds")


for (n in 1:length(grl_timing)){
  genome(grl_timing[[n]]) <- 'hg38'
  seqlevels(grl_timing[[n]], pruning.mode = 'tidy') =  paste0('chr', c(1:22, "X", "Y"))
  seqlevels(grl_timing[[n]]) = paste0('chr', c(1:22, "X", "Y"))
  grl_timing[[n]] <- get_mut_type(grl_timing[[n]], type = "snv")
}
mut_mat_all <- mut_matrix(vcf_list = grl_timing, ref_genome = ref_genome) + 0.001

plot_96_profile(mut_mat_all[])
ggsave("~/surfdrive/Shared/pmc_vanboxtel/projects/Lymphoma/3_Output/Biobank/2023_04/MTimeR_time_split_profile_noise_merged.png")
strict_refit <- fit_to_signatures_strict(mut_mat_all, signatures_subset, max_delta = 0.01)
fit_res_strict <- strict_refit$fit_res
p5D <- plot_contribution(fit_res_strict$contribution,
                         coord_flip = FALSE,
                         mode = "relative"
)
ggsave("~/surfdrive/Shared/pmc_vanboxtel/projects/Lymphoma/3_Output/Biobank/Report_figures/MT_SNVs_timing_signatures.pdf", width = 4, height = 4)

ggsave("~/surfdrive/Shared/pmc_vanboxtel/projects/Lymphoma/3_Output/Biobank/2023_04/MTimeR_time_split_refit_noise_merged.png")


## Extract certain mutations and refit signatures
#signatures = get_known_signatures(incl_poss_artifacts = T) #use selected signatures only

## Take early and late over all samples <- no large difference
gr_all <- GRanges()
gr_early <- GRanges()
gr_late <- GRanges()
## for all samps (maxd=0.05)
for (n in 1:length(mt_vcfs)){
  gr_early <- c(gr_early, rowRanges(mt_vcfs[n][[1]])[which(info(mt_vcfs[n][[1]])$CLS == "clonal [early]")])
  gr_late <- c(gr_late, rowRanges(mt_vcfs[n][[1]])[which(info(mt_vcfs[n][[1]])$CLS == "clonal [late]")])
}
## for burkitt samps (maxd=0.005)
for (burkitt in burkitt_samps){
  gr_early <- c(gr_early, rowRanges(mt_vcfs[burkitt][[1]])[which(info(mt_vcfs[burkitt][[1]])$CLS == "clonal [early]")])
  gr_late <- c(gr_late, rowRanges(mt_vcfs[burkitt][[1]])[which(info(mt_vcfs[burkitt][[1]])$CLS == "clonal [late]")])
}
## for DLBCL samps
for (dlbcl in dlbcl_samps){
  gr_early <- c(gr_early, rowRanges(mt_vcfs[dlbcl][[1]])[which(info(mt_vcfs[dlbcl][[1]])$CLS == "clonal [early]")])
  gr_late <- c(gr_late, rowRanges(mt_vcfs[dlbcl][[1]])[which(info(mt_vcfs[dlbcl][[1]])$CLS == "clonal [late]")])
}


genome(gr_early) = 'hg38'
seqlevels(gr_early, pruning.mode = 'tidy') =  paste0('chr', c(1:22, "X", "Y"))
seqlevels(gr_early) = paste0('chr', c(1:22, "X", "Y"))
genome(gr_late) = 'hg38'
seqlevels(gr_late, pruning.mode = 'tidy') =  paste0('chr', c(1:22, "X", "Y"))
seqlevels(gr_late) = paste0('chr', c(1:22, "X", "Y"))

gr_early <- get_mut_type(gr_early, type = "snv")
gr_late <- get_mut_type(gr_late, type = "snv")

mut_mat_early <- mut_matrix(vcf_list = gr_early, ref_genome = ref_genome)
mut_mat_late <- mut_matrix(vcf_list = gr_late, ref_genome = ref_genome)
colnames(mut_mat_early) <- 'early'
colnames(mut_mat_late) <- 'late'
mut_mat_all <- cbind(mut_mat_early, mut_mat_late) +0.001


plot_96_profile(mut_mat_all[])

strict_refit <- fit_to_signatures_strict(mut_mat_all, signatures, max_delta = 0.01)
fit_res_strict <- strict_refit$fit_res
plot_contribution(fit_res_strict$contribution,
                  coord_flip = FALSE,
                  mode = "relative"
)
ggsave("~/surfdrive/Shared/pmc_vanboxtel/projects/Lymphoma/3_Output/Biobank/2023_04/MTimeR_dlbcl_early_late.png")



## THIS WHEN YOU ONLY TAKE ONE SAMPLE AS EXAMPLE (LATE, EARLY, ALL)
## Take CHZ as example to compare signatures for all, green and purple -> seems to be a slight difference, SBS9 and SBS84
gr_all <- rowRanges(mt_vcfs["CHZ"][[1]]) 
genome(gr_all) = 'hg38'
seqlevels(gr_all, pruning.mode = 'tidy') =  paste0('chr', c(1:22, "X", "Y"))
seqlevels(gr_all) = paste0('chr', c(1:22, "X", "Y"))
gr_early <- gr_all[which(info(mt_vcfs["CHZ"][[1]])$CLS == "clonal [early]")]
gr_late <- gr_all[which(info(mt_vcfs["CHZ"][[1]])$CLS == "clonal [late]")]
gr_all <- get_mut_type(gr_all, type = "snv")
gr_early <- get_mut_type(gr_early, type = "snv")
gr_late <- get_mut_type(gr_late, type = "snv")

mut_mat_early <- mut_matrix(vcf_list = gr_early, ref_genome = ref_genome)
mut_mat_late <- mut_matrix(vcf_list = gr_late, ref_genome = ref_genome)
mut_mat_all <- mut_matrix(vcf_list = gr_all, ref_genome = ref_genome)
colnames(mut_mat_early) <- 'early'
colnames(mut_mat_late) <- 'late'
colnames(mut_mat_all) <- 'all'
mut_mat_all <- cbind(mut_mat_early, mut_mat_late, mut_mat_all) +0.001

plot_96_profile(mut_mat_all[])
ggsave("~/surfdrive/Shared/pmc_vanboxtel/projects/Lymphoma/3_Output/Biobank/2023_04/MTimeR_CHZ_early_late_96_profile.png")

strict_refit <- fit_to_signatures_strict(mut_mat_all, signatures, max_delta = 0.01)
fit_res_strict <- strict_refit$fit_res
plot_contribution(fit_res_strict$contribution,
                  coord_flip = FALSE,
                  mode = "relative"
)
ggsave("~/surfdrive/Shared/pmc_vanboxtel/projects/Lymphoma/3_Output/Biobank/2023_04/MTimeR_CHZ_early_late_cosmic_refit.png")


# Do bootstrap refit
contri_boots <- fit_to_signatures_bootstrapped(mut_mat_all,
                                               signatures,
                                               n_boots = 100,
                                               method = "strict"
)
bootstrap_means <- data.frame(matrix(ncol=ncol(contri_boots),nrow=3))
colnames(bootstrap_means) <- names(contri_boots[1,])
rownames(bootstrap_means) <- colnames(mut_mat_all)

# Calculate the mean contribution values over all (100) bootstrap runs
for (n in 1:3){
  sample <- contri_boots[seq(n, 300, by = 3),]
  bootstrap_means[n,] <- colSums(sample)/100
}
#write.csv(bootstrap_means, "~/surfdrive/Shared/pmc_vanboxtel/projects/Lymphoma/3_Output/Biobank/2023_04/bootstrap_means_extractedcosmic_absolute_2023_04.csv", row.names=TRUE)

bootstrap_means_t <- t(bootstrap_means)
plot <- plot_contribution(bootstrap_means_t,
                          coord_flip = FALSE,
                          mode = "relative"
)
plot <- plot + theme(axis.text.x = element_text(angle = 90))
plot

## Get GRanges of mutations filtered on time annotation
rowRanges(mts[[2]])[which(info(mts[[2]])$CLS == "clonal [NA]")]




## Overlaps with hotspot mutations to determine timed hotspots
grl_hsm <- grl_hotspots ## GRangesList HotSpot Mutations
names(grl_hsm) <- samples_names
seqlevels(grl_hsm) = paste0('chr', seqlevels(grl_hsm))
genes_gr_chr <- genes_gr
seqlevels(genes_gr_chr) = paste0('chr', seqlevels(genes_gr_chr))

timed_hsms <- data.frame(matrix(ncol=6,nrow=0))
colnames(timed_hsms) <- c('sample', 'chr' , 'pos', 'gene', 'distance_kb', 'time')
## Extract mutations based on timing estimation relative to cna event as well as timing of cna event, also immediately do an overlaps with hotspot mutations to determine timed hotspots
#grl_timing <- c()
n_cna_time_splits <- 3 # make n groups based on timing of cna event
# for (n in 1:(2*3)){ # prepare list of granges for various timing groups
#   grl_timing <- c(grl_timing, GRanges())
# }
for (n in 1:length(mt_granges)){ #loop through cnas per sample
  print(n)
  samp <- names(mt_granges)[n]
  for (m in 1:length(mt_granges[[n]])){
    cna_grange <- mt_granges[[n]][m]
    cna_time <- mt_granges[[n]][m]$time
    if (!is.na(cna_time)){ # only take events with a timing, as these potentially have timed mutations
      vcf_grange <- rowRanges(mt_vcfs[n][[1]])
      vcf_hits <- queryHits(findOverlaps(vcf_grange, cna_grange)) # find the index of the mutations in this event
      for (hit in vcf_hits){
        if (info(mt_vcfs[n][[1]])$CLS[hit] == 'clonal [early]' && !is.na(info(mt_vcfs[n][[1]])$CLS[hit])){ # check timing of mutation
          for (split in 1:n_cna_time_splits){
            if (((split-1) / n_cna_time_splits) <= cna_time && cna_time <= (split / n_cna_time_splits)){
              #grl_timing[[split]] <- c(grl_timing[[split]], vcf_grange[hit])
              timed_hsm <- grl_hsm[[samp]][queryHits(findOverlaps(grl_hsm[[samp]], vcf_grange[hit]))] ## vcf_grange[hit] overlaps met drivers van die sample
              if (length(timed_hsm)>0){
                #print(timed_hsm)
                timed_hsms[nrow(timed_hsms)+1,] <- c(names(grl_hsm)[n], gsub("chr", "", strsplit(names(timed_hsm[1]), ':')[[1]][1]), start(ranges(timed_hsm[1])), '', '', split)
                gene <- paste(mcols(genes_gr_chr[queryHits(findOverlaps(genes_gr_chr, timed_hsm[1]))])$gene, sep = ' ', collapse = ' ')
                if (gene == ""){
                  timed_hsms[nrow(timed_hsms), 'gene'] <- genes_gr_chr[nearest(timed_hsm[1], genes_gr_chr)]$gene[1]
                  timed_hsms[nrow(timed_hsms), 'distance_kb'] <- round(mcols(distanceToNearest(timed_hsm[1], genes_gr_chr))$distance /1000, digits=1)
                }
                else{
                  timed_hsms[nrow(timed_hsms), 'gene'] <- gene
                  timed_hsms[nrow(timed_hsms), 'distance_kb'] <- 0.0
                }
              }
            }
          }
        } 
        else if (info(mt_vcfs[n][[1]])$CLS[hit] == 'clonal [late]' && !is.na(info(mt_vcfs[n][[1]])$CLS[hit])){
          for (split in 1:n_cna_time_splits){
            if (((split-1) / n_cna_time_splits) <= cna_time && cna_time <= (split / n_cna_time_splits)){
              #grl_timing[[n_cna_time_splits+split]] <- c(grl_timing[[n_cna_time_splits+split]], vcf_grange[hit])
              timed_hsm <- grl_hsm[[samp]][queryHits(findOverlaps(grl_hsm[[samp]], vcf_grange[hit]))] ## vcf_grange[hit] overlaps met drivers van die sample
              if (length(timed_hsm)>0){
                #print(timed_hsm)
                timed_hsms[nrow(timed_hsms)+1,] <- c(names(grl_hsm)[n], gsub("chr", "", strsplit(names(timed_hsm[1]), ':')[[1]][1]), start(ranges(timed_hsm[1])), '', '', n_cna_time_splits+split)
                gene <- paste(mcols(genes_gr_chr[queryHits(findOverlaps(genes_gr_chr, timed_hsm[1]))])$gene, sep = ' ', collapse = ' ')
                if (gene == ""){
                  timed_hsms[nrow(timed_hsms), 'gene'] <- genes_gr_chr[nearest(timed_hsm[1], genes_gr_chr)]$gene[1]
                  timed_hsms[nrow(timed_hsms), 'distance_kb'] <- round(mcols(distanceToNearest(timed_hsm[1], genes_gr_chr))$distance /1000, digits=1)
                }
                else{
                  timed_hsms[nrow(timed_hsms), 'gene'] <- gene
                  timed_hsms[nrow(timed_hsms), 'distance_kb'] <- 0.0
                }
              }
            }
          }
        }
      }
    }
  }  
}
## 360/5005 hotspot mutations timed

write.csv(timed_hsms, "~/surfdrive/Shared/pmc_vanboxtel/projects/Lymphoma/3_Output/Biobank/2023_04/MTimeR_timed_hotspot_mutations.csv", row.names=F)

meta <- read.delim("~/surfdrive/Shared/pmc_vanboxtel/projects/Lymphoma/1_Input/wgs_vcfs/meta_data_wgs_vcfs_edited.tsv", header=TRUE, sep='\t')
newnames = setNames(
  c("Hodgkin lymphoma, nodular sclerosis, NOS",
    "Burkitt lymphoma, NOS (Includes all variants, see also M-9826/3)",
    "Malignant lymphoma, large B-cell, diffuse, NOS",
    "Hodgkin lymphoma, NOS",
    "Verdenking Maligniteit",
    "Hodgkin lymphoma, mixed cellularity, NOS",
    "Malignant lymphoma, non-Hodgkin, NOS",
    "Hodgkin lymphoma, nodular lymphocyte predominance",
    "Hodgkin lymphoma, lymphocyte depletion, NOS",
    "Precursor B-cell lymphoblastic lymphoma (see also M-9836/3)",
    "Mediastinal (thymic) large B-cell lymphoma (C38.3)",
    "B lymphoblastic leukemia/lymphoma, NOS",
    "Primary mediastinal large B-cell lymphoma"),
  c("HL", "Burkitt", "DLBCL", 'HL', 'other', 'HL', 'other', 'HL', 'HL', 'PBCLL', 'DLBCL', 'B-ALL','PMLBL'))

timed_hsms[, ncol(timed_hsms)+1] <- ''
colnames(timed_hsms)[ncol(timed_hsms)] <- 'tumor_types'
for (n in 1:nrow(timed_hsms)){
  cont_name <- timed_hsms[n,'sample']
  for (metarow in 1:nrow(meta)){
    full_names <- meta[metarow,'Biomaterial_Ids']
    if ((substring(full_names, 9, 11) == cont_name) || (substring(full_names, 21, 23) == cont_name)){
      timed_hsms[n,'tumor_types'] <- meta[metarow,'X01_Tumor_type']
    }
  }
}
for (n in seq(length(newnames))) { timed_hsms[,'tumor_types'][timed_hsms[,'tumor_types'] == unname(newnames[n])] = names(newnames)[n] }

timed_hsms_ordered <- timed_hsms[order(timed_hsms$tumor_types),]
p5F <- ggplot(timed_hsms_ordered[!(timed_hsms_ordered$sample) %in% c('EQD', 'GQG'), ], aes(x=as.numeric(time), y=sample, colour=tumor_types)) + geom_count() + scale_x_continuous(breaks=seq(1,6,1), limits = c(1, 6)) + theme_light() + xlab('Timing') + ylab('Hotspot Mutation (Sample)') + scale_y_discrete(limits = unique(timed_hsms_ordered[!(timed_hsms_ordered$sample) %in% c('EQD', 'GQG'), ]$sample)) + labs(colour = "Type")
ggsave("~/surfdrive/Shared/pmc_vanboxtel/projects/Lymphoma/3_Output/Biobank/Report_figures/MT_SNVs_timing_signatures.pdf", width = 5, height = 4)
MT_hotspot_mutations_timing.pdf

timed_hsms_copy <- timed_hsms
timed_hsms <- timed_hsms_copy

timed_hsms_ordered[!(timed_hsms_ordered$sample) %in% c('EQD', 'GQG'), ] %>% filter(time > 3) %>% nrow()

mid_row <- plot_grid(p5B, NULL, p5D, nrow = 1, scale = c(0.8, 1, 0.9))
bottom_row <- plot_grid(p5E, p5F, nrow = 1)
plot_grid(p5A, mid_row, bottom_row, nrow = 3, scale = c(0.9, 1, 1))
ggsave("~/surfdrive/Shared/pmc_vanboxtel/projects/Lymphoma/3_Output/Biobank/Report_figures/Figure5ABDEF.pdf", width = 9, height = 12)
