library(GenomicRanges)
library(IRanges)
library(trackViewer)

# Set working directory
setwd("~/surfdrive/Shared/pmc_vanboxtel/projects/Burkitt_github/3_Output/MYC_translocation/")

## 1) Set the viewing window: stop at end of exon 2
region <- GRanges("chr8", IRanges(127734000, 127739020))

## 2) Define exons (only exon 1 and exon 2)
features <- GRanges(
  "chr8",
  IRanges(
    start = c(127736231, 127738248),
    end   = c(127736624, 127739020),
    names = c("Exon1", "Exon2")
  )
)
features$fill <- c("#000000", "#000000")

## 3) Define SNPs + labels (same list as before)
SNP <- c(
  127735043, 127735583, 127735702, 127735743, 127735788, 
  127735815, 127735863, 127736314, 
  127736451, 127736556, 127736573, 127736588, 127736696, 127736757, 127736858, 127736922, 127737393
)

SNP_labels <- c(
  "PJBU", "PWSE", "PRN4", "POFO", "PHUI", "PIA9", "PCUU", 
  "PVA9", "PPWW", "P3G6", "P856", "P2PS", "PC1A", "POB0", 
  "PEVR", "P2RW", "PJU1"
)

sample.gr <- GRanges("chr8", IRanges(SNP, width=1, names=SNP_labels))

## 3b) Keep only SNPs inside the new viewing window (optional but recommended)
sample.gr <- subsetByOverlaps(sample.gr, region)

## 4) Apply colors (default #6C5B7B, exception for 127736696)
SNP_pos <- start(sample.gr)
sample.gr$color  <- ifelse(SNP_pos == 127736696, "#C6BFD0", "#6C5B7B")
sample.gr$border <- sample.gr$color
sample.gr$alpha  <- 1  # fully opaque

## 5) Plot
pdf("lollipop_myc.pdf", width = 7, height = 10)
lolliplot(sample.gr, features, ranges = region, cex = 1)
dev.off()

#NB: 5' chromosome location of breakpoint, breakpoint ranged from 2 to 36 length