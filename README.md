Code and data for manuscript "Single-cell whole-genome sequencing reveals convergent evolution in Burkitt lymphoma".

## Overview
1) Raw sequencing data can be obtained from EGA (EGAS5000000168).
2) Processed data (including PTATO-filtered vcf files) can be obtained from Mendeley (Steemers, Alexander (2026), “Single-cell whole-genome sequencing reveals convergent evolution in Burkitt lymphoma”, DOI: 10.17632/phk5jhhm7d.1).
3) Supporting files needed in scripts can be found under /Data.
4) Supporting functions needed in script can be found under /Functions.

## Code structure
**01_Preprocessing**
* Obtain the SBS and INDEL mutational load (01_MutLoad_sc_SBS.R & 02_MutLoad_sc_INDELS.R)
* Perform QC steps to remove low quality samples (03_QC_step1.R & 04_QC_step2.R) *(Supp. Fig. 14a-c,e)*
* Get clonal mutations from bulk samples (05_CFF_bulk.R)

**02_MIXCR**
* Plot B-cell receptor rearrangment (BCR) output from MiXCR *(Supp. Fig. 2)*
  
**03_CellPhy**
* Plot phylogenetic tree output from CellPhy *(Fig. 2,4a-b & Supp. Fig. 3,9,10,17)*
  
**04_Drivers**
* Filter putative driver mutations and generate oncoplots for both bulk and single-cell samples *(Fig. 1b, Supp. Fig. 4)*
  
**05_MutationalTimeR**
* Time copy number gains and copy-neutral loss of heterozygosity regions using MutationalTimeR *(Supp. Fig. 7)*

**06_CNA**
* Generate circosplot for bulk samples (Circosplot_bulk.R) *(Fig. 1d)*
* Generate CNV plots for all single-cells (CNV_heatmap.R) *(Supp. Fig. 6)*
* Get CNV information for bulk samples (CNV_analysis_bulk.R)
* Get CNV information for single-cell samples (CNV_analysis_single_cells.R)

**07_dNdS**
* Perform dN/dS analysis in bulk samples (dNdS_ratio_bulks.R) *(Fig. 1c)*
* Perform dN/dS analysis in single-cell data (dNdS_ratio.R) *(Fig. 3b)*
* Perform transcriptional strand bias analysis (TCR_analysis.R) *(Supp. Fig. 12)*

**08_MutSigs**
* Perform mutational patterns using own data + publicly available data from Machado et. al 2020 and PCAWG (MutationalPatterns_PMC_Machado_pcawg.R) *(Supp. Fig. 18,19)*
* Assign probability of mutational signature causing a driver mutation (SigProb.R) *(Supp. Fig. 11)*
* Compare mutational signatures between normal and tumour cells (wt_pre_post_comparison.R) *(Fig. 4c, 5a)*
  
**09_OEratio**
* Plot the observed vs expected mutational load in different cell types by taking into account all SBS mutations (OE_ratio.R) *(Fig. 3a, Supp. Fig. 8)*
* Plot the observed vs expected mutational load in different cell types by taking into account only SBS1 and SBSblood mutations (OE_ratio_SBS1_SBSblood.R) *(Fig. 5b)*

**10_Rtreefit**
* Generate ultrametric trees for all 6 single-cell patients *(Fig. 5c, Supp. Fig. 13,20)*
* Generate correlations plots between latency and age as well as latency and mutational rate (Correlations.R) *(Fig. 5d-e)*





