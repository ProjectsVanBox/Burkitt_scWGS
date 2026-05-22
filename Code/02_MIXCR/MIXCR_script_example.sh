#!/bin/bash

#SBATCH -t 2:0:0
#SBATCH --mem=25G
#SBATCH -c 8
#SBATCH --mail-type=FAIL,END


### Load modules 
module load samtools
module load Java/1.8.0_60


### Variables 
BAM=${1} # path to bam file
OutDir=${2} # path to output directory
Receptor="/hpc/pmc_vanboxtel/tools/mixcr/Receptor_subset.bed" # path to receptor bed file
PICARD=/hpc/local/CentOS7/pmc_vanboxtel/bin/picard-2.24.1/picard.jar
BGZIP=/hpc/local/CentOS7/pmc_vanboxtel/bin/samtools-bcftools-htslib-1.0_x64-linux/bin/bgzip
PWD=`pwd`/
#SampleName=$(basename $BAM | sed 's|_dedup.bam||g')
SampleName=$(basename $BAM | cut -f1 -d'_' | cut -f1 -d'.')

### Create directories 
mkdir -p $OutDir/input/bam
mkdir -p $OutDir/input/fastq/Merged
mkdir -p $OutDir/Output/$SampleName 


### Create a bamslice 
samtools view -L $Receptor -bo $OutDir/input/bam/${SampleName}_dedup.bam $BAM


### Create fastq files 
cd $OutDir/input/fastq
java -Xmx60G -jar $PICARD SamToFastq \
    INPUT=$OutDir/input/bam/${SampleName}_dedup.bam \
    OUTPUT_DIR=$PWD \
    RG_TAG=ID \
    OUTPUT_PER_RG=true
### Merge fastq files 
cat $OutDir/input/fastq/${SampleName}_*_1.fastq > $OutDir/input/fastq/Merged/${SampleName}_R1.fastq
cat $OutDir/input/fastq/${SampleName}_*_2.fastq > $OutDir/input/fastq/Merged/${SampleName}_R2.fastq
cd -


### Check MIXCR for the TCR receptor
/hpc/pmc_vanboxtel/tools/mixcr4.3.2/mixcr analyze rnaseq-full-length \
    --species hsa --dna \
    $OutDir/input/fastq/Merged/${SampleName}_R1.fastq \
    $OutDir/input/fastq/Merged/${SampleName}_R2.fastq \
    $OutDir/Output/$SampleName/${SampleName}_Receptor

# Check how much can be aligned
/hpc/pmc_vanboxtel/tools/mixcr4.3.2/mixcr exportQc align \
    $OutDir/Output/$SampleName/${SampleName}_Receptor.contigs.clns \
    $OutDir/Output/$SampleName/${SampleName}_Receptor.contigs.align.pdf 

# Check which chains are used
/hpc/pmc_vanboxtel/tools/mixcr4.3.2/mixcr exportQc chainUsage \
    $OutDir/Output/$SampleName/${SampleName}_Receptor.contigs.clns \
    $OutDir/Output/$SampleName/${SampleName}_Receptor.contigs.chainUsage.pdf

# Check what the coverage is per detected element 
/hpc/pmc_vanboxtel/tools/mixcr4.3.2/mixcr exportQc coverage \
    $OutDir/Output/$SampleName/${SampleName}_Receptor.vdjca \
    $OutDir/Output/$SampleName/${SampleName}_Receptor.coverage.pdf





