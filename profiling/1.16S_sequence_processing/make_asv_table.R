# Data Analysis NER2020, Maize and Sorghum samples
# here I go from raw fastq sequences to an ASV table as a phyloseq object


#================================================================================================
#================================================================================================
# SCRIPT 1 – make ASV table from raw reads - make_asv.slurm
#================================================================================================
#================================================================================================

#!/bin/sh
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=8
#SBATCH --time=24:00:00
#SBATCH --mem=30gb
#SBATCH --job-name=ASV
#SBATCH --error=./LOG/ASV_%J.err
#SBATCH --output=./LOG/ASV_%J.out

module load R/3.5
Rscript R_scripts/make_asv.R


#================================================================================================
# SCRIPT 1 – make_asv.R
#================================================================================================

#================================================================================================
### load necessary packages
#================================================================================================

library("knitr")
library("BiocStyle")
library("ggplot2")
library("gridExtra")
library("dada2")
library("phyloseq")
library("DECIPHER")
library("phangorn")



### filepath to (unzipped) f and r fastq files
miseq_path <- "/work/jyanglab/mmeier/NER2/NER2_seqs"

#================================================================================================
### filter and trim fastq reads
#================================================================================================

# Sort ensures forward/reverse reads are in same order
fnFs <- sort(list.files(miseq_path, pattern="_R1_001.fastq"))
fnRs <- sort(list.files(miseq_path, pattern="_R2_001.fastq"))
# Extract sample names, assuming filenames have format: SAMPLENAME_XXX.fastq
sampleNames <- sapply(strsplit(fnFs, "_S"), `[`, 1)
# Specify the full path to the fnFs and fnRs
fnFs <- file.path(miseq_path, fnFs)
fnRs <- file.path(miseq_path, fnRs)


### make 'filtered' folder and define the filenames for the filtered fastq.gz files.
filt_path <- file.path(miseq_path, "filtered") # Place filtered files in filtered/ subdirectory
if(!file_test("-d", filt_path)) dir.create(filt_path)
filtFs <- file.path(filt_path, paste0(sampleNames, "_F_filt.fastq.gz"))
filtRs <- file.path(filt_path, paste0(sampleNames, "_R_filt.fastq.gz"))

### filter and trim using dada2
##### THIS TAKES ~1h per 100 f/r samples
out <- filterAndTrim(fnFs, filtFs, fnRs, filtRs, truncLen=c(240,200),
  maxN=0, maxEE=c(2,2), truncQ=2, rm.phix=TRUE,
  compress=TRUE, multithread=FALSE)

#================================================================================================
### infer sequence variants from forward and reverse fastq reads
### derep is fast but takes a lot of memory. used 96G for 576 samples
### mergePairs takes a few hrs probably
#================================================================================================

# dereplication
derepFs <- derepFastq(filtFs, verbose=TRUE)
derepRs <- derepFastq(filtRs, verbose=TRUE)
# Name the derep-class objects by the sample names
names(derepFs) <- sampleNames
names(derepRs) <- sampleNames

# filter sequence variants that are probably PCR or sequencing errors and not biological
errF <- learnErrors(filtFs, multithread=TRUE)
errR <- learnErrors(filtRs, multithread=TRUE)

# independent inference by sample
dadaFs <- dada(derepFs, err=errF, pool = FALSE, multithread=TRUE)
dadaRs <- dada(derepRs, err=errR, pool = FALSE, multithread=TRUE)

# merge together the inferred forward and reverse sequences
mergers <- mergePairs(dadaFs, derepFs, dadaRs, derepRs)
save(mergers, file = "mergers.RData")




#================================================================================================
### merge fastq reads
#================================================================================================

library("knitr")
library("BiocStyle")
library("ggplot2")
library("gridExtra")
library("dada2")
library("phyloseq")
library("DECIPHER")
library("phangorn")

### filepath to (unzipped) f and r fastq files
miseq_path <- "/work/jyanglab/mmeier/NER2/NER2_seqs"


# Sort ensures forward/reverse reads are in same order
fnFs <- sort(list.files(miseq_path, pattern="_R1_001.fastq"))
fnRs <- sort(list.files(miseq_path, pattern="_R2_001.fastq"))
# Extract sample names, assuming filenames have format: SAMPLENAME_XXX.fastq
sampleNames <- sapply(strsplit(fnFs, "_S"), `[`, 1)

### make 'filtered' folder and define the filenames for the filtered fastq.gz files.
filt_path <- file.path(miseq_path, "filtered") # Place filtered files in filtered/ subdirectory
filtFs <- file.path(filt_path, paste0(sampleNames, "_F_filt.fastq.gz"))
filtRs <- file.path(filt_path, paste0(sampleNames, "_R_filt.fastq.gz"))

#================================================================================================
### infer sequence variants from forward and reverse fastq reads
### derep is fast but takes a lot of memory. used 96G for 576 samples
### mergePairs takes a few hrs probably
#================================================================================================

# dereplication
derepFs <- derepFastq(filtFs, verbose=TRUE)
derepRs <- derepFastq(filtRs, verbose=TRUE)
# Name the derep-class objects by the sample names
names(derepFs) <- sampleNames
names(derepRs) <- sampleNames

save(derepFs, file = "derepFs.rda")
save(derepRs, file = "derepRs.rda")
print("derep completed")

# filter sequence variants that are probably PCR or sequencing errors and not biological
errF <- learnErrors(filtFs, multithread=TRUE)
errR <- learnErrors(filtRs, multithread=TRUE)

save(errF, file = "errF.rda")
save(errR, file = "errR.rda")
print("error learning completed")


# independent inference by sample
dadaFs <- dada(derepFs, err=errF, pool = FALSE, multithread=TRUE)
dadaRs <- dada(derepRs, err=errR, pool = FALSE, multithread=TRUE)

save(dadaFs, file = "dadaFs.rda")
save(dadaRs, file = "dadaRs.rda")
print("dada inference completed")






# merge reads



# merge together the inferred forward and reverse sequences
mergers <- mergePairs(dadaFs, derepFs, dadaRs, derepRs)
save(mergers, file = "mergers.RData")
print("mergers completed")







#================================================================================================
#================================================================================================
# SCRIPT 2 – Construct ASV table,remove chimeras, assign taxonomy
#================================================================================================
#================================================================================================

#!/bin/sh
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=8
#SBATCH --time=20:00:00
#SBATCH --mem=32gb
#SBATCH --job-name=table
#SBATCH --error=./LOG/table_%J.err
#SBATCH --output=./LOG/table_%J.out

module load R/3.5
Rscript R_scripts/table.R



#================================================================================================
# SCRIPT 2 – table.R
#================================================================================================
#================================================================================================
### load necessary packages
#================================================================================================

library("knitr")
library("BiocStyle")
library("ggplot2")
library("gridExtra")
library("dada2")
library("phyloseq")
library("DECIPHER")
library("phangorn")


load("mergers.RData")


seqtabAll <- makeSequenceTable(mergers)
save(seqtabAll, file = "seqtabAll.RData")

# remove chimeric sequences by comparing each inferred sequence to the others in the table
seqtabNoC <- removeBimeraDenovo(seqtabAll)
save(seqtabNoC, file = "seqtabNoC.RData")


### Assign taxonomy

## filepath to SILVA train set for DADA 2 https://zenodo.org/record/1172783#.XhN6UxdKh24
fastaRef <- "/work/jyanglab/mmeier/NER2/silva_nr99_v138.1_wSpecies_train_set.fa.gz"
taxTab <- assignTaxonomy(seqtabNoC, refFasta = fastaRef, multithread=TRUE)
unname(head(taxTab))
## save taxTab
save(taxTab, file = "taxTab.RData")




#================================================================================================
#================================================================================================
# SCRIPT 3 – make phyloseq object
#================================================================================================
#================================================================================================

#!/bin/sh
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1
#SBATCH --time=4:00:00
#SBATCH --mem=16gb
#SBATCH --job-name=make_ps
#SBATCH --error=./LOG/make_ps_%J.err
#SBATCH --output=./LOG/make_ps_%J.out

module load R/3.5
Rscript R_scripts/make_ps.R


#================================================================================================
# SCRIPT 3 – make_ps.R
#================================================================================================


library("phyloseq")
library("dplyr")

### ASV table
load("seqtabNoC.RData")
### taxonomy table
load("taxTab.RData")
### metadata
samdf <- read.csv("NER2_sample_data.csv", header = TRUE)
rownames(samdf) <- samdf$Sample_ID ## use sample_ID as rownames, otherwise ps won't accept it

#=======================================================================
#### make phyloseq object
#=======================================================================

ps <- phyloseq(otu_table(seqtabNoC, taxa_are_rows=FALSE),
  sample_data(samdf),
  tax_table(taxTab))

### make padded numbers for asvs
asvnum <- paste0("asv_", formatC(c(1:ncol(otu_table(ps))), width = 6, format = "d", flag = "0"))

### remember which sequence belongs to which asv number in sequences.RData" file
names(asvnum) <- taxa_names(ps)
save(asvnum, file = "sequences.rda")

## use asv numbers as index instead of sequences
taxa_names(ps) <- asvnum
save(ps, file = "ps.rda") ## save phyloseq object





