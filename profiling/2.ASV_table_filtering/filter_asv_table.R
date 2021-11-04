
library("tidyverse")
library("phyloseq")

load("data/ps.rda")



sum(otu_table(ps))

#### 3) remove any ASVs not in Kingdom bacteria or archaea ####
# sort(unique(as.character(tax_table(ps)[, 1]))) ### Kingdom
ps_ab <- subset_taxa(ps, Kingdom %in% c("Bacteria", "Archaea"))
save(ps_ab, file = "cache/ps_ab.rda")
#sum(otu_table(ps_ab))


#### 4) remove Chloroplasts  ####
# sort(unique(as.character(tax_table(ps)[, 4]))) ### Order "Chloroplast"
ps_noC <- subset_taxa(ps_ab, Order!="Chloroplast")
save(ps_noC, file = "cache/ps_noC.rda")
#sum(otu_table(ps_noC))


#### 5) remove  remove Mitochondria  ####
# sort(unique(as.character(tax_table(ps)[, 5]))) ### family "Mitochondria"
ps_noMC <- subset_taxa(ps_noC, Family!="Mitochondria")
save(ps_noMC, file = "cache/ps_noMC.rda")
#sum(otu_table(ps_noMC))




#### 6) remove ASVs observed in less than 5 samples (5%, 3306*0.05) ####

pa <- otu_table(ps_noMC)
pa[pa > 0] <- 1

sample_counts <- sort(colSums(pa))

#length(sample_counts[sample_counts >= 166])

asvs_to_keep <- names(sample_counts[sample_counts >= 5])

ps_core <- prune_taxa(asvs_to_keep, ps_noMC)

save(ps_core, file = "cache/ps_core.rda")
#sum(otu_table(ps_core))




#### 7) remove taxa (genus or family) with less than 5 observed ASVs ####

# fill in family where genus unknown
taxtab <- as.data.frame(tax_table(ps_core))
taxtab$taxa <- ifelse(is.na(taxtab$Genus), paste0("f_",taxtab$Family), taxtab$Genus)
tax_table(ps_core) <- as.matrix(taxtab)

## get taxa frequency table
taxfreq <- plyr::count(taxtab, 'taxa') %>%
  arrange(desc(freq))

taxa_to_keep <- taxfreq %>%
  filter(freq >= 5)

removed_taxa <- taxfreq %>%
  filter(freq < 5)



asvs_to_keep <- taxtab %>%
  filter(taxa %in% taxa_to_keep$taxa)
asvs_to_keep <- rownames(asvs_to_keep)

ps_core5 <- prune_taxa(asvs_to_keep, ps_core)


save(ps_core5, file = "cache/ps_core5.rda")

sum(otu_table(ps_core5))

#load("cache/ps_core5.rda")




#### Construct phylogenetic tree from 3626 taxa ####


library("seqinr")
library("ape")

## get sequences from original table
load("data/seqtabNoC.RData")
seq <- colnames(seqtabNoC)

### make padded numbers for asvs and rename asvs (numbers instead of sequences)
asvs <- paste0("asv_", formatC(c(1:ncol(otu_table(ps))), width = 6, format = "d", flag = "0"))

### save ASVs ids and sequences for later
names(asvs) <- seq

sequences <- asvs
save(sequences, file = "cache/sequences.rda")

#find sequences of filtered ps object (the ones to keep)
headers <- asvs[asvs %in% taxa_names(ps_core5)]
seqs <- names(headers)

### write sequences to make tree from to fasta file:
write.fasta(as.list(seqs), headers, "cache/unaligned_ps_core5.fasta" , open = "w", nbchar = 60, as.string = FALSE)



### noe copy fasta file to HCC to run mafft, fasttree

" Shell scripts

module load mafft/7.407
mafft unaligned_ps_core5.fasta > aligned_ps_core5.fasta

module load fasttree/2.1
FastTree -gtr -nt aligned_ps_core5.fasta > ps_core5.tre

"



# combine ps object with newly generated tree

# attach new tree file to ps object
tre <- read.tree("cache/ps_core5.tre")
ps_core5t <- merge_phyloseq(ps_core5, tre)


save(ps_core5t, file = "cache/ps_core5t.rda")


