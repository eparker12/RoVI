#########################################
##### Taxonomic assignment in dada2 #####
#########################################
library(dada2); packageVersion("dada2")
library(phyloseq)

# extract and import even lines from fasta file (containing sequences)
system("awk 'NR % 2 == 0' dada2/merged_repseqs.fasta > dada2/merged_repseqs.csv")
seqs = as.character(read.csv("dada2/merged_repseqs.csv",header=F)[,1])

# assign taxonomy and add species using Silva
silva = assignTaxonomy(seqs, "silva_nr_v132_train_set.fa.gz", multithread=TRUE)
silva_plus = data.frame(addSpecies(silva, "silva_species_assignment_v132.fa.gz", verbose=TRUE))
silva_plus$sequence = rownames(silva_plus)

# extract sequence identifier
system("awk 'NR % 2 == 1' dada2/merged_repseqs.fasta > dada2/merged_repseqIDs.csv")
system("sed 's/^.//g' dada2/merged_repseqIDs.csv > dada2/merged_repseqIDs_clean.csv")
seq_IDs = as.character(read.csv("dada2/merged_repseqIDs_clean.csv",header=F)[,1])
rownames(silva_plus) = seq_IDs
write.csv(silva_plus, "dada2/silva_plus.csv")

# repeat allowing multiple assignments
silva_multiple = data.frame(addSpecies(silva, "silva_species_assignment_v132.fa.gz", verbose=TRUE, allowMultiple=3))
rownames(silva_multiple) = seq_IDs
write.csv(silva_multiple, "dada2/silva_multiple.csv")

#########################
#### End of pipeline ####
#########################
