### DADA2 Phyloseq Object
require(Biostrings)

seqtab <- readRDS("data/dada2/seqtab_nochim.rds")
otu_tbl <- otu_table(seqtab, taxa_are_rows = FALSE)

## Rep sequences
sv_seqs <- colnames(otu_tbl)
names(sv_seqs) <- paste0("SV",1:ncol(otu_tbl))

## Rename features
colnames(otu_tbl) <- paste0("SV",1:ncol(otu_tbl))

## Taxa data
taxa <- readRDS("data/dada2/taxa.rds")
rownames(taxa) <- names(sv_seqs)[match(sv_seqs, rownames(taxa))]

## Metadata
meta_df <-  column_to_rownames(mgtstMetadata, var = "sample_id")

## Tree data 
tree_dat <- readRDS("data/dada2/dada_tree_GTR.rds")

## Seq data 
seq_dat <- readDNAStringSet("data/dada2/sv_seqs.fasta")

## Create phyloseq object
dada_ps <- phyloseq(otu_tbl,
                    sample_data(meta_df),
                    tax_table(taxa))

phy_tree(dada_ps) <- tree_dat$tree
dada_ps@refseq <- seq_dat

## Removing 0 entry samples
# dada_samples <- sample_names(dada_ps)
# dada_nonzero_sample <- dada_samples[sample_sums(dada_ps) != 0]
# dada_ps <- prune_samples(dada_nonzero_sample, dada_ps) 
ProjectTemplate::cache("dada_ps")
