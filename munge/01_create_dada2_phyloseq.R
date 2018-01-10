### DADA2 Phyloseq Object
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

## Create phyloseq object
dada_ps <- phyloseq(otu_tbl,
                    sample_data(meta_df),
                    tax_table(taxa))

## Removing 0 entry samples
dada_samples <- sample_names(dada_ps)
dada_nonzero_sample <- dada_samples[sample_sums(dada_ps) != 0]
dada_ps <- prune_samples(dada_nonzero_sample, dada_ps) 

ProjectTemplate::cache("dada_ps")
