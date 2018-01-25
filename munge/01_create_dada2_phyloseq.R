### DADA2 Phyloseq Object
# Generate phlyoseq object using DADA2 pipeline results
require(Biostrings)

###################### Load DADA2 Pipeline Results #############################


## Metadata --------------------------------------------------------------------
# mgtstMetadata loaded with project `ProjectTemplate::load.project()`
# tsv file with mgtstMetadata in `data` directory
meta_df <-  column_to_rownames(mgtstMetadata, var = "sample_id")

## Count Table -----------------------------------------------------------------
otu_tbl <- readRDS("data/dada2/seqtab_nochim.rds") %>%
    otu_table(taxa_are_rows = FALSE)

## Taxa data -------------------------------------------------------------------
taxa <- readRDS("data/dada2/taxa.rds")

## Rep sequences ---------------------------------------------------------------
sv_seqs <- colnames(otu_tbl)

## Rename features -------------------------------------------------------------
colnames(otu_tbl) <- paste0("SV",1:ncol(otu_tbl))

names(sv_seqs) <- paste0("SV",1:ncol(otu_tbl))

rownames(taxa) <- names(sv_seqs)[match(sv_seqs, rownames(taxa))]


###################### Create Phyloseq Object ##################################

dada_ps <- phyloseq(otu_tbl, sample_data(meta_df), tax_table(taxa))

## Define Tree Slot ------------------------------------------------------------
# tree_dat <- readRDS("data/dada2/dada_tree_GTR.rds")
## Will need to update with final tree
tree_dat <- readRDS("data/dada2/tree_fit.rds")
phy_tree(dada_ps) <- tree_dat$tree

## Define Seq Slot -------------------------------------------------------------
dada_ps@refseq <- Biostrings::DNAStringSet(sv_seqs)

###################### Saving Phyloseq Object ##################################
saveRDS(dada_ps, "data/phyloseq_objects/dada_ps.rds")
