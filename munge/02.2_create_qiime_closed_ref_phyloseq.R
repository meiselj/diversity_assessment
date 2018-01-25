### QIIME
require(Biostrings)
require(ape)

#### QIIME Closed Reference Phyloseq Object

## Note no representative sequences for closed reference - can get
## from greengenes reference database
qiime_dir <- "data/qiime_closed_ref"
biom_file <- file.path(qiime_dir, "otu_table.biom")
tree_file <- file.path(qiime_dir, "97_otus.tree")

qiime_closed_ref_ps <- import_biom(BIOMfilename = biom_file,
                                   treefilename = tree_file)

## Update sample names to match metadata formated sample names
sample_names(qiime_closed_ref_ps) <- sample_names(qiime_closed_ref_ps) %>%
    str_replace("^N","nist_run") %>%
    str_replace("^J", "jhu_run") %>%
    str_replace("(?<=[:digit:])-(?=[:digit:])","_")


## Adding sample metadata
sample_data(qiime_closed_ref_ps) <- column_to_rownames(mgtstMetadata, var = "sample_id")

###################### Saving Phyloseq Object ##################################
saveRDS(qiime_closed_ref_ps, "data/phyloseq_objects/qiimeClosedRef_ps.rds")
