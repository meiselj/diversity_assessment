### QIIME  

#### QIIME Phyloseq Object
qiime_dir <- "data/qiime_open_ref"
biom_file <- file.path(qiime_dir, "otu_table_mc2_w_tax_no_pynast_failures.biom")
seq_file <- file.path(qiime_dir, "rep_set.fna")
tree_file <- file.path(qiime_dir, "rep_set.tre")

qiime_ps <- import_biom(BIOMfilename = biom_file)


## Adding sample metadata
# sample_data(qiime_ps) <- column_to_rownames(mgtstMetadata, var = "sample_id")

##__NOTE__
## The three samples removed by the QIIME pipeline are no template controls.

# meta_df %>% rownames_to_column() %>%
#     filter(!(rowname %in% sample_names(qiime_ps)))

ProjectTemplate::cache("qiime_ps")
