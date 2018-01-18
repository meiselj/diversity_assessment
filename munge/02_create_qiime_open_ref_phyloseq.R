### QIIME  
require(Biostrings)
require(ape)

#### QIIME Phyloseq Object

## Issue loading biom file convert to json using biom convert 
## ```
## biom convert -i otu_table_mc2_w_tax_no_pynast_failures.biom \
##              -o otu_table_mc2_w_tax_no_pynast_failures_json.biom \
##              --table-type "OTU table" --to-json
## ```
qiime_dir <- "data/qiime_open_ref" 
biom_file <- file.path(qiime_dir, "otu_table_mc2_w_tax_no_pynast_failures_json.biom")
seq_file <- file.path(qiime_dir, "rep_set.fna")
tree_file <- file.path(qiime_dir, "rep_set.tre")

qiime_open_ref_ps <- import_biom(BIOMfilename = biom_file)

## Update sample names to match metadata formated sample names
sample_names(qiime_open_ref_ps) <- sample_names(qiime_open_ref_ps) %>%  
    str_replace("^N","nist_run") %>% 
    str_replace("^J", "jhu_run") %>% 
    str_replace("(?<=[:digit:])-(?=[:digit:])","_")


## Adding sample metadata
sample_data(qiime_open_ref_ps) <- column_to_rownames(mgtstMetadata, var = "sample_id")

################ Adding seq data to phyloseq object ############################

seq_dat <- readDNAStringSet(seq_file)

###### Reformating OTU names
names(seq_dat) <- names(seq_dat) %>%
    str_replace(" .*","")

###### Filtering rep seq set
biom_otus <- taxa_names(qiime_open_ref_ps)
filtered_seq_dat <- seq_dat[names(seq_dat) %in% biom_otus]

## Defining seq slot
qiime_open_ref_ps@refseq <- filtered_seq_dat

################ Adding tree data to phyloseq object ###########################

tree_dat <- read.tree(tree_file)

## Defining tree slot
phy_tree(qiime_open_ref_ps) <- tree_dat

ProjectTemplate::cache("qiime_open_ref_ps")
