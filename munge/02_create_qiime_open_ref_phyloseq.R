### QIIME  
require(Biostrings)
require(ape)

#### QIIME Phyloseq Object
qiime_dir <- "data/qiime_open_ref"
biom_file <- file.path(qiime_dir, "otu_table_mc2_w_tax_no_pynast_failures.biom")
seq_file <- file.path(qiime_dir, "rep_set.fna")
tree_file <- file.path(qiime_dir, "rep_set.tre")

qiime_open_ref_ps <- import_biom(BIOMfilename = biom_file)


## Adding sample metadata
# sample_data(qiime_open_ref_ps) <- column_to_rownames(mgtstMetadata, var = "sample_id")

################ Renaming OTUs #################################################
##
## Addresses issue with OTU names interpreted as numeric and character strings 
##

# taxa_names(qiime_open_ref_ps) <- paste0("OTU_", taxa_names(qiime_open_ref_ps))
biom_otus <- taxa_names(qiime_open_ref_ps)

################ Adding seq data to phyloseq object ############################

seq_dat <- readDNAStringSet(seq_file)

###### Reformating OTU names
names(seq_dat) <- names(seq_dat) %>%
    str_replace(" .*","") #%>%
    # {paste0("OTU_",.)}

###### Filtering rep seq set
filtered_seq_dat <- seq_dat[names(seq_dat) %in% biom_otus]

## Defining seq slot
qiime_open_ref_ps@refseq <- filtered_seq_dat

################ Adding tree data to phyloseq object ###########################

tree_dat <- read.tree(tree_file)

###### Reformating OTU names
# tree_dat$tip.label <- paste0("OTU_", tree_dat$tip.label)

# tree_dat$tip.label

## Defining tree slot
phy_tree(qiime_open_ref_ps) <- tree_dat

ProjectTemplate::cache("qiime_open_ref_ps")
