### QIIME  
require(Biostrings)
require(ape)

#### QIIME Phyloseq Object

## Issue loading biom file convert to json using biom convert 
## ```
## biom convert -i otu_table.biom -o otu_table_json.biom --table-type "OTU table" --to-json
## ```
## Note seq data not included in phyloseq object due to minimize phyloseq object size
qiime_dir <- "data/qiime_de_novo" 
biom_file <- file.path(qiime_dir, "otu_table_json.biom")
tree_file <- file.path(qiime_dir, "rep_set.tre")

qiime_de_novo_ps <- import_biom(BIOMfilename = biom_file)

## Update sample names to match metadata formated sample names
sample_names(qiime_de_novo_ps) <- sample_names(qiime_de_novo_ps) %>%  
    str_replace("^N","nist_run") %>% 
    str_replace("^J", "jhu_run") %>% 
    str_replace("(?<=[:digit:])-(?=[:digit:])","_")


## Adding sample metadata
sample_data(qiime_de_novo_ps) <- column_to_rownames(mgtstMetadata, var = "sample_id")

################ Adding tree data to phyloseq object ###########################
tree_dat <- read_tree(tree_file)

## Excluding taxa from phyloseq object not in tree 
ps_taxa <- taxa_names(qiime_de_novo_ps)
tree_taxa <- tree_dat$tip.label

taxa_to_keep <- ps_taxa[ps_taxa %in% tree_taxa]

qiime_de_novo_ps <- prune_taxa(taxa_to_keep, qiime_de_novo_ps)

## Defining phyloseq object tree slot
phy_tree(qiime_de_novo_ps) <- tree_dat


###################### Caching Phyloseq Object #################################
ProjectTemplate::cache("qiime_de_novo_ps")
