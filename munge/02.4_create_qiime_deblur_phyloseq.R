### QIIME
require(Biostrings)
require(ape)

#### QIIME Phyloseq Object

## Issue loading biom file convert to json using biom convert
## ```
## biom convert -i all.biom \
##              -o all_json.biom \
##              --table-type "OTU table" --to-json
## ```
deblur_dir <- "data/deblur"
biom_file <- file.path(deblur_dir, "all_w_tax.biom")
seq_file <- file.path(deblur_dir, "all.seqs.fa")
tree_file <- file.path(deblur_dir, "all.seqs.tre")

qiime_deblur_ps <- import_biom(BIOMfilename = biom_file)

## Adding sample metadata
sample_names(qiime_deblur_ps) <- sample_names(qiime_deblur_ps)  %>%
    str_replace("^N","nist_run") %>%
    str_replace("^J", "jhu_run") %>%
    str_replace("(?<=[:digit:])-(?=[:digit:])","_") %>%
    str_replace("_S.*", "")

sample_data(qiime_deblur_ps) <- column_to_rownames(mgtstMetadata, var = "sample_id")

################ Adding seq data to phyloseq object ############################

seq_dat <- readDNAStringSet(seq_file)

###### Filtering rep seq set
biom_otus <- taxa_names(qiime_deblur_ps)
filtered_seq_dat <- seq_dat[names(seq_dat) %in% biom_otus]

## Defining seq slot
qiime_deblur_ps@refseq <- filtered_seq_dat

################ Adding tree data to phyloseq object ###########################
tree_dat <- read.tree(tree_file)

## Defining tree slot
phy_tree(qiime_deblur_ps) <- tree_dat

### ################## Rename OTUs #############################################
###
### Deblur uses the sequence variant as the OTU name, replacing names with SV1,
### SV2, ... to minimize feature id length
###
taxa_names(qiime_deblur_ps) <- paste0("SV", 1:length(taxa_names(qiime_deblur_ps)))

###################### Saving Phyloseq Object ##################################
saveRDS(qiime_deblur_ps, "data/phyloseq_objects/deblur_ps.rds")
