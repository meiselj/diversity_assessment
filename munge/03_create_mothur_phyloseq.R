### Mothur Phyloseq Object
require(ape)
require(Biostrings)

############################## Define Mothur Files #############################
## Defines file names and checks if the files exists  

file_basename <- paste0("mgtst.trim.contigs.good.unique.good.filter.",
                        "unique.precluster.pick.opti_mcc.unique_list.")
file_basename <- file.path("data", "mothur", file_basename)

shared_file <- paste0(file_basename, "shared")
if (!file.exists(shared_file)) stop("File not found - ", shared_file)

rep_fasta <- paste0(file_basename, "0.03.rep.fasta")
if (!file.exists(rep_fasta)) stop("File not found - ", rep_fasta)

constax_file <- paste0(file_basename, "0.03.cons.taxonomy")
if (!file.exists(constax_file)) stop("File not found - ", constax_file)

tree_file <- paste0(file_basename, "0.03.rep.phylip.tre")
if (!file.exists(tree_file)) stop("File not found - ", tree_file)

seq_file <- paste0(file_basename, "0.03.rep.ng.fasta")
if (!file.exists(seq_file)) stop("File not found - ", seq_file)




######################### Renaming OTUs ######################################## 
#  ids for seq and tree files do not match shared and tax file
otu_id_dat <- read_lines(rep_fasta) %>% 
    ## Only need read names
    {.[grepl("^>",x = .)]} %>% 
    ## Excluding information after first '|'
    str_replace("\\|.*", "")
    

## Generating a data frame with the representative sequence ids and otu ids 
otu_id_df <- data_frame(ids = otu_id_dat) %>% 
    separate(ids, c("seq_id","otu_id"), "\t") %>% 
    mutate(seq_id = str_replace(seq_id, ">",""))

## Tree data -------------------------------------------------------------------
tree_dat <- read.tree(tree_file)
tree_dat$tip.label <- otu_id_df$otu_id[match(tree_dat$tip.label, otu_id_df$seq_id)]

## Seq data --------------------------------------------------------------------
seq_dat <- readDNAStringSet(seq_file)
names(seq_dat) <- otu_id_df$otu_id[match(names(seq_dat), otu_id_df$seq_id)]




################### Generate Phyloseq Object ###################################
mothur_ps <- import_mothur(mothur_shared_file = shared_file,
                           mothur_constaxonomy_file = constax_file)

## Adding sample metadata ------------------------------------------------------
## Modifying sample names for consistency
sample_names(mothur_ps) <- sample_names(mothur_ps) %>% 
    str_replace("(?<=[:digit:])_(?=[:LETTER:])","-")

## Defining metadata slot ------------------------------------------------------
sample_data(mothur_ps) <- column_to_rownames(mgtstMetadata, var = "sample_id")

## Defining tree slot ----------------------------------------------------------
phy_tree(mothur_ps) <- tree_dat

## Defining refseq slot --------------------------------------------------------
mothur_ps@refseq <- seq_dat

###################### Caching Phyloseq Object #################################
ProjectTemplate::cache("mothur_ps")