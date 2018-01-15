### Mothur Phyloseq Object
require(ape)
require(Biostrings)

## Mothur files
shared_file <- "data/mothur/mgtst.trim.contigs.good.unique.good.filter.unique.precluster.pick.opti_mcc.unique_list.shared"

rep_fasta <- "data/mothur/mgtst.trim.contigs.good.unique.good.filter.unique.precluster.pick.opti_mcc.unique_list.0.03.rep.fasta"

constax_file <- "data/mothur/mgtst.trim.contigs.good.unique.good.filter.unique.precluster.pick.opti_mcc.unique_list.0.03.cons.taxonomy"

tree_file <- "data/mothur/mgtst.trim.contigs.good.unique.good.filter.unique.precluster.pick.opti_mcc.unique_list.0.03.rep.phylip.tre"

seq_file <- "data/mothur/mgtst.trim.contigs.good.unique.good.filter.unique.precluster.pick.opti_mcc.unique_list.0.03.rep.ng.fasta"

## Replacing rep seq name with OTU id 
otu_id_dat <- read_lines(rep_fasta)
otu_id_dat <- otu_id_dat[grepl("^>",x = otu_id_dat)] %>% 
    str_replace("\\|.*", "")

otu_id_df <- data_frame(ids = otu_id_dat) %>% 
    separate(ids, c("seq_id","otu_id"), "\t") %>% 
    mutate(seq_id = str_replace(seq_id, ">",""))

## Tree data
tree_dat <- read.tree(tree_file)

### Check for correct label matching
## label_test <- data_frame(seq_id = tree_dat$tip.label)

tree_dat$tip.label <- otu_id_df$otu_id[match(tree_dat$tip.label, otu_id_df$seq_id)]

## label_test$otu_id <- tree_dat$tip.label
## nrow(otu_id_df) == nrow(label_test)
##nrow(otu_id_df) == nrow(left_join(otu_id_df, label_test))

## Seq data
seq_dat <- readDNAStringSet(seq_file)
names(seq_dat) <- otu_id_df$otu_id[match(names(seq_dat), otu_id_df$seq_id)]

### Generate phyloseq object
mothur_ps <- import_mothur(mothur_shared_file = shared_file,
                           mothur_constaxonomy_file = constax_file)

## Adding sample metadata
## Modifying sample names for consistency
sample_names(mothur_ps) <- sample_names(mothur_ps) %>% 
    str_replace("(?<=[:digit:])_(?=[:LETTER:])","-")

sample_data(mothur_ps) <- column_to_rownames(mgtstMetadata, var = "sample_id")

## Adding Tree Data
phy_tree(mothur_ps) <- tree_dat

## Adding seq data 
mothur_ps@refseq <- seq_dat

ProjectTemplate::cache("mothur_ps")