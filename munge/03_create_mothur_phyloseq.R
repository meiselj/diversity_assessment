### Mothur Phyloseq Object

## Metadata
shared_file <- "data/mothur/mgtst.trim.contigs.good.unique.good.filter.unique.precluster.pick.opti_mcc.unique_list.shared"

constax_file <- "data/mothur/mgtst.trim.contigs.good.unique.good.filter.unique.precluster.pick.opti_mcc.unique_list.0.03.cons.taxonomy"

tree_file <- "data/mothur/mgtst.trim.contigs.good.unique.good.filter.unique.precluster.pick.opti_mcc.unique_list.0.03.rep.phylip.tre"

## TODO - Fix tree ids to match OTU ids. 
mothur_ps <- import_mothur(mothur_shared_file = shared_file,
                           mothur_constaxonomy_file = constax_file)

## Adding sample metadata
sample_data(mothur_ps) <- column_to_rownames(mgtstMetadata, var = "sample_id")

ProjectTemplate::cache("mothur_ps")