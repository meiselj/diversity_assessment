## Normalize phyloseq objects

################## Rarefaction #################################################
## Reasoning for rarefaction levels used in analysis
## 5000 in most analysis
## Weiss et al 2017? compared 2000, 5000, and 10000
## McMurdie and Holmes 2014 used the 15th percentile library size
apply_rare_methods_ps <- function(ps){
    ## Dropping refseq: reduce memory requirements, data not needed for analysis
    ps@refseq <- NULL 
    
    ## Defining rarifying levels
    rare_levels <- list(rare_2K = 2000, rare_5K = 5000, rare_10K = 10000)
    
    rare_ps <- rare_levels %>% 
        map(~rarefy_even_depth(ps, rngseed = 531, sample.size = .))
    
    ## Rarefying at 15th percentile library size
    q15 <- floor(quantile(sample_sums(ps), 0.15, na.rm = TRUE))
    rare_ps$rare_15P <- rarefy_even_depth(ps, rngseed = 531, sample.size = q15)
    
    ## Returns list of rarifyed phyloseq objects
    rare_ps
}

################# Normalization ################################################
normalize_counts <- function(method, ps) {
    if (class(ps) != "phyloseq") {
        stop("ps not a phyloseq object")
    }
    
    print(paste("Normalizing count using method", method))
    
    require(matrixStats)
    ## Extract count matrix from phyloseq object
    count_mat <- as(otu_table(ps), "matrix")
    
    ## extract normalizaed counts
    if (method == "RAW") {
        ## Raw counts no normalization applied
        norm_factors <- rep(1, ncol(count_mat))
    } else if (method == "UQ") {
        ## Upper quartile normalization
        count_mat[count_mat == 0] <- NA
        norm_factors <- colQuantiles(count_mat, p = 0.75 ,na.rm = TRUE)
    } else if (method == "CSS") {
        ## Cumulative sum scaling Paulson et al. 2013
        norm_factors <- metagenomeSeq::calcNormFactors(as(count_mat, "matrix"), p = 0.75)
        norm_factors <- norm_factors$normFactors
    } else if ( method == "TSS") {
        ## Total sum scaling
        norm_factors <- colSums(count_mat)
    } else if (method %in% c("RLE","TMM")) {
        ## To prevent log(0) issues
        count_mat = count_mat + 1
        ## EdgeR RLE and TMM normalization methods
        norm_factors <- edgeR::calcNormFactors(count_mat, method = method)
        norm_factors <- norm_factors * colSums(count_mat)
    } else {
        warning("Normalization method not defined")
    }
    
    ## Normalizing counts
    norm_mat <- sweep(count_mat, 2, norm_factors,'/')
    if ( sum(row.names(norm_mat) == taxa_names(ps)) == length(taxa_names(ps))) { 
        norm_tbl <- otu_table(norm_mat, taxa_are_rows = TRUE)
    } else if ( sum(colnames(norm_mat) == taxa_names(ps)) == length(taxa_names(ps))) {
        norm_tbl <- otu_table(norm_mat, taxa_are_rows = FALSE)
    } else {
        stop("taxa_names to not match row or column names")
    }
    
    
    ## Returning phyloseq object with normalized counts 
    ## Note: refseq data not included in returned phyloseq object 
    phyloseq(norm_tbl, sample_data(ps), tax_table(ps), phy_tree(ps))
    
}

apply_norm_methods_ps <- function(ps){
    ## Normalizing counts using different methods/ normlization factors
    norm_methods <- list(RAW = "RAW", ## No normalization
                         RLE = "RLE", ## EdgeR - relative log expression
                         TMM = "TMM", ## EdgeR - weighted trim mean of M-values
                         UQ = "UQ",   ## EdgeR - upperquartile
                         CSS = "CSS", ## metagenomeSeq - cumulative sum scaling, with p = 0.75
                         TSS = "TSS") ## Total sum scaling (proportions)
    
    map(norm_methods, normalize_counts,ps) 
}



## Removing no template controls
remove_ntc <- function(ps){
    non_ntc_samples <- sample_data(ps)$biosample_id != "NTC" 
    
    prune_samples(non_ntc_samples, ps)
}

## Removing samples with no reads
remove_no_read_samples <- function(ps) prune_samples(sample_sums(ps) > 0,  ps)

## List of phyoseq objects from different pipelines being evaluated
ps_list <- list(dada = dada_ps, 
                mothur = mothur_ps, 
                qiimeOpenRef = qiime_open_ref_ps,
                qiimeClosedRef = qiime_closed_ref_ps,
                qiimeDeNovo = qiime_de_novo_ps) 

## Removing NTC and samples with no counts
ps_no_ntc_list <- ps_list %>% map(remove_ntc) %>% map(remove_no_read_samples)

####################### Caching Normalized Count Data ##########################
## cache objects are nested list with 
##  - different pipeline phyloseq objects as the first level,
##  - and normalization methods as the second level normalization
## 
## methods split into two objects as only rareified counts are used to evaluate
## unweighted metrics

## Caching list of phyloseq objects with rarified counts  
ProjectTemplate::cache("rarified_ps_list",{map(ps_no_ntc_list, apply_rare_methods_ps)})

## Caching list of phyloseq objects with normalized counts
ProjectTemplate::cache("normalized_ps_list",{map(ps_no_ntc_list, apply_norm_methods_ps)})
