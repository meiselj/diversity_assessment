# ## Normalization
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
    phyloseq(norm_tbl, sample_data(ps), tax_table(ps), phy_tree(ps))

}
