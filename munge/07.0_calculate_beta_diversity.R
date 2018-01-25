## Calculating Beta Diversity Metrics for Normalized count data
calc_beta_div <- function(rds_file, div_metric){
    norm_ps <- readRDS(rds_file)

    div_file <- rds_file %>%
        str_replace("norm_data","diversity_data") %>%
        str_replace(".rds", paste0("_", div_metric, ".rds"))

    ## Check for diversity rds file
    if (file.exists(div_file)) {
        print("Skipping calculation diversity estimate already calculated")
        return("")
    }

    print(paste("- Calculating", div_metric))

    ## Used to determine if dataset too large for parallel diversity estimates
    n_taxa <- ntaxa(norm_ps)

    ## calc beta
    if (div_metric == "wunifrac") {
        if (n_taxa > 50000) {
            div_results <- safe_unifrac(norm_ps, weighted = TRUE)
        } else {
            div_results <- safe_unifrac(norm_ps, weighted = TRUE, parallel = TRUE)
        }
    } else if (div_metric == "unifrac") {
        if (n_taxa > 50000) {
            div_results <- safe_unifrac(norm_ps, weighted = FALSE)
        } else {
            div_results <- safe_unifrac(norm_ps, weighted = FALSE, parallel = TRUE)
        }
    } else if (div_metric == "jaccard") {
        div_results <- safe_dist(norm_ps, method = "jaccard")
    } else if (div_metric == "bray") {
        div_results <- safe_dist(norm_ps, method = "bray")
    } else {
        stop("Diversity metric not one of the expected values: wunifrac, unifrac, jaccard, or bray")
    }

    ## Save beta rds
    saveRDS(div_results, div_file)
}

## Calculate diversity metrics
norm_to_div <- function(rds_file) {
    print(paste("Processing:", rds_file))

    ## File check
    if (!file.exists(rds_file)) stop(paste(rds_file, "does not exist"))

    ## Weighted Metrics - calculate for all normalization methods
    calc_beta_div(rds_file, "wunifrac")
    calc_beta_div(rds_file, "bray")

    ## Unweighted Metrics - only calculated for rarified and raw data
    if (grepl("rare|RAW", x = rds_file)) {
        calc_beta_div(rds_file, "unifrac")
        calc_beta_div(rds_file, "jaccard")
    }

}

## Processing normalized count data
rds_files <- list.files("data/norm_data", full.names = TRUE) %>% 
    walk(norm_to_div)
