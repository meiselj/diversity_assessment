## Calculating Beta Diversity Metrics for Normalized count data

## Requirements for parallelizing distance calculations
# require(doParallel) 
# registerDoParallel(cores = 7)

## Safely calculate distance metrics Instead of stopping when there is an error
## when calculating the distance metric provides NULL results with error message
safe_dist <- safely(phyloseq::distance)
safe_unifrac <- safely(phyloseq::UniFrac) 

############## Unweighted metrics - rarified data only #########################
## Unweighted metrics include unweighted unifrac and jaccard 

## Generate Tidy data frame with phyloseq objects as a column
rare_df <- transpose(rarified_ps_list) %>%
    as_data_frame() %>%
    add_column(pipe = names(rarified_ps_list)) %>%
    gather(key = "method",value = "ps_obj", -pipe)

## Calculate Unweighted Metrics ------------------------------------------------
unweighted_beta_df <- rare_df %>%
    mutate(unifrac_dist = map(ps_obj, safe_unifrac, 
                              ## UniFrac parameters
                              weighted = FALSE, parallel = TRUE)) %>%
    mutate(jaccard_dist = map(ps_obj, safe_dist, method = "jaccard")) %>%
    select(-ps_obj)

### Caching unweighted metrics results -----------------------------------------
ProjectTemplate::cache("unweighted_beta_df")

########### Weighted Metrics ###################################################
## Weighted metrics evaluated compared include weighted unifrac and bray-curtis

## Generate Tidy data frame with phyloseq objects as a column - including both
## normalized and rarified counts
norm_df <- transpose(normalized_ps_list) %>%
    as_data_frame() %>%
    add_column(pipe = names(normalized_ps_list)) %>%
    gather(key = "method",value = "ps_obj", -pipe) %>%
    bind_rows(rare_df)

##  Calculating Weighted Metrics -----------------------------------------------
weighted_beta_df <- norm_df %>%
    mutate(wunifrac_dist = map(ps_obj,  safe_unifrac, 
                               ## UniFrac parameters
                               weighted = TRUE, parallel = TRUE)) %>%
    mutate(bray_dist = map(ps_obj, safe_dist, method = "bray")) %>%
    select(-ps_obj)

### Caching weighted metrics results -------------------------------------------
ProjectTemplate::cache("weighted_beta_df")
