## Calculating Beta Diversity Metrics

## Safely calculate distance metrics  
safe_dist <- safely(phyloseq::distance)

## Unweighted metrics - rarified data only
## Generate Tidy data frame with phyloseq objects as a column
rare_df <- transpose(rarified_ps_list) %>% 
    as_data_frame() %>% 
    add_column(pipe = names(rarified_ps_list)) %>% 
    gather(key = "method",value = "ps_obj", -pipe)

## Unweighted Metrics
unweighted_beta_df <- rare_df %>% 
    mutate(unifrac_dist = map(ps_obj, phyloseq::distance, method = "unifrac")) %>% 
    mutate(jaccard_dist = map(ps_obj, phyloseq::distance, method = "jaccard")) %>% 
    select(-ps_obj)


ProjectTemplate::cache("unweighted_beta_df")

## Weighted metrics normized and rarified data
## Generate Tidy data frame with phyloseq objects as a column
norm_df <- transpose(normalized_ps_list) %>% 
    as_data_frame() %>% 
    add_column(pipe = names(normalized_ps_list)) %>% 
    gather(key = "method",value = "ps_obj", -pipe) %>% 
    bind_rows(rare_df)

## Weighted Metrics
weighted_beta_df <- rare_df %>% 
    mutate(wunifrac_dist = map(ps_obj, safe_dist, method = "wunifrac")) %>%
    mutate(bray_dist = map(ps_obj, safe_dist, method = "bray")) %>% 
    select(-ps_obj)

ProjectTemplate::cache("weighted_beta_df")
w_dist <- map(weighted_beta_df$wunifrac_dist,1)
w_dist
devtools::install_github("joey711/phyloseq")
