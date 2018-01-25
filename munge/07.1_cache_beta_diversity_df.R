## Generating Data Frame with weighted and unweighted diversity metrics 
diversity_rds <- list.file("munge/diversity_data", full.names = TRUE)

diversity_names <- basename(diversity_rds) 

diversity_df <- data_frame(div_info, rds_file = diversity_rds) %>% 
    separate(div_info, c("pipe","method","dist_method")) %>% 
    mutate(dist_method = paste0(dist_method, "_dist")) %>% 
    select(-div_info) 

############## Unweighted metrics - rarified data only #########################
## Unweighted metrics include unweighted unifrac and jaccard
ProjectTemplate::cache("unweighted_beta_df", {
    diversity_df %>%
        filter(dist_method %in% c("unifrac", "jaccard")) %>%
        mutate(dist_results = map(rds_file, readRDS)) %>% 
        select(-rds_file)
})

########### Weighted Metrics ###################################################
## Weighted metrics evaluated compared include weighted unifrac and bray-curtis
ProjectTemplate::cache("weighted_beta_df", {
    diversity_df %>%
        filter(dist_method %in% c("wunifrac", "bray")) %>%
        mutate(dist_results = map(rds_file, readRDS)) %>% 
        select(-rds_file)
})