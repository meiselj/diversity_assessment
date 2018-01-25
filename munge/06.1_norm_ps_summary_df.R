## Rarified Count data summary 
ProjectTemplate::cache("norm_summary_df", {
    tibble(rds_file = list.files("munge/norm_data", full.names = TRUE)) %>% 
        mutate(norm_ps = map(rds_file, readRDS),
               n_taxa = map_dbl(norm_ps, ntaxa),
               n_samples = map_dbl(norm_ps, nsamples)) %>% 
        select(-norm_ps)
})
