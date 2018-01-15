## Normalize phyloseq objects

## Rarefaction
## 5000 in most analysis
## Weiss et al 2017? compared 2000, 5000, and 10000
## McMurdie and Holmes 2014 used the 15th percentile library size
## Fix number after rerunning pipeline with modified filtering methods. 
apply_rare_methods_ps <- function(ps){
    rare_levels <- list(rare_2K = 200, rare_5K = 500, rare_10K = 1000)
    
    rare_ps <- rare_levels %>% 
        map(~rarefy_even_depth(ps, rngseed = 531, sample.size = .))
    
    q15 <- floor(quantile(sample_sums(ps), 0.15, na.rm = TRUE))
    rare_ps$rare_15P <- rarefy_even_depth(ps, rngseed = 531, sample.size = q15)
    
    rare_ps
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

ps_list <- list(dada = dada_ps, mothur = mothur_ps, qiimeOpenRef = dada_ps)
ps_no_ntc_list <- ps_list %>% map(remove_ntc) %>% map(remove_no_read_samples)

ProjectTemplate::cache("rarified_ps_list",{map(ps_no_ntc_list, apply_rare_methods_ps)})

ProjectTemplate::cache("normalized_ps_list",{map(ps_no_ntc_list, apply_norm_methods_ps)})
