## Calculate alpha diversity from phyloseq objects
## Not sure whether to use vegan for iNext

## ## Function for generating a data.frame with alpha diversity metrics from an MRexperiment object
## Formatted for mgtst sample data
calc_alpha <- function(mrexp){
    count_mat <- mrexp@assayData$counts
    sam_dat <- pData(mrexp) %>% rownames_to_column(var = "samID") %>% 
        mutate(pcr_16S_plate = as.character(pcr_16S_plate))
    shan_div <- vegan::diversity(count_mat,index = "shannon",2)
    simp_div <- vegan::diversity(count_mat, index = "simpson",2)
    spec_div <- vegan::specnumber(count_mat, MARGIN = 2)  
    data_frame(samID = colnames(count_mat), 
               Shannon = shan_div, Simpson = simp_div, Observed = spec_div) %>% 
        left_join(sam_dat) %>% gather("div_index","value",Shannon:Observed)
}

## Calculate Breakaway diversity metrics