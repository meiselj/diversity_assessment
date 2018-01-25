setwd("~/diversity_assessment/")
library(ProjectTemplate)
library(reshape2)
library(plyr)
load.project()

# Extract on only the unmixed pre samples (t_fctr==20) or unmixed post samples (t_fctr==0)
map_sub<-subset(mgtstMetadata, 
                mgtstMetadata$t_fctr %in% c(0, 20) & mgtstMetadata$biosample_id != "NTC")
# Fix factoring of some columns
map_sub<-droplevels(map_sub)

# Look at all combinations of samples
sample_comparisons<-as.data.frame(t(combn(as.character(map_sub$sample_id), 2)))
colnames(sample_comparisons)<-c("sample_a", "sample_b")

# Determine what characteristics are shared by each sample pair
sample_comparisons<-merge(sample_comparisons, map_sub, by.x="sample_a", by.y="sample_id")
colnames(sample_comparisons)[3:ncol(sample_comparisons)]<-
  paste0( "sample_a_",colnames(sample_comparisons)[3:ncol(sample_comparisons)])
sample_comparisons<-merge(sample_comparisons, map_sub, by.x="sample_b", by.y="sample_id")
colnames(sample_comparisons)[9:ncol(sample_comparisons)]<-
  paste0("sample_b_",colnames(sample_comparisons)[9:ncol(sample_comparisons)])

# Same individual
sample_comparisons$same_individual<-
  (sample_comparisons$sample_a_biosample_id==sample_comparisons$sample_b_biosample_id)
sample_comparisons$sample_a_biosample_id<-NULL
sample_comparisons$sample_b_biosample_id<-NULL

# Same timepoint
sample_comparisons$same_timepoint<-
  (sample_comparisons$sample_a_t_fctr==sample_comparisons$sample_b_t_fctr)
sample_comparisons$sample_a_t_fctr<-NULL
sample_comparisons$sample_b_t_fctr<-NULL

# Same PCR product 
sample_comparisons$same_pcr_product<-
  ((sample_comparisons$sample_a_pcr_16S_plate==sample_comparisons$sample_b_pcr_16S_plate) & 
     (sample_comparisons$sample_a_pos == sample_comparisons$sample_b_pos))
sample_comparisons$sample_a_pcr_16S_plate<-NULL
sample_comparisons$sample_b_pcr_16S_plate<-NULL
sample_comparisons$sample_a_pos<-NULL
sample_comparisons$sample_b_pos<-NULL

# Same sequencing lab
sample_comparisons$same_seq_lab<-
  (sample_comparisons$sample_a_seq_lab==sample_comparisons$sample_b_seq_lab)
sample_comparisons$sample_a_seq_lab<-NULL
sample_comparisons$sample_b_seq_lab<-NULL

# Same sequencing run
sample_comparisons$same_seq_run<-
  (sample_comparisons$sample_a_seq_run==sample_comparisons$sample_b_seq_run)
sample_comparisons$sample_a_seq_run<-NULL
sample_comparisons$sample_b_seq_run<-NULL

# Classify variation (biological or technical) and add variation labels (between timepoint variation, etc.)
## biological variation
# between timepoint variation
between_timepoint<-subset(sample_comparisons, sample_comparisons$same_timepoint==FALSE &
                                     sample_comparisons$same_seq_lab==TRUE &
                                     sample_comparisons$same_seq_run==TRUE)
between_timepoint$variation<-"biological"
between_timepoint$variation_label<-"between_timepoint"

# between individual variation
between_individual<-subset(sample_comparisons, sample_comparisons$same_individual==FALSE & 
                            sample_comparisons$same_timepoint==TRUE &
                            sample_comparisons$same_seq_lab==TRUE &
                            sample_comparisons$same_seq_run==TRUE)
between_individual$variation<-"biological"
between_individual$variation_label<-"between_individuals_within_timepoint"

# within individual timepoint variation
within_individual_timepoint<-subset(sample_comparisons, sample_comparisons$same_individual==TRUE & 
                            sample_comparisons$same_timepoint==FALSE &
                            sample_comparisons$same_seq_lab==TRUE &
                            sample_comparisons$same_seq_run==TRUE)
within_individual_timepoint$variation<-"biological"
within_individual_timepoint$variation_label<-"within_individuals_between_timepoint"

variation<-rbind(between_timepoint, between_individual)
variation<-rbind(variation, within_individual_timepoint)
rm(between_timepoint, between_individual, within_individual_timepoint)

## technical variation
# between seq labs
between_seq_labs<-subset(sample_comparisons, sample_comparisons$same_individual==TRUE & 
                                      sample_comparisons$same_timepoint==TRUE &
                                      sample_comparisons$same_pcr_product==TRUE &
                                      sample_comparisons$same_seq_lab==FALSE)
between_seq_labs$variation<-"technical"
between_seq_labs$variation_label<-"between_seq_labs"

# within seq labs run variation
within_seq_lab_runs<-subset(sample_comparisons, sample_comparisons$same_individual==TRUE & 
                           sample_comparisons$same_timepoint==TRUE &
                           sample_comparisons$same_seq_lab==TRUE &
                             sample_comparisons$same_pcr_product==TRUE &
                             sample_comparisons$same_seq_run==FALSE)
within_seq_lab_runs$variation<-"technical"
within_seq_lab_runs$variation_label<-"within_seq_lab_runs"

# within seq labs PCR product variation
within_seq_lab_pcr<-subset(sample_comparisons, sample_comparisons$same_individual==TRUE & 
                              sample_comparisons$same_timepoint==TRUE &
                              sample_comparisons$same_seq_lab==TRUE &
                              sample_comparisons$same_seq_run==TRUE &
                              sample_comparisons$same_pcr_product==FALSE)
within_seq_lab_pcr$variation<-"technical"
within_seq_lab_pcr$variation_label<-"within_seq_lab_pcr"

variation<-rbind(variation, between_seq_labs)
variation<-rbind(variation, within_seq_lab_runs)
variation<-rbind(variation, within_seq_lab_pcr)
rm(between_seq_labs, within_seq_lab_runs, within_seq_lab_pcr)

# Order variation types
variation$variation_label<-factor(variation$variation_label, 
                                  levels=c("between_individuals_within_timepoint", 
                                           "within_individuals_between_timepoint", 
                                           "between_timepoint", "within_seq_lab_pcr", 
                                           "within_seq_lab_runs","between_seq_labs"))
# Plot what comparisons are being made
variation_tmp<-unique(variation[,c(3:7,9)])
# Clean up and merge duplicates in same variation label
variation_tmp[1,1]<-"TRUE OR FALSE"
variation_tmp<-variation_tmp[-2,]

variation_tmp[4,5]<-"TRUE OR FALSE"
variation_tmp<-variation_tmp[-5,]

variation_comparisons<-melt(variation_tmp, id.vars = "variation_label")
rm(variation_tmp)

variation_comparisons$value<-factor(variation_comparisons$value, levels = c("TRUE", "TRUE OR FALSE", "FALSE"), ordered = TRUE)

## Generate table for metrics
# data is unweighted_beta_df$unifrac_dist
compute_diversity_comparisons<-function(data, dist_mat, map, metric, variation_tests)
{
  beta_div_comp<- lapply(1:length(data$pipe), function(i){
    if(metric=="bray_curtis" | metric=="weighted_unifrac"){
      beta_div<-as.matrix(dist_mat[[i]][[1]])
    } else if (metric == "unweighted_unifrac" | metric == "jaccard"){
      beta_div<-as.matrix(dist_mat[[i]])
    }
    
    beta_div_m<-as.data.frame(beta_div)[which(row.names(beta_div) %in% as.character(map$sample_id)), 
                                        which(colnames(beta_div) %in% as.character(map$sample_id))]
    beta_div_m$sample<-row.names(beta_div_m)
    beta_div_m2<-melt(beta_div_m, id.vars=c("sample"))
    colnames(beta_div_m2)<-c("sample_a", "sample_b", "value")
    beta_div_m2$pipe<-data$pipe[i]
    beta_div_m2$normalization<-data$method[i]
    beta_div_m2$metric<-metric
    return(beta_div_m2)
  })
  output<-do.call("rbind", beta_div_comp)
  output2<-merge(variation_tests, output, by=c("sample_a", "sample_b"))
  output2$normalization<-factor(output2$normalization, levels = c("rare_2K", "rare_5K", "rare_10K", "rare_15P"), ordered = T)
  return(output2)
}

unweighted_unifrac_comp<-compute_diversity_comparisons(unweighted_beta_df, unweighted_beta_df$unifrac_dist, map_sub, "unweighted_unifrac", variation)
jaccard_comp<-compute_diversity_comparisons(unweighted_beta_df, unweighted_beta_df$jaccard_dist, map_sub, "jaccard", variation)
weighted_unifrac_comp<-compute_diversity_comparisons(weighted_beta_df, weighted_beta_df$wunifrac_dist, map_sub, "weighted_unifrac", variation)
bray_curtis_comp<-compute_diversity_comparisons(weighted_beta_df, weighted_beta_df$bray_dist, map_sub, "bray_curtis", variation)

# Only need to calculate for one metric since it should be same number of overall comparisons for all
comparison_numbers<-ddply(weighted_unifrac_comp, c("variation", "variation_label", "pipe", "normalization"), summarize, number_of_comparisons=length(variation_label))

comparisons<-rbind(bray_curtis_comp, jaccard_comp)
comparisons<-rbind(comparisons, unweighted_unifrac_comp)
comparisons<-rbind(comparisons, weighted_unifrac_comp)
comparisons_summary<-ddply(comparisons, c("variation", "variation_label", "pipe", "normalization", "metric"), summarize, mean_value=mean(value))
rm(comparisons)
comparisons_summary$pipeline_labels<-paste0(comparisons_summary$pipe, "_", comparisons_summary$normalization)
comparisons_summary$pipeline_labels<-factor(comparisons_summary$pipeline_labels, 
                                            levels = c("dada_rare_2K", "dada_rare_5K", "dada_rare_10K", "dada_rare_15P",
                                                       "mothur_rare_2K", "mothur_rare_5K", "mothur_rare_10K", "mothur_rare_15P",
                                                       "qiimeOpenRef_rare_2K", "qiimeOpenRef_rare_5K", "qiimeOpenRef_rare_10K", "qiimeOpenRef_rare_15P"),
                                            ordered = T)


pdf("./reports/biological_v_technical_variation_diversity_plots.pdf", height=8, width=12)

ggplot(variation_comparisons, aes(y=variation_label, x=variable, fill=value))+
  geom_tile()+theme_bw()+theme(axis.text.x=element_text(angle = -45, hjust = 0))+
  scale_fill_manual(values = c("black", "gray", "white"))+xlab("")+ylab("Type of Variation")

ggplot(comparison_numbers, aes(x=variation_label, y=number_of_comparisons, fill=variation))+
  geom_bar(stat="identity")+theme_bw()+ theme(axis.text.x=element_text(angle = -45, hjust = 0))+
  facet_grid(pipe~normalization)+xlab("Type of Variation")+
  scale_fill_manual(values = c("#d53e4f", "#3288bd"))

ggplot(unweighted_unifrac_comp, aes(x=variation_label, y=value, fill=variation))+
  geom_boxplot()+theme_bw()+ theme(axis.text.x=element_text(angle = -45, hjust = 0))+
  facet_grid(pipe~normalization)+xlab("Type of Variation") +ylab(unique(unweighted_unifrac_comp$metric))+
  scale_fill_manual(values = c("#d53e4f", "#3288bd"))


ggplot(jaccard_comp, aes(x=variation_label, y=value, fill=variation))+
  geom_boxplot()+theme_bw()+ theme(axis.text.x=element_text(angle = -45, hjust = 0))+
  facet_grid(pipe~normalization)+xlab("Type of Variation") +ylab(unique(jaccard_comp$metric))+
  scale_fill_manual(values = c("#d53e4f", "#3288bd"))


ggplot(weighted_unifrac_comp, aes(x=variation_label, y=value, fill=variation))+
  geom_boxplot()+theme_bw()+ theme(axis.text.x=element_text(angle = -45, hjust = 0))+
  facet_grid(pipe~normalization)+xlab("Type of Variation") +ylab(unique(weighted_unifrac_comp$metric))+
  scale_fill_manual(values = c("#d53e4f", "#3288bd"))


ggplot(bray_curtis_comp, aes(x=variation_label, y=value, fill=variation))+
  geom_boxplot()+theme_bw()+ theme(axis.text.x=element_text(angle = -45, hjust = 0))+
  facet_grid(pipe~normalization)+xlab("Type of Variation") +ylab(unique(bray_curtis_comp$metric))+
  scale_fill_manual(values = c("#d53e4f", "#3288bd"))


ggplot(comparisons_summary, aes(x=normalization, y=variation_label, fill=mean_value))+
  geom_tile()+theme_bw()+ theme(axis.text.x=element_text(angle = -45, hjust = 0))+
  facet_grid(metric~pipe, scales="free")+ scale_fill_distiller(palette = "YlOrRd", direction = 1)

dev.off()
