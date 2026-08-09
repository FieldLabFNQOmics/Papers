###########  

library(ConsensusMetaDA)

biom_file <- "./Fecal_feature-table-tax.biom"

sample_table_file <- "./Fecal_sample_file.txt"

fecal_bcg <- build_OTU_counts(biom = biom_file, sample_table = sample_table_file,  
                                     abundance_threshold = 10, 
                                     prevalence_threshold = 0, 
                                     rarity_threshold = 0, 
                                     variance_threshold = 0  )

biom_file <- "./Lung_feature-table-tax.biom"

sample_table_file <- "./Lung_sample_file.txt"

lung_bcg <- build_OTU_counts(biom = biom_file, sample_table = sample_table_file,  
                                     abundance_threshold = 0, 
                                     prevalence_threshold = 0, 
                                     rarity_threshold = 0, 
                                     variance_threshold = 0  )


plots_fecal_bcg <- ConsensusMetaDA::OTUs_plots(fecal_bcg)

plots_fecal_bcg2 <- ConsensusMetaDA::OTUs_plot(fecal_bcg)



