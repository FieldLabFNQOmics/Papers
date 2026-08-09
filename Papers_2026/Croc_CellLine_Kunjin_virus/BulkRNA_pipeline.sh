
##### Example Read processing and Alignment #######

fastqc ./Rawreads/LV1KunV96hrsR1_R1.fastq.gz --outdir ./Rawreads/

./programs/bbmap/bbduk.sh -Xmx20g threads=10 ref=./TruSeq3-SE.fa in1=./Rawreads/LV1KunV96hrsR1_R1.fastq.gz out1=./Trimmed_reads//LV1KunV96hrsR1_bbduk_trim_R1.fastq.gz ktrim=r k=23 mink=11 hdist=1 tbo qtrim=r trimq=20 minlen=50

fastqc ./Trimmed_reads//LV1KunV96hrsR1_bbduk_trim_R1.fastq.gz --outdir ./Trimmed_reads//

./programs/STAR/bin/Linux_x86_64_static/STAR --outSAMattrRGline ID:LV1KunV96hrsR1 SM:LV1KunV96hrsR1_l1 PL:ILLUMINA --sjdbOverhang 99 --readFilesCommand zcat --runThreadN 28 --genomeDir /g/data/pq84/software/STAR/Croc_ens/100bp --outFileNamePrefix ./Trimmed_reads//LV1KunV96hrsR1_ --outSAMtype BAM SortedByCoordinate --quantMode GeneCounts --readFilesIn ./Trimmed_reads//LV1KunV96hrsR1_bbduk_trim_R1.fastq.gz

samtools index ./Trimmed_reads//LV1KunV96hrsR1_Aligned.sortedByCoord.out.bam

###### Rscript ########

library(SummarizedExperiment)
library(knitr)
library(dplyr)
library(ggforce)
library(scales)
library(reshape2)  # for melt
library(cowplot)   # for plot_grid
library(umap)
library(Rtsne)
library(pheatmap)
library(ggrepel)
library(stringr)
library(consensusDE)

file_list <- list.files(path=getwd(),pattern = ".bam$",full=TRUE)

sample_table <- data.frame("file"=basename(file_list),"group"=c("LV1KunV192hrs","LV1KunV192hrs","LV1KunV192hrs","LV1KunV48hrs","LV1KunV48hrs","LV1KunV48hrs","Control","Control","Control","LV1KunV8hrs","LV1KunV8hrs","LV1KunV8hrs","LV1KunV96hrs","LV1KunV96hrs","LV1KunV96hrs"))

bam_dir <- as.character(gsub(basename(file_list)[1], "", file_list[1]))

summarized_Croc <- buildSummarized(sample_table = sample_table, bam_dir = bam_dir, read_format = "single", output_log="/g/data/pq84/rnaseq/Sarker//bam/", gtf="/g/data/pq84/software/STAR/Croc_ens/Croc_ens.gtf")

saveRDS(summarized_Croc, file = "summarized_Croc.rds")

summarized_Croc <- readRDS("./bamData/summarized_Croc.rds")

sample_name_column <- colnames(summarized_Croc)

sample_name_column <- c("LV1KunV192hrsR1","LV1KunV192hrsR2","LV1KunV192hrsR3","LV1KunV48hrsR1","LV1KunV48hrsR2","LV1KunV48hrsR3","LV1KunV8hrsC1","LV1KunV8hrsC2","LV1KunV8hrsC3","LV1KunV8hrsR1","LV1KunV8hrsR2","LV1KunV8hrsR3","LV1KunV96hrsR1","LV1KunV96hrsR2","LV1KunV96hrsR3")

#group_name_column <-  c("LV1KunV192hrs","LV1KunV192hrs","LV1KunV192hrs","LV1KunV48hrs","LV1KunV48hrs","LV1KunV48hrs","LV1KunV8hrs","LV1KunV8hrs","LV1KunV8hrs","LV1KunV8hrs","LV1KunV8hrs","LV1KunV8hrs","LV1KunV96hrs","LV1KunV96hrs","LV1KunV96hrs")

colData(summarized_Croc)$file <-  c("LV1KunV192hrsR1","LV1KunV192hrsR2","LV1KunV192hrsR3","LV1KunV48hrsR1","LV1KunV48hrsR2","LV1KunV48hrsR3","LV1KunV8hrsC1","LV1KunV8hrsC2","LV1KunV8hrsC3","LV1KunV8hrsR1","LV1KunV8hrsR2","LV1KunV8hrsR3","LV1KunV96hrsR1","LV1KunV96hrsR2","LV1KunV96hrsR3")

#colData(summarized_Croc)$group <- c("LV1KunV192hrs","LV1KunV192hrs","LV1KunV192hrs","LV1KunV48hrs","LV1KunV48hrs","LV1KunV48hrs","LV1KunV8hrs","LV1KunV8hrs","LV1KunV8hrs","LV1KunV8hrs","LV1KunV8hrs","LV1KunV8hrs","LV1KunV96hrs","LV1KunV96hrs","LV1KunV96hrs")

colData(summarized_Croc)$pairs <- c("LV1KunV192hrs","LV1KunV192hrs","LV1KunV192hrs","LV1KunV48hrs","LV1KunV48hrs","LV1KunV48hrs","Control","Control","Control","LV1KunV8hrs","LV1KunV8hrs","LV1KunV8hrs","LV1KunV96hrs","LV1KunV96hrs","LV1KunV96hrs")

colnames(summarized_Croc) <- sample_name_column

summarized_Croc$file
?edgeR::DGEList()
summarized_Croc_filter <- buildSummarized(summarized = summarized_Croc,
                                        sample_table = sample_table,
                                        filter = TRUE,
                                        output_log = "./bamData/")

summarized_Croc_filter$group

saveRDS(summarized_Croc_filter, file = "summarized_Croc_filter.rds")

summarized_Croc_filter <- readRDS("summarized_Croc_filter.rds")

all_pairs_DE_filter_ruv <- multi_de_pairs(summarized = summarized_Croc_filter,
                                          paired = "unpaired",
                                          ruv_correct = TRUE,
                                         # plot_dir = "./bamData/multi_DE_outputs/plots/",
                                        #  output_voom = "./bamData/multi_DE_outputs/voom_results/",
                                         # output_edger = "./bamData/multi_DE_outputs/edger_results/",
                                        #  output_deseq = "./bamData/multi_DE_outputs/deseq_results/",
                                        #  output_combined = "./bamData/multi_DE_outputs/combined_results2/"
                                          
)

?multi_de_pairs()

all_pairs_DE_filter_ruv_noLabel <- multi_de_pairs(summarized = summarized_Croc_filter,
                                          paired = "unpaired",
                                          ruv_correct = TRUE,
                                          plot_dir = "./bamData/multi_DE_outputs/plot_noLabel/",
                                          label = FALSE
)

colData(summarized_Croc_filter)$pairs <- c("pair2","pair2","pair2","pair2","pair2","pair2","pair1","pair1","pair1","pair2","pair2","pair2","pair2","pair2","pair2")

all_pairs_DE_filter_ruv_noLabel_paired <- multi_de_pairs(summarized = summarized_Croc_filter,
                                                  paired = "paired",
                                                  ruv_correct = TRUE,
                                                  plot_dir = "./bamData/multi_DE_outputs/plot_noLabel_paired/",
                                                  label = FALSE
)

all_pairs_DE_filter_ruv$merged

uniq_comaprisons <- names(all_pairs_DE_filter_ruv$merged)



### p_intersect

setwd("./bamData//")


dir.create("p_intersect3")

setwd("./bamData/p_intersect3/")

pdf(file="Volcano_comparisons_p_intersect.pdf")

pair <- "LV1KunV192hrs-Control"


##########
library(biomaRt)

listEnsembl(version=113)

ensembl13 = useEnsembl(biomart="ensembl",version=113)

croc_ensembl <- useEnsembl(biomart="ensembl", dataset="cporosus_gene_ensembl")

# List available attributes (columns in the dataset)
attributes <- listAttributes(croc_ensembl)
print(head(attributes))


croc_genes <- getBM(
  attributes = c(
    "ensembl_gene_id",
    "external_gene_name",
    "description"
  ),
  mart = croc_ensembl
)

# View the first few rows
head(croc_genes)

write.csv(croc_genes, "crocodile_genes.csv", row.names = FALSE)
########

for (pair in uniq_comaprisons) {
  pairs <- str_split_fixed(pair,"-",2)
  x_label = paste("Enriched in ",pairs[,2], "<- log2(FC) -> Enriched in ", pairs[,1])
  results2_sub <- subset(all_pairs_DE_filter_ruv$merged[[pair]] )
  results2_sub$Color[results2_sub$p_intersect < 0.05] <- "p_intersect < 0.05"
  results2_sub$Color[results2_sub$p_intersect < 0.005] <- "p_intersect < 0.005"
  results2_sub$Color[results2_sub$p_intersect < 0.001] <- "p_intersect < 0.001"
  results2_sub$Color[results2_sub$p_intersect >= 0.05  ] <- "p_intersect  >= 0.05"
  results2_sub$Color[abs(results2_sub$LogFC) < 1] <- "Log2FoldChange < 1"
  results2_sub$Color <- factor(results2_sub$Color,
                               levels = c("p_intersect < 0.05",
                                          "p_intersect < 0.005", "p_intersect < 0.001","p_intersect  >= 0.05", "Log2FoldChange < 1"))
  # merge gene annotation
  results2_sub <- merge(
    results2_sub, 
    croc_genes, 
    by.x = "ID", 
    by.y = "ensembl_gene_id", 
    all.x = TRUE
  )
  # significant label subset
  label_data <- subset(
    results2_sub,
    p_intersect < 0.001 & abs(LogFC) > 1
  )
  
  # choose label column
  label_data$Label <- ifelse(
    is.na(label_data$external_gene_name) | label_data$external_gene_name=="",
    label_data$external_gene_name
  )
  
  plot <- ggplot(results2_sub, aes(x = LogFC, y = -log10(p_intersect),color = Color, label = ID)) +
    geom_vline(xintercept = c(1, -1), lty = "dashed") +
    geom_point() +
    labs(x = x_label,y = "Significance, -log10(P)",color = "Significance") +
    scale_color_manual(values = c("p_intersect < 0.001" = "dodgerblue","p_intersect < 0.005" = "lightblue", "p_intersect < 0.05" = "orange2", "p_intersect  >= 0.05" = "yellow", "Log2FoldChange < 2" = "gray"),guide = guide_legend(override.aes = list(size = 4))) +
   # geom_text_repel(data = subset(results2_sub, results2_sub$p_intersect < 0.001 & (results2_sub$LogFC > 1 | results2_sub$LogFC < -1) ),size = 3, point.padding = 0.15, color = "black",min.segment.length = .1, box.padding = .2,max.overlaps = 10) +
    geom_text_repel(data = label_data,  aes(label = Label),size = 3, point.padding = 0.15, color = "black",min.segment.length = .1, box.padding = .2,max.overlaps = 10) +
    coord_cartesian(ylim = c(0, 60)) + theme_bw(base_size = 6) +
    theme(legend.position = "bottom",legend.text = element_text(size = 6), panel.border = element_blank(), panel.grid.major = element_blank(),
          panel.grid.minor = element_blank(), axis.line = element_line(colour = "black"))
  print(plot)
  
  filename1 = paste(pairs[,1],"_vs_",pairs[,2],"_Enriched_",pairs[,1],".tsv",sep="")
  pdffilename = paste( pairs[,1],"_vs_",pairs[,2],"_Enriched_",pairs[,1],".pdf",sep="")
  ggsave(pdffilename, plot = plot, device = "pdf")
  write.table(results2_sub,file=filename1,sep="\t",row.names = FALSE)
}

dev.off()

#############################

setwd("./bamData/p_intersect3/")


#pair <- "LV1KunV192hrs-Control"


pdf(file="Volcano_comparisons_p_intersect_log5.pdf")

#max((-log10(results2_sub$p_intersect)))
#min((-log10(results2_sub$p_intersect)))

for (pair in uniq_comaprisons) {
  pairs <- str_split_fixed(pair,"-",2)
  x_label = paste("Enriched in ",pairs[,2], "<- log2(FC) -> Enriched in ", pairs[,1])
  results2_sub <- subset(all_pairs_DE_filter_ruv$merged[[pair]] )
  results2_sub$Color[results2_sub$LogFC > 1 & results2_sub$p_intersect < 0.05] <- "p_intersect < 0.05 & Log2FoldChange > 1"
  # results2_sub$Color[results2_sub$p_intersect < 0.005] <- "p_intersect < 0.005"
  #  results2_sub$Color[results2_sub$p_intersect < 0.001] <- "p_intersect < 0.001"
  results2_sub$Color[results2_sub$LogFC < -1 & results2_sub$p_intersect < 0.05] <- "p_intersect < 0.05 & Log2FoldChange < -1"
  # results2_sub$Color[results2_sub$p_intersect >= 0.05  ] <- "p_intersect  >= 0.05"
  results2_sub$Color[abs(results2_sub$LogFC) < 1 | results2_sub$p_intersect >= 0.05 ] <- "p_intersect  >= 0.05 & Log2FoldChange < 1"
  results2_sub$Color <- factor(results2_sub$Color,
                               levels = c("p_intersect < 0.05",
                                          "p_intersect  >= 0.05 & Log2FoldChange < 1",  "p_intersect < 0.05 & Log2FoldChange < -1", "p_intersect < 0.05 & Log2FoldChange > 1"))
  # merge gene annotation
  results2_sub <- merge(
    results2_sub, 
    croc_genes, 
    by.x = "ID", 
    by.y = "ensembl_gene_id", 
    all.x = TRUE
  )
  # significant label subset
  label_data <- subset(
    results2_sub,
    p_intersect < 0.001 & abs(LogFC) > 5
  )
  
  # choose label column
  label_data$Label <- label_data$external_gene_name
  list(label_data$Label)
#  label_data$Label <- ifelse(
 #   is.na(label_data$external_gene_name) | label_data$external_gene_name=="",
#    label_data$external_gene_name
 # )
  
  
  plot <- ggplot(results2_sub, aes(x = LogFC, y = -log10(p_intersect),color = Color, label = ID)) +
    geom_vline(xintercept = c(1, -1), lty = "dashed") +
    geom_hline(yintercept = c(-log10(0.05)), lty = "dashed") +
    geom_point(size = 2) +
    labs(x = x_label,y = "Significance, -log10(P)",color = "Significance : ") +
    # scale_color_manual(values = c("p_intersect < 0.001" = "dodgerblue","p_intersect < 0.005" = "lightblue", "p_intersect < 0.05" = "orange2", "p_intersect  >= 0.05" = "yellow", "Log2FoldChange < 2" = "gray"),guide = guide_legend(override.aes = list(size = 4))) +
    scale_color_manual(values = c("p_intersect < 0.05 & Log2FoldChange > 1" = "dodgerblue", "p_intersect < 0.05 & Log2FoldChange < -1" = "red", "p_intersect  >= 0.05 & Log2FoldChange < 1" = "darkgrey"),guide = guide_legend(override.aes = list(size = 3))) +
    #geom_text_repel(data = subset(results2_sub, results2_sub$p_intersect < 0.0001 & (results2_sub$LogFC > 3 | results2_sub$LogFC < -3) ),size = 2, point.padding = 0.15, color = "black",min.segment.length = .1, box.padding = .2,max.overlaps = 50) +
    
    geom_text_repel(data = label_data,  aes(label = Label),size = 3, point.padding = 0.15, color = "black",min.segment.length = .1, box.padding = .2,max.overlaps = 10) +
    
    coord_cartesian(ylim = c(0, max((-log10(results2_sub$p_intersect))))) + theme_bw(base_size = 12) + 
    theme(legend.position = "bottom",legend.text = element_text(size = 8), panel.border = element_blank(), panel.grid.major = element_blank(),
          panel.grid.minor = element_blank(), axis.line = element_line(colour = "black"),   axis.text = element_text(size = 8) )
 # print(plot)
  
 # filename1 = paste(pairs[,1],"_vs_",pairs[,2],"_log5.tsv",sep="")
  pdffilename = paste( pairs[,1],"_vs_",pairs[,2],"_log5.pdf",sep="")
  ggsave(pdffilename, plot = plot, device = "pdf", width = 10, height = 8)
#  write.table(results2_sub,file=filename1,sep="\t",row.names = FALSE)
}
dev.off()


######## GO/Over representation analyses ######

library(biomaRt)

library(clusterProfiler)


ensembl = useEnsembl(biomart="ensembl", dataset="hsapiens_gene_ensembl")

ensembl = useEnsembl(biomart="ensembl", dataset="hsapiens_gene_ensembl")


uniq_comaprisons <- names(all_pairs_DE_filter_ruv$merged)

pair <- "LV1KunV192hrs-Control"

res_DE <- all_pairs_DE_filter_ruv$merged[[pair]]

############################## LV1KunV192hrs-Control######## KEGG Plot ###

pair <- "LV1KunV192hrs-Control" 
pair <- "LV1KunV48hrs-Control" 
pair <- "LV1KunV8hrs-Control" 
pair <- "LV1KunV96hrs-Control" 
pair <- "LV1KunV192hrs-LV1KunV48hrs" 
pair <- "LV1KunV192hrs-LV1KunV8hrs"  
pair <- "LV1KunV192hrs-LV1KunV96hrs" 
pair <- "LV1KunV48hrs-LV1KunV8hrs"   
pair <- "LV1KunV48hrs-LV1KunV96hrs"  
pair <- "LV1KunV8hrs-LV1KunV96hrs"

for (pair in uniq_comaprisons) {
  pairs <- str_split_fixed(pair,"-",2)
  res_DE <- all_pairs_DE_filter_ruv$merged[[pair]]
  
#LPS Up : LogFC > 1
croc_ensembl <- useEnsembl(biomart="ensembl", dataset="cporosus_gene_ensembl")

#topTable_pVal <- dplyr::filter(res_DE, p_intersect < 0.05 & (LogFC >1 | LogFC < -1   )) %>% dplyr::select(ID, LogFC)

topTable_pVal <- dplyr::filter(res_DE, p_intersect < 0.05 & (LogFC >1 | LogFC < -1    )) %>% dplyr::select(ID, LogFC)

#print(length(topTable_pVal))
#table(topTable_pVal)

ensembl_gene_id <- topTable_pVal$ID

print(length(ensembl_gene_id))

croc_entrez_ids <- getBM(
  attributes = c("ensembl_gene_id", "entrezgene_id"),
  filters = "ensembl_gene_id",
  values = unique(ensembl_gene_id),
  mart = croc_ensembl
)

# Remove missing Entrez IDs
croc_entrez_ids <- na.omit(croc_entrez_ids)

kk <- enrichKEGG(
  gene = unique(croc_entrez_ids$entrezgene_id),
  organism = "cpoo",  
  pvalueCutoff = 0.05
)

# Plot KEGG Pathways
kk_plot <- dotplot(kk, showCategory = 20) + ggtitle("Top 20 Over-Represented KEGG Pathways")

#ggsave("KEGG_Pathways_control_vs_192hrs.pdf")


pdffilename = paste("KEGG_Pathways_", pairs[,1],"_vs_",pairs[,2],".pdf",sep="")
ggsave(pdffilename, plot = kk_plot, device = "pdf", width = 10, height = 8)


}


########## test 3: Top GO Terms by Category" #######

for (pair in uniq_comaprisons) {
  str_split_fixed(pair, "-", 2)
  res_DE <- all_pairs_DE_filter_ruv$merged[[pair]]
  
  croc_ensembl <- useEnsembl(biomart="ensembl", dataset="cporosus_gene_ensembl")
  
  # Filter significant genes and add regulation direction
  topTable_pVal <- dplyr::filter(res_DE, p_intersect < 0.05 & (LogFC > 1 | LogFC < -1)) %>% 
    dplyr::select(ID, LogFC) %>%
    dplyr::mutate(regulation = ifelse(LogFC > 1, "Upregulated", "Downregulated"))
  
  ensembl_gene_id <- topTable_pVal$ID
  
  # Get GO data with categories
  croc_genes <- getBM(
    attributes = c("ensembl_gene_id", "external_gene_name", "go_id", "name_1006", "namespace_1003"),
    mart = croc_ensembl
  )
  
  # Clean data and remove empty entries
  croc_genes_clean <- croc_genes[!is.na(croc_genes$go_id) & 
                                   croc_genes$go_id != "" & 
                                   !is.na(croc_genes$namespace_1003) &
                                   croc_genes$namespace_1003 != "", ]
  
  # Merge with regulation information
  croc_genes_with_reg <- merge(croc_genes_clean, 
                               topTable_pVal[, c("ID", "regulation")], 
                               by.x = "ensembl_gene_id", 
                               by.y = "ID", 
                               all.x = TRUE)
  
  # Filter to only include genes in our significant list
  croc_genes_with_reg <- croc_genes_with_reg[!is.na(croc_genes_with_reg$regulation), ]
  
  # Count GO terms by category and regulation
  go_summary <- croc_genes_with_reg %>%
    group_by(namespace_1003, name_1006, regulation) %>%
    summarise(gene_count = n_distinct(ensembl_gene_id), .groups = "drop") %>%
    group_by(namespace_1003, regulation) %>%
    top_n(10, gene_count) %>%
    arrange(namespace_1003, regulation, desc(gene_count)) %>%
    ungroup()
  
  # Create mirrored data: negative values for downregulated, positive for upregulated
  go_summary_mirrored <- go_summary %>%
    mutate(
      gene_count_mirrored = ifelse(regulation == "Downregulated", -gene_count, gene_count),
      abs_count = abs(gene_count_mirrored)
    )
  
  # Get top terms by absolute count for better ordering
  top_terms <- go_summary_mirrored %>%
    group_by(namespace_1003, name_1006) %>%
    summarise(max_count = max(abs_count), .groups = "drop") %>%
    group_by(namespace_1003) %>%
    top_n(15, max_count) %>%
    arrange(namespace_1003, max_count)
  
  # Filter to only top terms
  go_plot_data <- go_summary_mirrored %>%
    semi_join(top_terms, by = c("namespace_1003", "name_1006"))
  
  # Create the mirrored plot
  p <- ggplot(go_plot_data, aes(x = reorder(name_1006, abs_count), 
                                y = gene_count_mirrored, 
                                fill = regulation)) +
    geom_bar(stat = "identity") +
    coord_flip() +
    facet_wrap(~namespace_1003, scales = "free_y", ncol = 1) +
    scale_fill_manual(values = c("Upregulated" = "#d73027", "Downregulated" = "#4575b4")) +
    scale_y_continuous(
      labels = abs,  # Show absolute values on axis
      breaks = function(x) pretty(c(-max(abs(x)), max(abs(x))))
    ) +
    theme_minimal() +
    theme(
      axis.text.y = element_text(size = 8),
      axis.text.x = element_text(size = 8),
      strip.text = element_text(size = 10, face = "bold"),
      legend.position = "bottom",
      panel.grid.minor = element_blank(),
      strip.background = element_blank()
    ) +
    labs(
      title = paste("GO Terms by Regulation Direction -", pair),
      #subtitle = "Left: Downregulated | Right: Upregulated",
      x = "GO Term",
      y = "Gene Count",
      fill = "Regulation"
    ) +
    geom_vline(xintercept = 0, linetype = "dashed", alpha = 0.5)
  
  print(p)
  
  # Optional: Save the plot
   ggsave(paste0("top10_GO_Terms_", pair, ".pdf"), p, width = 12, height = 10, dpi = 300)
}
