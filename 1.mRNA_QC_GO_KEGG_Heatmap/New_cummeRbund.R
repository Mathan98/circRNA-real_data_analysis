#R

# load cummeRbund
library(cummeRbund)
library(dplyr)
library(GGally)
library(tidyr)
library(magrittr)
library(tidyverse)
library(stringr)


## H1N1_human_tophat
setwd('C:/Users/matha/Documents/Cufflinks/H1N1_human_allvsall_cuffdiff_tophat')

## H1N1 human
cuff <- readCufflinks(gtfFile="C:/Users/matha/Documents/Cufflinks/human_gencode_vch38.gtf",
                      genome="C:/Users/matha/Documents/Cufflinks/Homo_sapiens.fa")

##############################################################################################

## H1N1_human_STAR
setwd('C:/Users/matha/Documents/Cufflinks/H1N1_human_allvsall_cuffdiff_STAR')

## H1N1 human
cuff <- readCufflinks(gtfFile="C:/Users/matha/Documents/Cufflinks/human_gencode_vch38.gtf",
                      genome="C:/Users/matha/Documents/Cufflinks/Homo_sapiens.fa")


##################################################################################################################

## H1N1 human_virus_TopHat
setwd('C:/Users/matha/Documents/Cufflinks/H1N1_human_virus_allvsall_cuffdiff_tophat')

## loading the Cuffdiff results

## H1N1 human_virus
cuff <- readCufflinks(gtfFile="C:/Users/matha/Documents/Cufflinks/human_virus.gtf",
                      genome="C:/Users/matha/Documents/Cufflinks/human_virus.fa")



## H1N1 human_virus_STAR
# setwd('C:/Users/matha/Documents/Cufflinks/H1N1_human_virus_allvsall_cuffdiff_STAR')

## H1N1_human_virus maxFrags
#setwd('C:/Users/matha/Documents/H1N1_human_virus_allvsall_cuffdiff_STAR_maxFrags')


## loading the Cuffdiff results

## H1N1 human_virus
 # cuff <- readCufflinks(gtfFile="C:/Users/matha/Documents/H1N1_human_virus_allvsall_cuffdiff_STAR/human_virus.gtf",
 #                       genome="C:/Users/matha/Documents/H1N1_human_virus_allvsall_cuffdiff_STAR/human_virus.fa")

#################################################################################################################
 
## H3N2 human
 setwd('C:/Users/matha/Documents/Cufflinks/H3N2_human_allvsall_cuffdiff')

## H3N2_human_tophat
cuff <- readCufflinks(gtfFile="C:/Users/matha/Documents/Cufflinks/H3N2_human_allvsall_cuffdiff/human_gencode_vch38.gtf",
                      genome="C:/Users/matha/Documents/Cufflinks/H3N2_human_allvsall_cuffdiff/Homo_sapiens.fa")


##################################################################################################################


## H5N1_1 human
 setwd('C:/Users/matha/Documents/Cufflinks/H5N1_1_human_allvsall_cuffdiff')

## H5N1_1human_tophat
 cuff <- readCufflinks(gtfFile="C:/Users/matha/Documents/Cufflinks/H5N1_1_human_allvsall_cuffdiff/human_gencode_vch38.gtf",
                       genome="C:/Users/matha/Documents/Cufflinks/H5N1_1_human_allvsall_cuffdiff/Homo_sapiens.fa")


#############################################################################################

 
 ## H5N1_2 human
 setwd('C:/Users/matha/Documents/Cufflinks/H5N1_2_human_allvsall_cuffdiff')
 
 ## H5N1_2human_tophat
 cuff <- readCufflinks(gtfFile="C:/Users/matha/Documents/Cufflinks/H5N1_1_human_allvsall_cuffdiff/human_gencode_vch38.gtf",
                       genome="C:/Users/matha/Documents/Cufflinks/H5N1_1_human_allvsall_cuffdiff/Homo_sapiens.fa")
 
 
############################################################################################
 ################################################################################################################################
 
 ## GTF dataframe
 
 gtf <- rtracklayer::import("C:/Users/matha/Documents/Simulated_reads/human_gencode_vch38.gtf") %>% as.data.frame()
 gtf_pc <- gtf %>% filter(gene_type == "protein_coding" | gene_type == 'transcribed_unprocessed_pseudogene' |
                            gene_type == "transcribed_processed_pseudogene" | gene_type == "processed_transcript" | 
                            gene_type == "transcribed_unitary_pseudogene" | gene_type == "IG_V_pseudogene" |
                            gene_type == "IG_C_gene" | gene_type == 'TR_C_gene' |
                            gene_type == "TR_V_gene" | gene_type == "TR_D_gene" | 
                            gene_type == "IG_J_pseudogene" | gene_type == "IG_pseudogene" |
                            gene_type == "unprocessed_pseudogene" | gene_type == 'processed_pseudogene' |
                            gene_type == "unitary_pseudogene" | gene_type == "polymorphic_pseudogene" | 
                            gene_type == "IG_V_gene" | gene_type == "IG_J_gene" |
                            gene_type == "TR_J_gene" | gene_type == 'TR_V_pseudogene' |
                            gene_type == "IG_C_pseudogene" | gene_type == "TR_J_pseudogene" | 
                            gene_type == "IG_D_gene")
 gtf_geneNames <- unique(gtf_pc$gene_name) %>% as.data.frame()
 
 ################################################################################################################################
 
## Get diff genes
gene_diff_data <- diffData(genes(cuff))
# isoform_diff_data <- diffData(isoforms(cuff))

## Add names genes
# isoform.features <- annotation(isoforms(cuff))
gene.features<-annotation(genes(cuff))
gene_diff_data$gene_name <- gene.features$gene_short_name

# ##Add names transcripts
# isoform.features<-annotation(isoforms(cuff))
# isoform_diff_data$isoform_name <- isoform.features$nearest_ref_id
# isoform_diff_data$gene_name <- isoform.features$gene_short_name

### Filter for human or virus genes
## For human
# ViralGeneIds<-read.delim2("C:/Users/matha/Documents/Cufflinks/H1N1_human_virus_allvsall_cuffdiff_STAR/virus_geneIDs.txt", header = F) %>% as.matrix()
# gene_diff_data <- gene_diff_data %>% filter(!(gene_name %in% ViralGeneIds))

## For virus
# ViralGeneIds<-read.delim2("C:/Users/matha/Documents/Cufflinks/H1N1_human_virus_allvsall_cuffdiff_STAR/virus_geneIDs.txt", header = F) %>% as.matrix()
# gene_diff_data <- gene_diff_data %>% filter(gene_name %in% ViralGeneIds)



## Extract comparisons with their respective time points
gene_diff = list()
for ( i in c(3, 6, 12, 18)) {
  gene_diff[[i]] <- gene_diff_data %>% 
    filter(sample_1 == paste("C", i, "H", sep="")) %>% 
    filter(sample_2 == paste("T", i, "H", sep=""))
}

gene_diff_data <- do.call(rbind, gene_diff)
gene_diff_data <- gene_diff_data %>% 
  mutate(gene_name = strsplit(gene_name, ",")) %>%
  unnest_longer(gene_name) 

gene_diff_data <- gene_diff_data[(gene_diff_data$gene_name %in% gtf_geneNames$.),] 
gene_diff_data <- gene_diff_data %>% distinct(gene_id, gene_name, sample_1, .keep_all = T)
gene_diff_data <- gene_diff_data %>% 
  group_by_at(vars(-gene_name)) %>%
  summarise(.groups = "keep", 
            gene_name = str_c(gene_name, collapse=",")) %>%
  as.data.frame()


# isoform_diff = list()
# for ( i in c(3, 6, 12, 18)) {
#   isoform_diff[[i]] <- isoform_diff_data %>% 
#     filter(sample_1 == paste("C", i, "H", sep="")) %>% 
#     filter(sample_2 == paste("T", i, "H", sep=""))
# }
# 
# isoform_diff_data <- do.call(rbind, isoform_diff)


################################################################################################################################

## remove lowly expressed
# rep FPKM
rep_FPKM <- repFpkmMatrix(genes(cuff))
rep_FPKM$gene_name <- gene.features$gene_short_name

rep_FPKM_long <- rep_FPKM %>% pivot_longer(!gene_name, names_to = "Samples", values_to = "count")
p<-ggplot(rep_FPKM_long, aes(x=log10(count), fill=Samples)) +
  geom_density(alpha=0.2)
p


rep_FPKM_1 <- rep_FPKM %>%
  is_greater_than(0.4) %>%
  rowSums() %>%
  is_greater_than(2) %>%
  as.data.frame() %>%
  filter(. == "TRUE")

rep_FPKM <- rep_FPKM[rownames(rep_FPKM) %in% rownames(rep_FPKM_1), ]
rep_FPKM$gene_id <- rownames(rep_FPKM)

rep_FPKM_unnest <- rep_FPKM %>% 
  mutate(gene_name = strsplit(gene_name, ",") ) %>%
  unnest_longer(gene_name)

rep_FPKM_unnest <- rep_FPKM_unnest[(rep_FPKM_unnest$gene_name %in% gtf_geneNames$.),]

rep_FPKM <- rep_FPKM_unnest %>% as.data.frame() %>%
  group_by_at(vars(-gene_name)) %>%
  summarise(.groups = "keep", gene_name = str_c(gene_name, collapse=",")) %>%
  as.data.frame()

rownames(rep_FPKM) <- rep_FPKM$gene_id
rep_FPKM <- rep_FPKM %>% select(-gene_id)



## Find DE
gene_diff_data_DE <- gene_diff_data %>% filter(gene_id %in% rownames(rep_FPKM)) %>% filter(abs(log2_fold_change) >= 2, q_value <= 0.05)

## Timepoint
Timepoint_3hrs <- gene_diff_data_DE %>% filter(sample_1 == "C3H")
Timepoint_6hrs <- gene_diff_data_DE %>% filter(sample_1 == "C6H")
Timepoint_12hrs <- gene_diff_data_DE %>% filter(sample_1 == "C12H")
Timepoint_18hrs <- gene_diff_data_DE %>% filter(sample_1 == "C18H")

## Summarise gene table
mRNA_table <- gene_diff_data %>% group_by(sample_1) %>%
  summarise(n = n())

## Summarise DE table
mRNA_DE_table <- gene_diff_data_DE %>% group_by(sample_1) %>%
  summarise(n = n())

################################################################################################################################
## Distribution plot

rep_FPKM_long <- rep_FPKM %>% pivot_longer(!gene_name, names_to = "Samples", values_to = "count")
rep_FPKM_long$Samples <- factor(rep_FPKM_long$Samples, levels = c('C3H_0', 'C3H_1', 'T3H_0', 'T3H_1', 'C6H_0', 'C6H_1', 'T6H_0', 'T6H_1', 'C12H_0', 'C12H_1', 'T12H_0', 'T12H_1', 'C18H_0', 'C18H_1', 'T18H_0', 'T18H_1'))

p<-ggplot(rep_FPKM_long, aes(x=log10(count), fill=Samples)) +
  geom_density(alpha=0.5) +
  theme(axis.text.x = element_text(size = 15),
        axis.title.x = element_text(size = 15, face = "bold"),
        axis.text.y = element_text(size = 15),
        axis.title.y = element_text(size = 15, face = "bold"))
p
# Add mean lines
p+geom_vline(data=rep_FPKM_long, aes(xintercept=grp.mean, color=Samples),
             linetype="dashed")

################################################################################################################################
## Boxplot linear mRNAs


b <- csBoxplot(genes(cuff), replicates =T)
b


sample_fpkm <- b[["data"]]
################################################################################################################################
## Select FPKM data

library(ggpubr)

# colnames(rep_FPKM) <- c('C3H_0', 'C3H_1', 
#                         'T3H_0', 'T3H_1', 
#                         'C6H_0', 'C6H_1', 
#                         'T6H_0', 'T6H_1', 
#                         'C12H_0', 'C12H_1', 
#                         'T12H_0', 'T12H_1', 
#                         'C18H_0', 'C18H_1', 
#                         'T18H_0', 'T18H_1')

# rep_FPKM$isoforms <- isoform.features$gene_short_name
# rep_FPKM_longer <- pivot_longer(rep_FPKM[,-c(17)], cols = everything())
# rep_FPKM_longer$name <- factor(rep_FPKM_longer$name, levels = c('C3H_0', 'C3H_1', 'T3H_0', 'T3H_1', 'C6H_0', 'C6H_1', 'T6H_0', 'T6H_1', 'C12H_0', 'C12H_1', 'T12H_0', 'T12H_1', 'C18H_0', 'C18H_1', 'T18H_0', 'T18H_1'))


rep_fpkm_3 <- pivot_longer(rep_FPKM[,-c(17)], cols = everything(), names_to = "condition", values_to = "fpkm")
rep_fpkm_3$condition <- factor(rep_fpkm_3$condition, levels = c('C3H_0', 'C3H_1', 'T3H_0', 'T3H_1', 'C6H_0', 'C6H_1', 'T6H_0', 'T6H_1', 'C12H_0', 'C12H_1', 'T12H_0', 'T12H_1', 'C18H_0', 'C18H_1', 'T18H_0', 'T18H_1'))



## Boxplot
Boxplot_1 = ggplot(rep_fpkm_3 , aes(x=condition, y= log10(as.numeric(fpkm)), fill=condition)) +  
  labs(x = "Samples", y = "log10(FPKM)", title = "FPKM Distribution of Human Genes", fill = "Sample_names") + 
  geom_boxplot() + theme(axis.text.x = element_text(size = 20, angle=90),
                         axis.title.x = element_blank(),
                         axis.text.y = element_text(size = 20),
                         axis.title = element_text(size = 20),
                         plot.title = element_text(size = 25, hjust = 0.5), 
                         legend.position = "none") 
Boxplot_1


# boxplot(rep_FPKM_2[,-17] %>% +1 %>% log10(), xlab="", ylab="Log2 counts per million",las=2)

################################################################################################################################

# logFPKM_data<- gene_diff_data %>% select(sample_1, sample_2, value_1, value_2) %>% as_tibble()
# 
# ## melt data to factorise
# melted_FPKM_1 <- melt(logFPKM_data, id = c("sample_1", "sample_2"))
# melted_FPKM_2 <- melt(logFPKM_data, id = c("value_1", "value_2"))
# melted_FPKM <- cbind(melted_FPKM_2$value, melted_FPKM_1$value) %>% as_tibble()
# melted_FPKM$V1 <- factor(melted_FPKM$V1, levels = c('C3H', 'T3H', 'C6H', 'T6H', 'C12H', 'T12H', 'C18H', 'T18H'))
# 
# ## Boxplot
# Boxplot_1 = ggplot(melted_FPKM , aes(x=V1, y= log10(as.numeric(V2)), fill=V1)) +  
#   labs(x = "Samples", y = "log10(FPKM)", title = "FPKM Distribution of Human Genes", fill = "Sample_names") + 
#   geom_boxplot() + theme(axis.text = element_text(size = 20), axis.title = element_text(size = 20),
#                          plot.title = element_text(size = 25, hjust = 0.5), 
#                          legend.title = element_text(size = 15),
#                          legend.text = element_text(size = 13)) 
# Boxplot_1
# 
# ## Violin Plot
# Violin_1 = ggplot(melted_FPKM, aes(x=V1, y= log10(as.numeric(V2) + 1), fill=V1)) +  
#   labs(x = "Timepoints", y = "log10(FPKM)", title = "Human FPKM Distribution") + geom_violin(trim=FALSE) +
#   theme(axis.text = element_text(size = 20), axis.title = element_text(size = 20),
#         plot.title = element_text(size = 25, hjust = 0.5), 
#         legend.title = element_text(size = 15),
#         legend.text = element_text(size = 13)) 
# Violin_1

############################################################################################################################

## Scatterplot

ggplot(Timepoint_3hrs, aes(x=log10(as.numeric(value_1 + 1)), y= log10(as.numeric(value_2+ 1)))) + geom_point() + geom_smooth(method=lm) + 
  labs(x = "C3H (Log10 FPKM + 1)", y = "T3H (Log10 FPKM + 1)", title = "Human Control 3hrs vs Treated 3hrs")

ggplot(Timepoint_6hrs, aes(x=log10(as.numeric(value_1 + 1)), y= log10(as.numeric(value_2+ 1)))) + geom_point() + geom_smooth(method=lm) + 
  labs(x = "C6H (Log10 FPKM + 1)", y = "T6H (Log10 FPKM + 1)", title = "Human Control 6hrs vs Treated 6hrs")

ggplot(Timepoint_12hrs, aes(x=log10(as.numeric(value_1 + 1)), y= log10(as.numeric(value_2+ 1)))) + geom_point() + geom_smooth(method=lm) + 
  labs(x = "C12H (Log10 FPKM + 1)", y = "T12H (Log10 FPKM + 1)", title = "Human Control 12hrs vs Treated 12hrs")

ggplot(Timepoint_18hrs, aes(x=log10(as.numeric(value_1 + 1)), y= log10(as.numeric(value_2+ 1)))) + geom_point() + geom_smooth(method=lm) + 
  labs(x = "C18H (Log10 FPKM + 1)", y = "T18H (Log10 FPKM + 1)", title = "Human Control 18hrs vs Treated 18hrs")

library(ggpubr)

## Scatterplot matrix
gene_diff_data$sample_1 <- factor(gene_diff_data$sample_1, levels = c('C3H', 'C6H', 'C12H', 'C18H'))
gene_diff_data$sample_2 <- factor(gene_diff_data$sample_2, levels = c('T3H', 'T6H', 'T12H', 'T18H'))

ggplot(gene_diff_data, aes(x=log10(as.numeric(value_1 + 1)), y=log10(as.numeric(value_2+ 1)))) +
  geom_point() + geom_smooth(method=lm) + stat_cor(method = "spearman", label.y.npc="top", label.x.npc = "left", size=3) +
  facet_grid(sample_2~sample_1, labeller = label_both, scale = "free") +
  ggtitle("Log10 (FPKM + 1) of Ctrl vs Treated") +
  scale_color_brewer(palette = "Dark2") +
  xlab("Control") + 
  ylab("Treated") +
  theme(text = element_text(size = 15),
        axis.text = element_text(size = 15), axis.title = element_text(size = 20),
        plot.title = element_text(size = 20)) 


#################################################################################################################################

## Plot Heatmap
# Load library 
library(pheatmap)
library(gplots)
library(dendextend)
library(colorspace)
library(ggplot2)
library(readxl)
library(stringr)

myGeneIds <- rownames(rep_FPKM)
myGenes<-getGenes(cuff,myGeneIds)

d <- csHeatmap(myGenes,cluster='both')
d

########################################################################################################

## make names bold
make_bold_names <- function(mat, rc_fun, rc_names) {
  bold_names <- rc_fun(mat)
  ids <- rc_names %>% match(rc_fun(mat))
  ids %>%
    walk(
      function(i)
        bold_names[i] <<-
        bquote(bold(.(rc_fun(mat)[i]))) %>%
        as.expression()
    )
  bold_names
}

########################################################################################################
## FPKM DE

# Log transform the FPKM
FPKM <- rep_FPKM %>%
  mutate(log10(across(where(is.numeric), ~ .x + 1))) %>% 
  as.data.frame()

## DE genes
FPKM_de <- FPKM[(rownames(FPKM) %in% gene_diff_data_DE$gene_id),]
rownames(FPKM_de) <- FPKM_de$gene_name
FPKM_de <- FPKM_de %>% select(-gene_name)


FPKM_logs<- array(0, dim=dim(FPKM_de))

for (i in 1:ncol(FPKM_de)) {FPKM_logs[,i] <- FPKM_de[,i]}
FPKM_logs[, 1:16] <- sapply(FPKM_logs[, 1:16], as.numeric)

colnames(FPKM_logs) <-colnames(FPKM_de)
# rownames(FPKM_logs) <- rep_FPKM$gene_name

# Arrange col order
col.order <- c("C3H_0", "C3H_1", "C6H_0", "C6H_1", "C12H_0", "C12H_1", "C18H_0", "C18H_1", 
               "T3H_0", "T3H_1", "T6H_0", "T6H_1", "T12H_0", "T12H_1", "T18H_0", "T18H_1")
FPKM_logs <- FPKM_logs[, col.order]

#  Define colours
hmcols <- rev(redgreen(2750));

# plot heatmap
pheatmap(FPKM_logs,
                       cluster_rows = TRUE,
                       cluster_cols = F, clustering_distance_rows = "euclidean",
                       clustering_distance_cols = "euclidean", clustering_method = "complete",
                       cutree_rows = 1, cutree_cols = 1, fontsize = 18, face="bold",
                       fontsize_row = 8, color=hmcols,scale="row", show_rownames = F,
                       labels_col = make_bold_names(FPKM_logs, colnames, col.order))

########################################################################################################
## FPKM All

# Log transform the FPKM
FPKM <- rep_FPKM %>%
  mutate(log10(across(where(is.numeric), ~ .x + 1))) %>% 
  as.data.frame()

rownames(FPKM) <- FPKM$gene_name
FPKM <- FPKM %>% select(-gene_name)


FPKM_logs<- array(0, dim=dim(FPKM))

for (i in 1:ncol(FPKM)) {FPKM_logs[,i] <- FPKM[,i]}
FPKM_logs[, 1:16] <- sapply(FPKM_logs[, 1:16], as.numeric)

colnames(FPKM_logs) <-colnames(FPKM)
# rownames(FPKM_logs) <- rep_FPKM$gene_name

# Arrange col order
col.order <- c("C3H_0", "C3H_1", "C6H_0", "C6H_1", "C12H_0", "C12H_1", "C18H_0", "C18H_1", 
               "T3H_0", "T3H_1", "T6H_0", "T6H_1", "T12H_0", "T12H_1", "T18H_0", "T18H_1")
FPKM_logs <- FPKM_logs[, col.order]

#  Define colours
hmcols <- rev(redgreen(2750));

# plot heatmap
pheatmap(FPKM_logs,
         cluster_rows = TRUE,
         cluster_cols = FALSE, clustering_distance_rows = "euclidean",
         clustering_distance_cols = "euclidean", clustering_method = "complete",
         cutree_rows = 1, cutree_cols = 1, fontsize = 18, face="bold",
         fontsize_row = 8, color=hmcols,scale="row", show_rownames = F,
         labels_col = make_bold_names(FPKM_logs, colnames, col.order))

##############################################################################################################################

### Correlativity of replicates

library(ggpubr)

FPKM_logs_cor <- as.data.frame(FPKM_logs)

## Control
# 3hrs
Ctrl_3hrs_rep <- ggplot(data = as.data.frame(FPKM_logs_cor), aes(x = as.numeric(C3H_0), y = as.numeric(C3H_1))) +
  geom_point() +
  geom_smooth(method = "lm", se = FALSE) +
  stat_cor(method = "pearson", label.y = max(as.numeric(FPKM_logs_cor$C3H_1)), label.x = 0, size = 6) +
  labs(title = "3hrs Ctrl Replicates") +
  theme(axis.title.x = element_blank(),
        axis.title.y = element_blank(),
        axis.text = element_text(size = 20),
        plot.title = element_text(size = 20))

# 6hrs
Ctrl_6hrs_rep <- ggplot(data = as.data.frame(FPKM_logs_cor), aes(x = as.numeric(C6H_0), y = as.numeric(C6H_1))) +
  geom_point() +
  geom_smooth(method = "lm", se = FALSE) +
  stat_cor(method = "pearson", label.y = max(as.numeric(FPKM_logs_cor$C6H_1)), label.x = 0, size = 6) +
  labs(title = "6hrs Ctrl Replicates") +
  theme(axis.title.x = element_blank(),
        axis.title.y = element_blank(),
        axis.text = element_text(size = 20),
        plot.title = element_text(size = 20))

# 12hrs
Ctrl_12hrs_rep <- ggplot(data = as.data.frame(FPKM_logs_cor), aes(x = as.numeric(C12H_0), y = as.numeric(C12H_1))) +
  geom_point() +
  geom_smooth(method = "lm", se = FALSE) +
  stat_cor(method = "pearson", label.y = max(as.numeric(FPKM_logs_cor$C12H_1)), label.x = 0, size = 6) +
  labs(title = "12hrs Ctrl Replicates") +
  theme(axis.title.x = element_blank(),
        axis.title.y = element_blank(),
        axis.text = element_text(size = 20),
        plot.title = element_text(size = 20))

# 18hrs
Ctrl_18hrs_rep <- ggplot(data = as.data.frame(FPKM_logs_cor), aes(x = as.numeric(C18H_0), y = as.numeric(C18H_1))) +
  geom_point() +
  geom_smooth(method = "lm", se = FALSE) +
  stat_cor(method = "pearson", label.y = max(as.numeric(FPKM_logs_cor$C18H_1)), label.x = 0, size=6) +
  labs(title = "18hrs Ctrl Replicates") +
  theme(axis.title.x = element_blank(),
        axis.title.y = element_blank(),
        axis.text = element_text(size = 20),
        plot.title = element_text(size = 20))


## Treated

# 3hrs
Treated_3hrs_rep <- ggplot(data = as.data.frame(FPKM_logs_cor), aes(x = as.numeric(T3H_0), y = as.numeric(T3H_1))) +
  geom_point() +
  geom_smooth(method = "lm", se = FALSE) +
  stat_cor(method = "pearson", label.y = max(as.numeric(FPKM_logs_cor$T3H_1)), label.x = 0, size=6) +
  labs(title = "3hrs IAV-infected Replicates") +
  theme(axis.title.x = element_blank(),
        axis.title.y = element_blank(),
        axis.text = element_text(size = 20),
        plot.title = element_text(size = 20))

# 6hrs
Treated_6hrs_rep <- ggplot(data = as.data.frame(FPKM_logs_cor), aes(x = as.numeric(T6H_0), y = as.numeric(T6H_1))) +
  geom_point() +
  geom_smooth(method = "lm", se = FALSE) +
  stat_cor(method = "pearson", label.y = max(as.numeric(FPKM_logs_cor$T6H_1)), label.x = 0, size=6) +
  labs(title = "6hrs IAV-infected Replicates") +
  theme(axis.title.x = element_blank(),
        axis.title.y = element_blank(),
        axis.text = element_text(size = 20),
        plot.title = element_text(size = 20))

# 12hrs
Treated_12hrs_rep <- ggplot(data = as.data.frame(FPKM_logs_cor), aes(x = as.numeric(T12H_0), y = as.numeric(T12H_1))) +
  geom_point() +
  geom_smooth(method = "lm", se = FALSE) +
  stat_cor(method = "pearson", label.y = max(as.numeric(FPKM_logs_cor$T12H_1)), label.x = 0, size=6) +
  labs(title = "12hrs IAV-infected Replicates") +
  theme(axis.title.x = element_blank(),
        axis.title.y = element_blank(),
        axis.text = element_text(size = 20),
        plot.title = element_text(size = 20))

# 18hrs
Treated_18hrs_rep <- ggplot(data = as.data.frame(FPKM_logs_cor), aes(x = as.numeric(T18H_0), y = as.numeric(T18H_1))) +
  geom_point() +
  geom_smooth(method = "lm", se = FALSE) +
  stat_cor(method = "pearson", label.y = max(as.numeric(FPKM_logs_cor$T18H_1)), label.x = 0, size=6) +
  labs(title = "18hrs IAV-infected Replicates") +
  theme(axis.title.x = element_blank(),
        axis.title.y = element_blank(),
        axis.text = element_text(size = 20),
        plot.title = element_text(size = 20))



ctrl_rep_cor <- (Ctrl_3hrs_rep + Ctrl_6hrs_rep)/(Ctrl_12hrs_rep + Ctrl_18hrs_rep)
treated_rep_cor <- (Treated_3hrs_rep + Treated_6hrs_rep)/(Treated_12hrs_rep + Treated_18hrs_rep)



### Treated vs Treated

## 3vs6
Treated_3vs6hrs_rep <- ggplot(data = as.data.frame(FPKM_logs_cor), aes(x = ((as.numeric(T3H_0) + as.numeric(T3H_1))/2), 
                                                                        y = ((as.numeric(T6H_0) + as.numeric(T6H_1))/2))) +
  geom_point() +
  geom_smooth(method = "lm", se = FALSE) +
  stat_cor(method = "pearson", label.y = max(as.numeric(FPKM_logs_cor$T6H_1)), label.x = 1.2) +
  labs(title = "3hrs vs 6hrs Treated Replicates") +
  labs(x ="Treated 3hrs", y = "Treated 6hrs")

## 3vs12
Treated_3vs12hrs_rep <- ggplot(data = as.data.frame(FPKM_logs_cor), aes(x = ((as.numeric(T3H_0) + as.numeric(T3H_1))/2), 
                                                                        y = ((as.numeric(T12H_0) + as.numeric(T12H_1))/2))) +
  geom_point() +
  geom_smooth(method = "lm", se = FALSE) +
  stat_cor(method = "pearson", label.y = max(as.numeric(FPKM_logs_cor$T12H_1)), label.x = 1.2) +
  labs(title = "3hrs vs 12hrs Treated Replicates") +
  labs(x ="Treated 3hrs", y = "Treated 12hrs")


## 3vs18
Treated_3vs18hrs_rep <- ggplot(data = as.data.frame(FPKM_logs_cor), aes(x = ((as.numeric(T3H_0) + as.numeric(T3H_1))/2), 
                                                                        y = ((as.numeric(T18H_0) + as.numeric(T18H_1))/2))) +
  geom_point() +
  geom_smooth(method = "lm", se = FALSE) +
  stat_cor(method = "pearson", label.y = max(as.numeric(FPKM_logs_cor$T18H_1)), label.x = 1.2) +
  labs(title = "3hrs vs 18hrs Treated Replicates") +
  labs(x ="Treated 3hrs", y = "Treated 18hrs")


treated_using3hrs_cor <- (Treated_3vs6hrs_rep + Treated_3vs12hrs_rep)/(Treated_3vs18hrs_rep)

##############################################################################################################################
## Volcano plot


gene_diff_data$diffexpressed <- "Not DE"
gene_diff_data$diffexpressed[gene_diff_data$log2_fold_change > 2 & gene_diff_data$q_value < 0.05] <- "DE Upregulated"
gene_diff_data$diffexpressed[gene_diff_data$log2_fold_change < -2 & gene_diff_data$q_value < 0.05] <- "DE Downregulated"




## Substituting the min and max volumes
VC_gene_diff_data <- gene_diff_data
for (i in c(3, 6, 12, 18)) {
  
  ## Substituting INF for max number of that timepoint
  VC_gene_diff_data[(VC_gene_diff_data$sample_1 %in% paste0("C", i, "H")) & (VC_gene_diff_data$log2_fold_change == "Inf"),][7] <- max(VC_gene_diff_data$log2_fold_change[is.finite(VC_gene_diff_data$log2_fold_change) & (VC_gene_diff_data$sample_1 %in% paste0("C", i, "H"))])

  ## Substituting -INF for min number of that timepoint
  VC_gene_diff_data[(VC_gene_diff_data$sample_1 %in% paste0("C", i, "H")) & (VC_gene_diff_data$log2_fold_change == "-Inf"),][7] <- min(VC_gene_diff_data$log2_fold_change[is.finite(VC_gene_diff_data$log2_fold_change) & (VC_gene_diff_data$sample_1 %in% paste0("C", i, "H"))])
  
}


# gene_diff_data$sample_1 <- factor(gene_diff_data$sample_1, levels = c('C3H', 'C6H', 'C12H', 'C18H'))
# gene_diff_data$sample_2 <- factor(gene_diff_data$sample_2, levels = c('T3H', 'T6H', 'T12H', 'T18H'))


# Change timepoint in %in%
VP <- list()
for (i in c(3, 6, 12, 18)) {
  VP[[i]] <- ggplot(VC_gene_diff_data[(VC_gene_diff_data$sample_1 %in% paste0("C", i, "H")),], aes(x=log2_fold_change, y=-log10(q_value), col=diffexpressed), size = 2) + 
                        geom_point() + scale_color_manual(values = c( "green3", "red2", "grey22", "grey22")) + 
                        geom_hline(yintercept = -log10(0.05),
                                   linetype = "dashed") + 
                        geom_vline(xintercept = c(2, -2),
                                   linetype = "dashed",
                                   colour = c("red2", "green3")) +
                        labs(x = "logFC", y = "-log10 p value" ,colour = "", title = paste0("Human ", i, "hrs Linear Volcano Plot"), 
                              fill = "") +
                        scale_alpha_manual(values=c(2,2,2)) +
                        xlim(-10,10) +
                        ylim(0,3.5) +
                        theme(axis.text = element_text(size = 20), 
                              axis.title.x = element_text(size = 15,  hjust = 0.42),
                              axis.title.y = element_text(size = 15),
                              plot.title = element_text(size = 19, face = "bold", hjust = 0.5), 
                              legend.text = element_text(size=15) ) 
  
}

library(patchwork)
ggp<- (VP[[3]] + VP[[6]])/(VP[[12]] + VP[[18]]) & theme(legend.position = "bottom")  
ggp + plot_layout(guides = "collect")


#  # Add all in one plot
# ggplot(gene_diff_data, aes(x=log2_fold_change, y=-log10(p_value), col=diffexpressed), size = 2) + 
#   facet_grid(sample_2~sample_1, labeller = label_both, scale = "free") +
#   geom_point() + scale_color_manual(values = c( "grey22", "green3", "red2", "red2")) + 
#   labs(x = "logFC", y = "-log10 p value" ,colour = "Expression", title = "Human 18hrs Linear Volcano Plot")
###########################################################################

d <- csDendro(genes(cuff),replicates=T)
d

## dendogram

dist <- dist(t(log10(rep_FPKM[, -17]+1)) , diag=TRUE)

# Hierarchical Clustering with hclust
hc <- hclust(dist)

# Plot the result
plot(hc)

###########################################################################
## MDS plot

dat<-repCountMatrix(genes(cuff)) %>% +1 %>% log10()
plotMDS(dat,dim=c(1,2))

## MDS plot
dat<-repCountMatrix(genes(cuff)) %>% +1 %>% log10()

samples <- c("red", "red", "orange", "orange",
             "darkgreen", "darkgreen", "chartreuse4", "chartreuse4",
             "darkblue", "darkblue","cornflowerblue", "cornflowerblue",
             "purple", "purple", "darkgoldenrod", "darkgoldenrod")

plotMDS(dat, col=samples)




# dat<-log10(dat+1)
# d<-JSdist(makeprobs(dat))
# fit <- cmdscale(d,eig=TRUE, k=2)
# res<-data.frame(names=rownames(fit$points),M1=fit$points[,1],M2=fit$points[,2])
# p <- ggplot(res)
# p <- p + geom_point(aes(x=M1,y=M2,color=names), size = 2) + 
#   # geom_text(aes(x=M1,y=M2,label=names,color=names), size = 5) +
#   geom_text_repel(aes(x=M1,y=M2,label=names, color = names), size = 5) + labs(color = "Sample names") +
#   theme_bw(base_size = 20)
# p


##############################################################################################################################
## PCA plot

# genes.PCA.rep<-PCAplot(isoforms(cuff),"PC1","PC2",replicates=T)
# genes.PCA.rep


data <- repFpkmMatrix(isoforms(cuff)) %>% +1 %>% log10()
pca <- prcomp(data, scale=T)

# trace(fviz, edit=TRUE)

fviz_pca_biplot(pca,
                label="var", labelsize = 6, arrowsize = 0.7,
                repel = T,
                col.var = factor(c("C3H_0", "C3H_1", "C6H_0", "C6H_1", "C12H_0", "C12H_1", "C18H_0", "C18H_1", 
                                   "T3H_0", "T3H_1", "T6H_0", "T6H_1", "T12H_0", "T12H_1", "T18H_0", "T18H_1")),
                xlab = "PC 1 (70%)",
                ylab = "PC 2 (25%)") +
  theme(axis.text = element_text(size = 15),
        text = element_text(size = 15))

fviz_eig(pca, addlabels = TRUE)

##############################################################################################################################

## Density  plot

# d <- csDensity(isoforms(cuff), replicates = T) 
# d

##############################################################################################################################

##Bar expression plot
library(tidyverse)
library(ggforce)

## Human
myGeneIds<-read.table("C:/Users/matha/Documents/custom_geneIDs.txt") %>% as.matrix()


## Virus
#myGeneIds<-read.delim2("C:/Users/matha/Documents/H1N1_human_virus_allvsall_cuffdiff_STAR/virus_geneIDs.txt") %>% as.matrix()


rep_FPKM <- repFpkmMatrix(genes(cuff))
gene.features<-annotation(genes(cuff))
rep_FPKM$gene_name <- gene.features$gene_short_name
rep_FPKM_1 <- rep_FPKM %>% transform(gene_name = str_split(gene_name, ",")) %>% unnest(gene_name)

rep_FPKM_custom <- rep_FPKM_1[(rep_FPKM_1$gene_name %in% myGeneIds),]
colnames(rep_FPKM_custom) <- c("C_3H_0", "C_3H_1", "T_3H_0", "T_3H_1", "C_6H_0", "C_6H_1", "T_6H_0", "T_6H_1", 
                               "C_12H_0", "C_12H_1", "T_12H_0", "T_12H_1", "C_18H_0", "C_18H_1", "T_18H_0", "T_18H_1",
                               "gene_name")

# Arrange the data and get sufficient information
X <- rep_FPKM_custom %>% 
  pivot_longer(cols = colnames(rep_FPKM_custom)[-17],
               names_to = c("Condition", "Time", "Rep"),
               names_sep = "_",
               values_to = "FPKM",
               values_transform = list(FPKM = as.numeric)) %>% 
  mutate(FPKM = FPKM) %>%
  mutate(across(Time, ~ as.numeric(gsub("H", "", .x)))) %>%
  group_by(gene_name, Condition, Time) %>%
  summarise(n = n(),
            mean = mean(FPKM),
            sd = sd(FPKM),
            se = sd / sqrt(n)) %>% 
  mutate(Time = str_trim(paste0(Time, "hrs", sep=""), side = "both"))  %>%
  mutate(Gene_con = paste0(Time, "_", Condition)) %>%
  ungroup()


## Plot barplots [ Change page number (up to 4 if ncol and nrow = 2) ]

X_1<-  X %>% mutate(Gene_con = fct_relevel(Gene_con, 
                            "3hrs_C", "3hrs_T", "6hrs_C", 
                            "6hrs_T", "12hrs_C", "12hrs_T", 
                            "18hrs_C", "18hrs_T"))


  
ggplot(X_1, aes(x = Gene_con, y = mean)) + 
  geom_col(aes(fill = Gene_con)) +
  geom_errorbar(aes(ymin = mean - se, ymax = mean + se), width=0.2) +
  xlab("Condition vs Treatment across Timepoints") +
  ylab("Expression Level (FPKM)") +
  theme(text = element_text(size=20),
        axis.text.y = element_text(size = 15),
        axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1, size = 15), 
        axis.title = element_text(size = 15, face = "bold"),
        plot.title = element_text(size = 10), legend.position="none") +
  facet_wrap_paginate(~gene_name, ncol = 3, nrow = 2, page = 3, scales = "free")

############################################################################################################

library(ggVennDiagram)
library(VennDiagram)

display_venn <- function(x, ...){
  library(VennDiagram)
  grid.newpage()
  venn_object <- venn.diagram(x, filename = NULL, ...)
  grid.draw(venn_object)
}

Timepoint_3hrs_DE <- Timepoint_3hrs %>% filter(abs(log2_fold_change) >= 2, q_value <= 0.05) %>% distinct(., .keep_all = T)

Timepoint_6hrs_DE <- Timepoint_6hrs %>% filter(abs(log2_fold_change) >= 2, q_value <= 0.05) %>% distinct(., .keep_all = T)

Timepoint_12hrs_DE <- Timepoint_12hrs %>% filter(abs(log2_fold_change) >= 2, q_value <= 0.05) %>% distinct(., .keep_all = T)

Timepoint_18hrs_DE <- Timepoint_18hrs %>% filter(abs(log2_fold_change) >= 2, q_value <= 0.05) %>% distinct(., .keep_all = T)

S = list(A = Timepoint_3hrs_DE$gene_name, 
         B = Timepoint_6hrs_DE$gene_name, 
         C = Timepoint_12hrs_DE$gene_name,
         D = Timepoint_18hrs_DE$gene_name)

display_venn(
  S,
  category.names = c("DE 3hrs","DE 6hrs", "DE 12hrs", "DE 18hrs"),
  # Circles
  lwd = 2,
  lty = 'blank',
  fill = c("#BFBF00", "dodgerblue4", "#B30000", "#006D2C"),
  print.mode=c("raw","percent"),
  # Numbers
  cex = 1.8,
  fontface = "bold",
  # Set names
  cat.cex = 2,
  cat.fontface = "bold",
  cat.default.pos = "outer"
)



# ggVennDiagram(S, category.names = c("DE 6hrs", "DE 12hrs", "DE 18hrs"), set_color = "black",
#               edge_size = 1) 
# 
# 
# 
# 
# scale_fill_manual(values = c("dodgerblue4", "#B30000", "#006D2C", 
#                                      "#B30000", "#006D2C", "dodgerblue4", "#B30000"))
# 
# 
# 
# 
# 
# venn <- Venn(S)
# data <- process_data(venn)
# L = as.data.frame(venn_region(data))
# L_1 = as.data.frame(L[,5])
# venn_region(data) %>%
#   mutate(name = c("6 hrs", "12 hrs", "18 hrs", "6 vs 12 hrs", "6 vs 18 hrs", "12 vs 18 hrs",
#                         "Common")) -> R
# 
# 
# R %>% ggplot() +
#   geom_sf(aes(fill=count), data = venn_region(data), colour = alpha(c("dodgerblue4", "#B30000", "#006D2C", 
#                                                  "#B30000", "#006D2C", "dodgerblue4", "#B30000"), 0.8), show.legend = F) +
#   geom_sf(aes(color=name), data = venn_region(data), show.legend = TRUE, size = 2) +
#   #geom_sf(aes(fill=Algorithms, colour = Algorithms), show.legend = T) +
#   #geom_sf(aes(color=id), data = venn_setedge(data), show.legend = F) +
#   geom_sf_text(aes(label= name), data = venn_setlabel(data), nudge_y = -10, nudge_x = -9, size = 6, fontface= "bold") +
#   geom_sf_label(aes(label = paste0(count, " (", scales::percent(count/sum(count), accuracy = 2), ")")), fontface = "bold", data = venn_region(data), size = 7, colour = "black", label.padding = unit(0.55, "lines")) +
#   scale_fill_manual(values = alpha(c("dodgerblue4", "#006D2C", "#00441B",
#                                     "#B30000", "#006D2C", "dodgerblue4", "#B30000"), 0.58)) + 
#                                       labs(title = " Overlap of DE mRNAs between \n6 hrs, 12 hrs and 18 hrs") +
#   theme_void() + theme(plot.title=element_text(family='', face='bold', size=16, hjust=0.5, vjust=1))
# 
# 
# venn <- Venn(S)
# data <- process_data(venn)

common_genes <- Timepoint_18hrs_DE[Timepoint_18hrs_DE$gene_name %in% Timepoint_12hrs_DE$gene_name,] 

#######3#
common_genes_up <- common_genes %>% filter(log2_fold_change >= 2) %>%
                dplyr::select(gene_name) %>% distinct(.)

common_genes_up_tx <- gtf[(gtf$gene_name %in% common_genes_up$gene_name), ] %>% filter(transcript_type == "protein_coding") %>% distinct(., transcript_id)

ucsc_3UTR <- read.delim2('C:/Users/matha/Documents/CircRNA_Analysis_Compile/CIRIquant_Superdepth_noStringency/circR_results/gene_ucsc_3_utr.bed', sep = '\t', header = F)

circr2_input <- ucsc_3UTR[(ucsc_3UTR$V4 %in% common_genes_up_tx$transcript_id),]


###########
common_genes_down <- common_genes %>% filter(log2_fold_change <= -2) %>%
  dplyr::select(gene_name) %>% distinct(.)



###########
write.table(common_genes, 'C:/Users/matha/Desktop/Thesis/Thesis_pics/commonGenes_12_18.txt', 
            row.names = F, col.names = F, quote = F)

write.table(common_genes_up, 'C:/Users/matha/Documents/CircRNA_Analysis_Compile/CIRIquant_Superdepth_noStringency/circR_results/commonGenes_12_18_up.txt', 
            row.names = F, col.names = F, quote = F)

write.table(common_genes_up_tx, 'C:/Users/matha/Desktop/Thesis/Thesis_pics/commonGenes_12_18_up_tx.txt', 
            row.names = F, col.names = F, quote = F)

write.table(circr2_input, 'C:/Users/matha/Documents/CircRNA_Analysis_Compile/CIRIquant_Superdepth_noStringency/circR_results/circr2_input.bed', 
            row.names = F, col.names = F, quote = F, sep = '\t')

write.table(common_genes_down, 'C:/Users/matha/Desktop/Thesis/Thesis_pics/commonGenes_12_18_down.txt', 
            row.names = F, col.names = F, quote = F)

############################################################################################################

### Number of mRNAs

mRNA_DE_table <- gene_diff_data_DE %>% group_by(sample_1) %>%
  summarise(n = n())


###########################################################################################################

## GO and KEGG

library(org.Hs.eg.db)
library(ggplot2)
library(stringr)
library(tidyr)
library(clusterProfiler)


###########################
## 3hrs
DE_3hrs <- gene_diff_data_DE %>% filter(sample_1 == "C3H") %>% transform(gene_name = strsplit(gene_name, ",")) %>% 
  unnest(gene_name)
DE_3hrs_gene_id <- gtf[(gtf$gene_name %in% DE_3hrs$gene_name),]
DE_3hrs_gene_id$gene_id <- sub("\\..*", "",DE_3hrs_gene_id$gene_id)
DE_3hrs_bitr <- bitr(DE_3hrs_gene_id$gene_id, fromType = "ENSEMBL", toType = c("SYMBOL","ENTREZID"), OrgDb = org.Hs.eg.db)

# Enrich GO
go_res_3hrs <- enrichGO(DE_3hrs_bitr$ENTREZID, keyType = "ENTREZID", 
                   org.Hs.eg.db, ont = "ALL", pvalueCutoff = 0.05, qvalueCutoff = 0.05, 
                   pAdjustMethod = "BH", readable = T)

# go_res_3hrs@result[338,3] <- sub("(.*),.*", "\\1", go_res_3hrs@result[338,3])


dotplot(go_res_3hrs,showCategory = 10, size = NULL,color = "p.adjust", font.size = 14, split = "ONTOLOGY") + 
  facet_grid(ONTOLOGY ~ ., scales = "free") +
  theme(strip.text = element_text(size = 12, face= "bold"),
        axis.text.y = element_text(size = 10, face= "bold"),
        axis.text.x = element_text(size = 10)) 



###########################
## 6hrs
DE_6hrs <- gene_diff_data_DE %>% filter(sample_1 == "C6H") %>% transform(gene_name = strsplit(gene_name, ",")) %>% 
  unnest(gene_name)
DE_6hrs_gene_id <- gtf[(gtf$gene_name %in% DE_6hrs$gene_name),]
DE_6hrs_gene_id$gene_id <- sub("\\..*", "",DE_6hrs_gene_id$gene_id)

DE_6hrs_bitr <- bitr(DE_6hrs_gene_id$gene_id, fromType = "ENSEMBL", toType = c("SYMBOL","ENTREZID"), OrgDb = org.Hs.eg.db)


# Enrich GO
go_res_6hrs <- enrichGO(DE_6hrs_bitr$ENTREZID, keyType = "ENTREZID", 
                   org.Hs.eg.db, ont = "ALL", pvalueCutoff = 0.05, qvalueCutoff = 0.05, 
                   pAdjustMethod = "BH", readable = T)

go_res_6hrs@result[368,3] <- sub("(.*),.*", "\\1", go_res_6hrs@result[368,3])


dotplot(go_res_6hrs,showCategory = 10, size = NULL,color = "p.adjust", font.size = 18, split = "ONTOLOGY") + 
  facet_grid(ONTOLOGY ~ ., scales = "free") +
  theme(strip.text = element_text(size = 12, face= "bold"),
        axis.text.y = element_text(size = 10, face= "bold"),
        axis.text.x = element_text(size = 10)) 



###########################
## 12hrs
DE_12hrs <- gene_diff_data_DE %>% filter(sample_1 == "C12H") %>% transform(gene_name = strsplit(gene_name, ",")) %>% 
  unnest(gene_name)
DE_12hrs_gene_id <- gtf[(gtf$gene_name %in% DE_12hrs$gene_name),]
DE_12hrs_gene_id$gene_id <- sub("\\..*", "",DE_12hrs_gene_id$gene_id)

DE_12hrs_bitr <- bitr(DE_12hrs_gene_id$gene_id, fromType = "ENSEMBL", toType = c("SYMBOL","ENTREZID"), OrgDb = org.Hs.eg.db)


# Enrich GO
go_res_12hrs <- enrichGO(DE_12hrs_bitr$ENTREZID, keyType = "ENTREZID", 
                   org.Hs.eg.db, ont = "ALL", pvalueCutoff = 0.05, qvalueCutoff = 0.05, 
                   pAdjustMethod = "BH", readable = T)

go_res_12hrs@result[522,3] <- sub("(.*),.*", "\\1", go_res_12hrs@result[522,3])

dotplot(go_res_12hrs,showCategory = 10, size = NULL,color = "p.adjust", font.size = 14, split = "ONTOLOGY") + 
  facet_grid(ONTOLOGY ~ ., scales = "free") +
  theme(strip.text = element_text(size = 12, face= "bold"),
        text = element_text(size = 20),
        axis.text.y = element_text(size = 12, face= "bold"),
        axis.text.x = element_text(size = 10)) 


###########################
## 18hrs

DE_18hrs <- gene_diff_data_DE %>% filter(sample_1 == "C18H") %>% transform(gene_name = strsplit(gene_name, ",")) %>% 
  unnest(gene_name)
DE_18hrs_gene_id <- gtf[(gtf$gene_name %in% DE_18hrs$gene_name),]
DE_18hrs_gene_id$gene_id <- sub("\\..*", "",DE_18hrs_gene_id$gene_id)

DE_18hrs_bitr <- bitr(DE_18hrs_gene_id$gene_id, fromType = "ENSEMBL", toType = c("SYMBOL","ENTREZID"), OrgDb = org.Hs.eg.db)


# Enrich GO
go_res_18hrs <- enrichGO(DE_18hrs_bitr$ENTREZID, keyType = "ENTREZID", 
                   org.Hs.eg.db, ont = "ALL", pvalueCutoff = 0.05, qvalueCutoff = 0.05, 
                   pAdjustMethod = "BH", readable = T)


dotplot(go_res_18hrs,showCategory = 10, size = NULL,color = "p.adjust", font.size = 14, split = "ONTOLOGY") + 
  facet_grid(ONTOLOGY ~ ., scales = "free") +
  theme(strip.text = element_text(size = 14, face= "bold"),
        axis.text.y = element_text(size = 10, face= "bold"),
        strip.text.x = element_text(size = 10)) 


######################
## Build go table

go_res_3hrs@result$Timepoint <- "03hrs"
go_res_6hrs@result$Timepoint <- "06hrs"
go_res_12hrs@result$Timepoint <- "12hrs"
go_res_18hrs@result$Timepoint <- "18hrs"

# Combine all results into one table
go_table <- rbind(go_res_3hrs@result %>% arrange(p.adjust) %>% filter(as.numeric(p.adjust) <= 0.05),
                  go_res_6hrs@result %>% arrange(p.adjust) %>% filter(as.numeric(p.adjust) <= 0.05),
                  go_res_12hrs@result %>% arrange(p.adjust) %>% filter(as.numeric(p.adjust) <= 0.05),
                  go_res_18hrs@result %>% arrange(p.adjust) %>% filter(as.numeric(p.adjust) <= 0.05))

# Arrange them by p.adjust value
go_table_adj <- go_table %>% arrange(p.adjust)

# Take top 10 from each Ontology
go_table_adj_top10_BP <- go_table_adj %>% filter(ONTOLOGY == 'BP') %>% distinct(., Description) %>% slice(1:10)
go_table_adj_top10_MF <- go_table_adj %>% filter(ONTOLOGY == 'MF') %>% distinct(., Description) %>% slice(1:10)
go_table_adj_top10_CC <- go_table_adj %>% filter(ONTOLOGY == 'CC') %>% distinct(., Description) %>% slice(1:10)


###########################################################################################################
## pheatmap bold names

make_bold_names <- function(mat, rc_fun, rc_names) {
  bold_names <- rc_fun(mat)
  ids <- rc_names %>% match(rc_fun(mat))
  ids %>%
    walk(
      function(i)
        bold_names[i] <<-
        bquote(bold(.(rc_fun(mat)[i]))) %>%
        as.expression()
    )
  bold_names
}

###################################################################################################################
## Gene LogFC

## Substituting the min and max volumes
gene_diff_data_logFC <- gene_diff_data
for (i in c(3, 6, 12, 18)) {
  
  ## Substituting INF for max number of that timepoint
  gene_diff_data_logFC[(gene_diff_data_logFC$sample_1 %in% paste0("C", i, "H")) & (gene_diff_data_logFC$log2_fold_change == "Inf"),][7] <- max(gene_diff_data_logFC$log2_fold_change[is.finite(gene_diff_data_logFC$log2_fold_change) & (gene_diff_data_logFC$sample_1 %in% paste0("C", i, "H"))])
  
  ## Substituting -INF for min number of that timepoint
  gene_diff_data_logFC[(gene_diff_data_logFC$sample_1 %in% paste0("C", i, "H")) & (gene_diff_data_logFC$log2_fold_change == "-Inf"),][7] <- min(gene_diff_data_logFC$log2_fold_change[is.finite(gene_diff_data_logFC$log2_fold_change) & (gene_diff_data_logFC$sample_1 %in% paste0("C", i, "H"))])
  
}


## Create a long df
GeneLogFC <- gene_diff_data_logFC %>%
  select(c("sample_1", "log2_fold_change", "gene_name", "gene_id")) %>%
  transform(gene_name = strsplit(gene_name, ",")) %>% 
  unnest(gene_name)

## Convert to wider df
GeneLogFC_mx <- GeneLogFC %>%
  pivot_wider(names_from = sample_1, values_from = log2_fold_change) %>%
  as.data.frame()



#### Extract DE genes from matrix
## DE gene list
DE <- gene_diff_data %>% filter(abs(log2_fold_change) >= 2, q_value <= 0.05) %>% select(gene_id)



DE_geneLogFC_mx <- GeneLogFC_mx %>%
  filter(gene_id %in% DE$gene_id)



# rownames(DE_geneLogFC_mx) <- DE_geneLogFC_mx$gene_name
# 
# DE_geneLogFC_mx <- DE_geneLogFC_mx %>% select(c(3:6))
# colnames(DE_geneLogFC_mx) <- c("3hrs", "6hrs", "12hrs", "18hrs")


###########################################################################################################
## HEATMAP FOR BP

go_table_manual <- go_table_adj[go_table_adj$Description %in% go_table_adj_top10_BP$Description,]
go_table_manual <- go_table_manual %>% select(Description, Timepoint, p.adjust)
go_table_manual$Timepoint <- factor(go_table_manual$Timepoint, levels = c("03hrs", "06hrs", "12hrs", "18hrs"))
go_table_manual_2 <- go_table_manual %>% pivot_wider(names_from = Timepoint, values_from = p.adjust)
go_table_manual_2[is.na(go_table_manual_2)] <- 1
go_table_manual_2 <- go_table_manual_2 %>% column_to_rownames("Description")
go_table_manual_2 <- go_table_manual_2[, rev(colnames(go_table_manual_2))] %>% select(1, 4, 3, 2)

## Heatmap 

pheatmap(go_table_manual_2 %>% as.matrix(.), 
                          #annotation_row = rownames(go_table_manual_2), 
                          cluster_rows = F,
                          cluster_cols = F, clustering_distance_rows = "euclidean",
                          clustering_distance_cols = "euclidean", clustering_method = "complete",
                          fontsize = 12, annotation_names_row = F,
                          fontsize_row = 15, fontsize_col = 20, annotation_names_col = F,
                          #color=hmcols,
                          scale="none", show_rownames = T,
                          main = "Top Significant GO BP terms", border_color='grey60',
         labels_row = make_bold_names(go_table_manual_2, rownames, c("defense response to virus", "regulation of innate immune response", "negative regulation of viral genome replication")),
         labels_col = make_bold_names(go_table_manual_2, colnames, c("3hrs", "6hrs", "12hrs", "18hrs")))

## We take "defense response to virus", "cytokine-mediated signaling pathway", "negative regulation of viral genome replication" after checking the go_table_adj

############################################################################
############################################################################

## Heatmap of Associated genes for BP

# Build database
bitr_db <- rbind(DE_3hrs_bitr, DE_6hrs_bitr, DE_12hrs_bitr, DE_18hrs_bitr)
gene_id_db <- rbind(DE_3hrs_gene_id, DE_6hrs_gene_id, DE_12hrs_gene_id, DE_18hrs_gene_id)

# GO terms of choice for BP
terms_list <- c("regulation of innate immune response", "defense response to virus", "negative regulation of viral genome replication")

# Unnest all genes from the GO terms of choice
gene_list <- go_table_adj %>% filter(Description %in% terms_list) %>% 
  transform(geneID = strsplit(geneID, "/")) %>% 
  unnest(geneID)

# Get symbol OF gene and then get ID
gene_bitr <- bitr_db[(bitr_db$SYMBOL %in% gene_list$geneID),]
gene_symbol_id <- gene_id_db[(gene_id_db$gene_id %in% gene_bitr$ENSEMBL),] %>% select(gene_id, gene_name) %>% unique()

# Merged to a big table
merged <- merge(gene_symbol_id, gene_bitr, by.x = "gene_id", by.y = "ENSEMBL") %>% unique()
merged_GO <- merge(merged, gene_list, by.x="SYMBOL", by.y="geneID")


## Get the logFC values
gene_diff_data_BP <- DE_geneLogFC_mx %>%
  filter(gene_name %in% gene_symbol_id$gene_name)


gene_diff_data_merged <- merge(gene_diff_data_BP, merged_GO, by.x = "gene_name", by.y = "gene_name")
gene_diff_data_BP_1 <- gene_diff_data_merged %>% select(gene_name, C3H, C6H, C12H, C18H, Description) %>% unique()
colnames(gene_diff_data_BP_1) <- c('gene_name', '3hrs', '6hrs', '12hrs', '18hrs', 'Description')


# Get GO term with expression value matrix
gene_diff_data_BP_2 <- gene_diff_data_BP_1 %>% select('gene_name', '3hrs', '6hrs', '12hrs', '18hrs') %>% unique()
rownames(gene_diff_data_BP_2) <- gene_diff_data_BP_2$gene_name
gene_diff_data_BP_2 <- gene_diff_data_BP_2%>% select(-gene_name)

# Get annotation terms for each gene
annotation <-  gene_diff_data_BP_1 %>% select(gene_name, Description) %>% as.data.frame() %>% unique() 
annotation <- annotation %>%
  pivot_wider(names_from = "Description", values_from = "Description") %>% 
  as.data.frame() %>%
  column_to_rownames("gene_name") %>% rev()
colnames(annotation) <- c("GO BP Term 2", "GO BP Term 1", "GO BP Term 3")



# Draw pheatmap
pheatmap(gene_diff_data_BP_2,
         annotation_row = annotation[,c(3,1,2)], 
         cluster_rows = T,
         cluster_cols = F, clustering_distance_rows = "euclidean",
         clustering_distance_cols = "euclidean", clustering_method = "complete",
         fontsize = 20, annotation_names_row = F,
         fontsize_row = 10, fontsize_col = 24, annotation_names_col = F,
         #color=hmcols,
         scale="none", show_rownames = T,
         main = "GO BP terms with Gene log2FC expression", border_color=NA,
         labels_col = make_bold_names(gene_diff_data_BP_2, colnames, c("3hrs", "6hrs", "12hrs", "18hrs")))


############################################################################
############################################################################
## HEATMAP FOR MF

go_table_manual <- go_table_adj[go_table_adj$Description %in% go_table_adj_top10_MF$Description,]
go_table_manual <- go_table_manual %>% select(Description, Timepoint, p.adjust)
go_table_manual$Timepoint <- factor(go_table_manual$Timepoint, levels = c("03hrs", "06hrs", "12hrs", "18hrs"))
go_table_manual_2 <- go_table_manual %>% pivot_wider(names_from = Timepoint, values_from = p.adjust)
go_table_manual_2[is.na(go_table_manual_2)] <- 1
go_table_manual_2 <- go_table_manual_2 %>% column_to_rownames("Description")
go_table_manual_2 <- go_table_manual_2[, rev(colnames(go_table_manual_2))] %>% select("03hrs", "06hrs", "12hrs", "18hrs")

## Heatmap 

pheatmap(go_table_manual_2 %>% as.matrix(.), 
         #annotation_row = rownames(go_table_manual_2), 
         cluster_rows = F,
         cluster_cols = F, clustering_distance_rows = "euclidean",
         clustering_distance_cols = "euclidean", clustering_method = "complete",
         fontsize = 12, annotation_names_row = F,
         fontsize_row = 15, fontsize_col = 20, annotation_names_col = F,
         #color=hmcols,
         scale="none", show_rownames = T,
         main = "Top Significant GO MF terms", border_color='grey60',
         labels_row = make_bold_names(go_table_manual_2, rownames, c("double-stranded RNA binding", "type I interferon receptor binding", "chemokine activity")),
         labels_col = make_bold_names(go_table_manual_2, colnames, c("3hrs", "6hrs", "12hrs", "18hrs")))

############################################################################
############################################################################

## Heatmap of Associated genes for MF

# Build database
bitr_db <- rbind(DE_3hrs_bitr, DE_6hrs_bitr, DE_12hrs_bitr, DE_18hrs_bitr)
gene_id_db <- rbind(DE_3hrs_gene_id, DE_6hrs_gene_id, DE_12hrs_gene_id, DE_18hrs_gene_id)

# GO terms of choice for MF
terms_list <- c("double-stranded RNA binding", "type I interferon receptor binding", "chemokine activity")

# Unnest all genes from the GO terms of choice
gene_list <- go_table_adj %>% filter(Description %in% terms_list) %>% 
  transform(geneID = strsplit(geneID, "/")) %>% 
  unnest(geneID)

# Get symbol OF gene and then get ID
gene_bitr <- bitr_db[(bitr_db$SYMBOL %in% gene_list$geneID),]
gene_symbol_id <- gene_id_db[(gene_id_db$gene_id %in% gene_bitr$ENSEMBL),] %>% select(gene_id, gene_name) %>% unique()

# Merged to a big table
merged <- merge(gene_symbol_id, gene_bitr, by.x = "gene_id", by.y = "ENSEMBL") %>% unique()
merged_GO <- merge(merged, gene_list, by.x="SYMBOL", by.y="geneID")


## Get the logFC values
gene_diff_data_MF <- DE_geneLogFC_mx %>%
  filter(gene_name %in% gene_symbol_id$gene_name)


gene_diff_data_merged <- merge(gene_diff_data_MF, merged_GO, by.x = "gene_name", by.y = "gene_name")
gene_diff_data_MF_1 <- gene_diff_data_merged %>% select(gene_name, C3H, C6H, C12H, C18H, Description) %>% unique()
colnames(gene_diff_data_MF_1) <- c('gene_name', '3hrs', '6hrs', '12hrs', '18hrs', 'Description')


# Get GO term with expression value matrix
gene_diff_data_MF_2 <- gene_diff_data_MF_1 %>% select('gene_name', '3hrs', '6hrs', '12hrs', '18hrs') %>% unique()
rownames(gene_diff_data_MF_2) <- gene_diff_data_MF_2$gene_name
gene_diff_data_MF_2 <- gene_diff_data_MF_2%>% select(-gene_name)

# Get annotation terms for each gene
annotation <-  gene_diff_data_MF_1 %>% select(gene_name, Description) %>% as.data.frame() %>% unique() 
annotation <- annotation %>%
  pivot_wider(names_from = "Description", values_from = "Description") %>% 
  as.data.frame() %>%
  column_to_rownames("gene_name") %>% rev()
colnames(annotation) <- c("GO MF Term 1", "GO MF Term 3", "GO MF Term 2")



# Draw pheatmap
pheatmap(gene_diff_data_MF_2,
         annotation_row = annotation[,c(2,3,1)], 
         cluster_rows = T,
         cluster_cols = F, clustering_distance_rows = "euclidean",
         clustering_distance_cols = "euclidean", clustering_method = "complete",
         fontsize = 20, annotation_names_row = F,
         fontsize_row = 10, fontsize_col = 24, annotation_names_col = F,
         #color=hmcols,
         scale="none", show_rownames = T,
         main = "GO MF terms with Gene log2FC expression", border_color=NA,
         labels_col = make_bold_names(gene_diff_data_MF_2, colnames, c("3hrs", "6hrs", "12hrs", "18hrs")))

############################################################################
############################################################################

## HEATMAP FOR CC

go_table_manual <- go_table_adj[go_table_adj$Description %in% go_table_adj_top10_CC$Description,]
go_table_manual <- go_table_manual %>% select(Description, Timepoint, p.adjust)
go_table_manual$Timepoint <- factor(go_table_manual$Timepoint, levels = c("06hrs", "12hrs", "18hrs"))
go_table_manual_2 <- go_table_manual %>% pivot_wider(names_from = Timepoint, values_from = p.adjust)
go_table_manual_2[is.na(go_table_manual_2)] <- 1
go_table_manual_2 <- go_table_manual_2 %>% column_to_rownames("Description")
go_table_manual_2 <- go_table_manual_2[, rev(colnames(go_table_manual_2))] %>% select("06hrs", "12hrs", "18hrs")


## Heatmap 

pheatmap(go_table_manual_2 %>% as.matrix(.), 
         #annotation_row = rownames(go_table_manual_2), 
         cluster_rows = F,
         cluster_cols = F, clustering_distance_rows = "euclidean",
         clustering_distance_cols = "euclidean", clustering_method = "complete",
         fontsize = 12, annotation_names_row = F,
         fontsize_row = 15, fontsize_col = 20, annotation_names_col = F,
         #color=hmcols,
         scale="none", show_rownames = T,
         main = "Top Significant go CC terms", border_color='grey60',
         labels_row = make_bold_names(go_table_manual_2, rownames, c("MCM complex", "CMG complex", "canonical inflammasome complex")),
         labels_col = make_bold_names(go_table_manual_2, colnames, c("3hrs", "6hrs", "12hrs", "18hrs")))


## We take "P-body", "inflammasome complex", "CENP-A containing nucleosome" after checking the go_table_adj



############################################################################
############################################################################

## Heatmap of Associated genes for CC

# Build database
bitr_db <- rbind(DE_3hrs_bitr, DE_6hrs_bitr, DE_12hrs_bitr, DE_18hrs_bitr)
gene_id_db <- rbind(DE_3hrs_gene_id, DE_6hrs_gene_id, DE_12hrs_gene_id, DE_18hrs_gene_id)

# GO terms of choice for CC
terms_list <- c("MCM complex", "CMG complex", "canonical inflammasome complex")

# Unnest all genes from the GO terms of choice
gene_list <- go_table_adj %>% filter(Description %in% terms_list) %>% 
  transform(geneID = strsplit(geneID, "/")) %>% 
  unnest(geneID)

# Get symbol OF gene and then get ID
gene_bitr <- bitr_db[(bitr_db$SYMBOL %in% gene_list$geneID),]
gene_symbol_id <- gene_id_db[(gene_id_db$gene_id %in% gene_bitr$ENSEMBL),] %>% select(gene_id, gene_name) %>% unique()

# Merged to a big table
merged <- merge(gene_symbol_id, gene_bitr, by.x = "gene_id", by.y = "ENSEMBL") %>% unique()
merged_GO <- merge(merged, gene_list, by.x="SYMBOL", by.y="geneID")


## Get the logFC values
gene_diff_data_CC <- DE_geneLogFC_mx %>%
  filter(gene_name %in% gene_symbol_id$gene_name)


gene_diff_data_merged <- merge(gene_diff_data_CC, merged_GO, by.x = "gene_name", by.y = "gene_name")
gene_diff_data_CC_1 <- gene_diff_data_merged %>% select(gene_name, C3H, C6H, C12H, C18H, Description) %>% unique()
colnames(gene_diff_data_CC_1) <- c('gene_name', '3hrs', '6hrs', '12hrs', '18hrs', 'Description')


# Get GO term with expression value matrix
gene_diff_data_CC_2 <- gene_diff_data_CC_1 %>% select('gene_name', '3hrs', '6hrs', '12hrs', '18hrs') %>% unique()
rownames(gene_diff_data_CC_2) <- gene_diff_data_CC_2$gene_name
gene_diff_data_CC_2 <- gene_diff_data_CC_2%>% select(-gene_name)

# Get annotation terms for each gene
annotation <-  gene_diff_data_CC_1 %>% select(gene_name, Description) %>% as.data.frame() %>% unique() 
annotation <- annotation %>%
  pivot_wider(names_from = "Description", values_from = "Description") %>% 
  as.data.frame() %>%
  column_to_rownames("gene_name") %>% rev()
colnames(annotation) <- c("GO CC Term 1", "GO CC Term 2", "GO CC Term 3")



# Draw pheatmap
pheatmap(gene_diff_data_CC_2,
         annotation_row = annotation[,c(3,2,1)], 
         cluster_rows = T,
         cluster_cols = F, clustering_distance_rows = "euclidean",
         clustering_distance_cols = "euclidean", clustering_method = "complete",
         fontsize = 12, annotation_names_row = F,
         fontsize_row = 12, fontsize_col = 20, annotation_names_col = F,
         #color=hmcols,
         scale="none", show_rownames = T,
         main = "GO CC terms with Gene log2FC expression", border_color=NA,
         labels_col = make_bold_names(gene_diff_data_CC_2, colnames, c("3hrs", "6hrs", "12hrs", "18hrs")))




############################################################################
############################################################################

terms_list <- c("defense response to virus", "cytokine-mediated signaling pathway", "negative regulation of viral genome replication")
terms_list <- c("guanyl nucleotide binding", "tumor necrosis factor receptor superfamily binding", "chemokine activity")
terms_list <- c("P-body", "inflammasome complex", "CENP-A containing nucleosome")
#####################################################################################################
#####################################################################################################
#####################################################################################################
###########################
# Enrich KEGG
library(R.utils)

R.utils::setOption("clusterProfiler.download.method","auto")

## 3hrs
# Enrich KEGG

keg_res_3hrs <- enrichKEGG(DE_3hrs_bitr$ENTREZID, organism= "hsa", keyType= "kegg", 
                      pvalueCutoff = 0.01, pAdjustMethod = "BH")
keg_res_3hrs <- setReadable(keg_res_3hrs, OrgDb = org.Hs.eg.db, keyType="ENTREZID")

dotplot(keg_res_3hrs, showCategory = 10, size = NULL,color = "p.adjust", font.size = 14) + 
  # facet_grid(ONTOLOGY ~ ., scales = "free") +
  theme(axis.text.x = element_text(size = 14, face= "bold"),
        axis.text.y = element_text(size = 14, face= "bold")) 



###########################
## 6hrs

keg_res_6hrs <- enrichKEGG(DE_6hrs_bitr$ENTREZID, organism= "hsa", keyType= "kegg", 
                      pvalueCutoff = 0.01, pAdjustMethod = "BH")
keg_res_6hrs <- setReadable(keg_res_6hrs, OrgDb = org.Hs.eg.db, keyType="ENTREZID")

dotplot(keg_res_6hrs, showCategory = 10, size = NULL,color = "p.adjust", font.size = 14) + 
  # facet_grid(ONTOLOGY ~ ., scales = "free") +
  theme(axis.text.x = element_text(size = 14, face= "bold"),
        axis.text.y = element_text(size = 14, face= "bold")) 



###########################
## 12hrs
keg_res_12hrs <- enrichKEGG(DE_12hrs_bitr$ENTREZID, organism= "hsa", keyType= "kegg", 
                      pvalueCutoff = 0.01, pAdjustMethod = "BH")
keg_res_12hrs <- setReadable(keg_res_12hrs, OrgDb = org.Hs.eg.db, keyType="ENTREZID")

dotplot(keg_res_12hrs, showCategory = 10, size = NULL,color = "p.adjust", font.size = 14) + 
  # facet_grid(ONTOLOGY ~ ., scales = "free") +
  theme(axis.text.x = element_text(size = 14, face= "bold"),
        axis.text.y = element_text(size = 14, face= "bold")) 



###########################
## 18hrs
keg_res_18hrs <- enrichKEGG(DE_18hrs_bitr$ENTREZID, organism= "hsa", keyType= "kegg", 
                      pvalueCutoff = 0.01, pAdjustMethod = "BH")
keg_res_18hrs <- setReadable(keg_res_18hrs, OrgDb = org.Hs.eg.db, keyType="ENTREZID")

dotplot(keg_res_18hrs, showCategory = 10, size = NULL,color = "p.adjust", font.size = 14) + 
  # facet_grid(ONTOLOGY ~ ., scales = "free") +
  theme(axis.text.x = element_text(size = 14, face= "bold"),
        axis.text.y = element_text(size = 14, face= "bold")) 

##################################################################################################
######################
## Build kegg table

keg_res_3hrs@result$Timepoint <- "3hrs"
keg_res_6hrs@result$Timepoint <- "6hrs"
keg_res_12hrs@result$Timepoint <- "12hrs"
keg_res_18hrs@result$Timepoint <- "18hrs"

# Combine all results into one table
kegg_table <- rbind(keg_res_3hrs@result %>% arrange(p.adjust) %>% filter(as.numeric(p.adjust) <= 0.05),
                  keg_res_6hrs@result %>% arrange(p.adjust) %>% filter(as.numeric(p.adjust) <= 0.05),
                  keg_res_12hrs@result %>% arrange(p.adjust) %>% filter(as.numeric(p.adjust) <= 0.05),
                  keg_res_18hrs@result %>% arrange(p.adjust) %>% filter(as.numeric(p.adjust) <= 0.05))

# Arrange them by p.adjust value
kegg_table_adj <- kegg_table %>% arrange(p.adjust) %>% distinct(., Description, .keep_all = T) %>% slice(1:10)


###########################################################################################################
##################################################################################################
## HEATMAP FOR KEGG

kegg_table_manual <- kegg_table[kegg_table$Description %in% kegg_table_adj$Description,]
kegg_table_manual <- kegg_table_manual %>% select(Description, Timepoint, p.adjust)
kegg_table_manual$Timepoint <- factor(kegg_table_manual$Timepoint, levels = c("3hrs", "6hrs", "12hrs", "18hrs"))
kegg_table_manual_2 <- kegg_table_manual %>% pivot_wider(names_from = Timepoint, values_from = p.adjust)
kegg_table_manual_2[is.na(kegg_table_manual_2)] <- 1
kegg_table_manual_2 <- kegg_table_manual_2 %>% column_to_rownames("Description")
kegg_table_manual_2 <- kegg_table_manual_2[, rev(colnames(kegg_table_manual_2))] %>% select("3hrs", "6hrs", "12hrs", "18hrs")


## Heatmap 

pheatmap(kegg_table_manual_2 %>% as.matrix(.), 
         #annotation_row = rownames(kegg_table_manual_2), 
         cluster_rows = F,
         cluster_cols = F, clustering_distance_rows = "euclidean",
         clustering_distance_cols = "euclidean", clustering_method = "complete",
         fontsize = 12, annotation_names_row = F,
         fontsize_row = 15, fontsize_col = 20, annotation_names_col = F,
         #color=hmcols,
         scale="none", show_rownames = T,
         main = "Top 10 Significant KEGG terms", border_color='grey60',
         labels_row = make_bold_names(kegg_table_manual_2, rownames, c("Influenza A", "Cytokine-cytokine receptor interaction", "RIG-I-like receptor signaling pathway")),
         labels_col = make_bold_names(kegg_table_manual_2, colnames, c("3hrs", "6hrs", "12hrs", "18hrs")))


## We take "P-body", "inflammasome complex", "CENP-A containing nucleosome" after checking the go_table_adj



############################################################################
############################################################################

## Heatmap of Associated genes for KEGG

# Build database
bitr_db <- rbind(DE_3hrs_bitr, DE_6hrs_bitr, DE_12hrs_bitr, DE_18hrs_bitr)
gene_id_db <- rbind(DE_3hrs_gene_id, DE_6hrs_gene_id, DE_12hrs_gene_id, DE_18hrs_gene_id)

# GO terms of choice for CC
terms_list <- c("Influenza A", "Cytokine-cytokine receptor interaction", "RIG-I-like receptor signaling pathway")

# Unnest all genes from the GO terms of choice
gene_list <- kegg_table %>% filter(Description %in% terms_list) %>% 
  transform(geneID = strsplit(geneID, "/")) %>% 
  unnest(geneID)

# Get symbol OF gene and then get ID
gene_bitr <- bitr_db[(bitr_db$SYMBOL %in% gene_list$geneID),]
gene_symbol_id <- gene_id_db[(gene_id_db$gene_id %in% gene_bitr$ENSEMBL),] %>% select(gene_id, gene_name) %>% unique()

# Merged to a big table
merged <- merge(gene_symbol_id, gene_bitr, by.x = "gene_id", by.y = "ENSEMBL") %>% unique()
merged_KEGG <- merge(merged, gene_list, by.x="SYMBOL", by.y="geneID")


## Get the logFC values
gene_diff_data <- DE_geneLogFC_mx %>%
  filter(gene_name %in% gene_symbol_id$gene_name)


gene_diff_data_merged <- merge(gene_diff_data, merged_KEGG, by.x = "gene_name", by.y = "gene_name")
gene_diff_data_1 <- gene_diff_data_merged %>% select(gene_name, C3H, C6H, C12H, C18H, Description) %>% unique()
colnames(gene_diff_data_1) <- c('gene_name', '3hrs', '6hrs', '12hrs', '18hrs', 'Description')


# Get GO term with expression value matrix
gene_diff_data_2 <- gene_diff_data_1 %>% select('gene_name', '3hrs', '6hrs', '12hrs', '18hrs') %>% unique()
rownames(gene_diff_data_2) <- gene_diff_data_2$gene_name
gene_diff_data_2 <- gene_diff_data_2%>% select(-gene_name)

# Get annotation terms for each gene
annotation <-  gene_diff_data_1 %>% select(gene_name, Description) %>% as.data.frame() %>% unique() 
annotation <- annotation %>%
  pivot_wider(names_from = "Description", values_from = "Description") %>% 
  as.data.frame() %>%
  column_to_rownames("gene_name") %>% rev()
colnames(annotation) <- c("KEGG Term 3", "KEGG Term 2", "KEGG Term 1")



# Draw pheatmap
pheatmap(gene_diff_data_2,
         annotation_row = annotation, 
         cluster_rows = T,
         cluster_cols = F, clustering_distance_rows = "euclidean",
         clustering_distance_cols = "euclidean", clustering_method = "complete",
         fontsize = 20, annotation_names_row = F,
         fontsize_row = 10, fontsize_col = 24, annotation_names_col = F,
         #color=hmcols,
         scale="none", show_rownames = T,
         main = "KEGG terms with Gene log2FC expression", border_color=NA,
         labels_col = make_bold_names(gene_diff_data_2, colnames, c("3hrs", "6hrs", "12hrs", "18hrs")))






##################################################################################################
### cOMMON 12hrs and 18hrs

DE_12_18hrs <- DE_18hrs[DE_18hrs$SYMBOL %in% DE_12hrs$SYMBOL,]


# Enrich GO terms
GOenrich <- enrichGO(DE_12_18hrs$ENTREZID, keyType = "ENTREZID", org.Hs.eg.db, ont = "All", 
                     pvalueCutoff = 0.05, pAdjustMethod = "BH", qvalueCutoff = 0.05)


dotplot(GOenrich,showCategory = 10, size = NULL,font.size = 8, color = "p.adjust",
        title = "Linear 12hrs and 18hrs GO  Enrichment", split = "ONTOLOGY") + facet_grid_sc(ONTOLOGY ~ ., scales = "free")



# Enrich KEGG terms
keg_res <- enrichKEGG(DE_12_18hrs$ENTREZID, organism="hsa", keyType="kegg", pvalueCutoff = 0.05, 
                      pAdjustMethod = "BH", qvalueCutoff = 0.05)

# Dotplot for KEGG
dotplot(keg_res,showCategory = 10, size = NULL,font.size = 8, title = "Linear 12hrs and 18hrs KEGG enrichment") 


############################################################################################################

### CeRNA 12hrs and 18hrs


ceRNA_enrichment <- read.delim2("C:/Users/matha/Documents/CircRNA_ceRNA/CIRIquant_Stringent/Gene_targets_12_18hrs.txt", sep = "\t", header = F)

ceRNA_enrichment_linear <- clusterProfiler::bitr(ceRNA_enrichment$V1, fromType = "SYMBOL", toType = "ENTREZID", OrgDb = org.Hs.eg.db)


# Enrich GO terms
GOenrich <- enrichGO(ceRNA_enrichment_linear$ENTREZID, keyType = "ENTREZID", org.Hs.eg.db, ont = "All", 
                     pvalueCutoff = 0.01, pAdjustMethod = "BH", qvalueCutoff = 0.05)


dotplot(GOenrich,showCategory = 10, size = NULL,font.size = 8, 
        title = "CeRNA Linear 12hrs and 18hrs GO  Enrichment", split = "ONTOLOGY") + facet_grid_sc(ONTOLOGY ~ ., scales = "free")



# Enrich KEGG terms
keg_res <- enrichKEGG(ceRNA_enrichment_linear$ENTREZID, organism="hsa", keyType="kegg", pvalueCutoff = 0.05, 
                      pAdjustMethod = "BH", qvalueCutoff = 0.05)

# Dotplot for KEGG
dotplot(keg_res,showCategory = 10, size = NULL,font.size = 8, title = "CeRNA Linear 12hrs and 18hrs KEGG enrichment") 





disp<-dispersionPlot(genes(cuff)) +
  theme(text = element_text(size=20))
disp


genes.MDS<-MDSplot(genes(cuff)) + theme(text = element_text(size=20))
genes.MDS

