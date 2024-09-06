
## Circular
library(ggplot2)
library(dplyr)
library(rtracklayer)
library(pheatmap)

## Read and combine data from different timepoints
setwd('C:/Users/matha/Documents/CircRNA_Analysis_Compile/CIRIquant_Superdepth_noStringency/Allhrs_Superdepth_NoStringency')

Circ_3hrs <- read.table(paste0("New_design_circRNA_03hrs.tsv", sep=""), sep=",", header = T) %>%
  mutate(Timepoint = "3hrs")
Circ_6hrs <- read.table(paste0("New_design_circRNA_06hrs.tsv", sep=""), sep=",", header = T) %>%
  mutate(Timepoint = "6hrs")
Circ_12hrs <- read.table(paste0("New_design_circRNA_12hrs.tsv", sep=""), header = T, sep = ",") %>% 
  mutate(Timepoint = "12hrs")
Circ_18hrs <- read.table(paste0("New_design_circRNA_18hrs.tsv", sep=""), sep=",", header = T) %>%
  mutate(Timepoint = "18hrs")

## Combine all timepoints into one df
Circ_SD <- rbind(Circ_3hrs, Circ_6hrs, Circ_12hrs, Circ_18hrs)
colnames(Circ_SD)[colnames(Circ_SD) == 'X'] <- 'CircRNA_ID'

## Get other informations from .gtf
Datasets <- paste0("SRR54446", 59:74, sep="")

GTF = list()
for (i in Datasets){
  GTF[[i]] <- rtracklayer::import(paste0("../", i, ".gtf", sep="")) %>% as.data.frame() %>%
    mutate(Dataset = paste0(i))
}

Circ_gtf <- do.call(rbind, GTF)
colnames(Circ_gtf)[colnames(Circ_gtf) == 'circ_id'] <- 'CircRNA_ID'

## Merge info
Circ_SD_Ns <- merge(Circ_SD, Circ_gtf, by="CircRNA_ID") %>% dplyr::select(CircRNA_ID, logFC, logCPM, LR, PValue, DE, FDR, junc_ratio, Timepoint,
                                                                          strand, circ_type, gene_id, gene_name, gene_type, Dataset)

## Rearrange row order
x <- c("3hrs", "6hrs", "12hrs", "18hrs")

Circ_SD_Ns <- Circ_SD_Ns %>% 
  mutate(Timepoint =  factor(Timepoint, levels = x)) %>% 
  arrange(Timepoint) %>% distinct(CircRNA_ID, logFC, FDR, Timepoint, .keep_all = T)

## Get DE
Circ_SD_DE <- Circ_SD_Ns %>% filter(round(abs(logFC),1) >= 2, round(FDR,2) <= 0.05)


## DE for timepoints

Circ_DE6hrs <- Circ_SD_DE %>% filter(Timepoint == "6hrs")
Circ_DE12hrs <- Circ_SD_DE %>% filter(Timepoint == "12hrs")
Circ_DE18hrs <- Circ_SD_DE %>% filter(Timepoint == "18hrs")


##################################################################################

### Number of circRNAs

CircRNA_DE_table <- Circ_SD_DE %>% group_by(Timepoint) %>% summarise(n = n())

CircRNA_table <- Circ_SD_Ns %>% group_by(Timepoint) %>%
  summarise(n = n())

#################################################################################
### Get Annotation

## Combine both database

# CircAtlas
circAtlas <- read.table("C:/Users/matha/Documents/circAtlas/human_bed_v3.0.txt", sep = "\t", header = T)
circAtlas <- circAtlas %>% mutate(circRNA_id = paste0(Chro, ":", Start, "|", End, Strand)) %>% dplyr::select(circRNA_id)

# CircBase
circBase <- read.table("C:/Users/matha/Documents/circBase/circBase_circRNA.txt", sep = "\t", header = F) %>% filter(grepl("hg38", V6))
circBase <- circBase %>% mutate(circRNA_id = paste0(V1, ":", V2 + 1, "|", V3, V5)) %>% dplyr::select(circRNA_id)

# Combine database
Circ_Database <- rbind(circAtlas, circBase) %>% distinct()

## CircRNA Annotation

Circ_SD_Ns$Annotation_status <- paste0(Circ_SD_Ns$CircRNA_ID, Circ_SD_Ns$strand) %in% Circ_Database$circRNA_id
Circ_SD_Ns$Annotation_status[(Circ_SD_Ns$Annotation_status == "FALSE")] <- "Not Annotated"
Circ_SD_Ns$Annotation_status[(Circ_SD_Ns$Annotation_status == "TRUE")] <- "Annotated"



##############
## Total
## Barplot
Annotation <- Circ_SD_Ns %>% 
  group_by(Timepoint, Annotation_status) %>%
  summarise(n = n()) %>%
  mutate(freq = (n / sum(n))) %>% 
  ggplot(aes(x = Timepoint, y = freq, fill = Annotation_status, label = paste0(n, "\n(", scales::percent(freq), ")"))) + 
  geom_col(position = 'dodge', width = .75) + 
  geom_text(position = position_dodge(width = .75),    # move to center of bars
            vjust = -0.2,    # nudge above top of bar
            size = 3.7,
            fontface = 2) + 
  scale_y_continuous(labels = scales::percent, expand = c(0,0), limits = c(0, 1.1)) +
  labs(x = "Timepoint", y = "Percentage Frequency", 
       title = "Annotation Status of Total CircRNAs Across Timepoints",
       fill = "Annotation status") +
  theme_classic() +
  theme(axis.text = element_text(size = 15), 
        axis.title = element_text(size = 12, face = "bold"),
        legend.title = element_text(size = 12, face="bold"),
        legend.text = element_text(size=11),
        plot.title=element_text(size = 15, face = "bold", hjust=0.5))

Annotation


## DE

Annotation_DE <- Circ_SD_Ns %>% filter(abs(logFC) >= 2, FDR <= 0.05) %>%
  group_by(Timepoint, Annotation_status) %>%
  summarise(n = n()) %>%
  mutate(freq = (n / sum(n))) %>% 
  ggplot(aes(x = Timepoint, y = freq, fill = Annotation_status, 
             label = paste0(n, "\n(", scales::percent(freq), ")"))) + 
  geom_col(position = 'dodge', width = .75) + 
  geom_text(position = position_dodge(width = .75),    # move to center of bars
            vjust = -0.5,    # nudge above top of bar
            size = 3) + 
  scale_y_continuous(labels = scales::percent, expand = c(0,0), 
                     limits = c(0, 1.05)) +
  labs(x = "Timepoint", y = "Percentage Frequency", 
       title = "Annotation Status of DE CircRNAs Across Timepoints",
       fill = "Annotation_status") +
  theme_classic() +
  theme(plot.title=element_text(hjust=0.5))


Annotation_DE


##################


################################################################################################################################

library(tidyverse)
library(tidyr)


circ_3hrs_cpm <- read.table("3hrs/circRNA_CPM.csv", row.names=1, sep=",", header = T)
colnames(circ_3hrs_cpm) <- c("T3H_0", "T3H_1", "C3H_0", "C3H_1" )

circ_6hrs_cpm <- read.table("6hrs/circRNA_CPM.csv", row.names=1, sep=",", header = T)
colnames(circ_6hrs_cpm) <- c("T6H_0", "T6H_1", "C6H_0", "C6H_1" )

circ_12hrs_cpm <- read.table("12hrs/circRNA_CPM.csv", row.names=1, sep=",", header = T)
colnames(circ_12hrs_cpm) <- c("T12H_0", "T12H_1", "C12H_0", "C12H_1" )

circ_18hrs_cpm <- read.table("18hrs/circRNA_CPM.csv", row.names=1, sep=",", header = T)
colnames(circ_18hrs_cpm) <- c("T18H_0", "T18H_1", "C18H_0", "C18H_1" )


cpm_list <- list(circ_3hrs_cpm, circ_6hrs_cpm, circ_12hrs_cpm, circ_18hrs_cpm)
CPM_merge_1 <- merge(circ_3hrs_cpm, circ_6hrs_cpm, by = "row.names", all = T)
CPM_merge_2 <- merge(CPM_merge_1, circ_12hrs_cpm, by.x = 'Row.names', by.y = 'row.names', all = T)
CPM_merge_3 <- merge(CPM_merge_2, circ_18hrs_cpm, by.x = 'Row.names', by.y = 'row.names', all = T)

CPM_long <- CPM_merge_3 %>% pivot_longer(cols = c("T3H_0", "T3H_1", "C3H_0", "C3H_1",
                                                  "T6H_0", "T6H_1", "C6H_0", "C6H_1",
                                                  "T12H_0", "T12H_1", "C12H_0", "C12H_1",
                                                  "T18H_0", "T18H_1", "C18H_0", "C18H_1"),
                                         names_to = 'Datasets',
                                         values_to = 'CPM')



## Boxplot
Boxplot_2 = ggplot(CPM_long , aes(x=Datasets, y= log10(CPM + 1), fill=Datasets)) +  
  labs(x = "Timepoints", y = "log10(CPM + 1)", title = "Human CircRNA Expression Across Timepoints") +
  theme(text = element_text(size=22),
        axis.text = element_text(angle = 90, vjust = 0.5, hjust=0.5, size = 18),
        axis.title.x = element_blank(),
        axis.title.y = element_text(size = 18, face = "bold"),
        plot.title = element_text(size = 20, hjust = 0.5, face="bold")) +
  geom_boxplot() 
Boxplot_2


## Violin Plot
Violin_1 = ggplot(CPM_long, aes(x=Datasets, y= log10(CPM), fill=Datasets)) +  
  labs(x = "Datasets", y = "CPM", title = "Human CPM Distribution") + geom_violin(trim=FALSE)
Violin_1



## Select FC data
logFC_circRNA<- Circ_SD_Ns %>% select(logFC, Timepoint) %>% as_tibble()

## melt data to factorise
logFC_circRNA$Timepoint <- factor(logFC_circRNA$Timepoint, levels = c('3hrs', '6hrs', '12hrs', '18hrs'))

## Boxplot
Boxplot_1 = ggplot(logFC_circRNA , aes(x=Timepoint, y= logFC, fill=Timepoint)) +  
  labs(x = "Timepoints", y = "log2(FC)", title = "Human CircRNA Expression Across Timepoints") +
  theme(text = element_text(size=22),
        axis.text = element_text(angle = 90, vjust = 0.5, hjust=0.5, size = 18),
        axis.title.x = element_blank(),
        axis.title.y = element_text(size = 18, face = "bold"),
        plot.title = element_text(size = 20, hjust = 0.5, face="bold")) +
  geom_boxplot() 
Boxplot_1



## Violin Plot
Violin_1 = ggplot(logFC_circRNA, aes(x=Timepoint, y= logFC, fill=Timepoint)) +  
  labs(x = "Timepoints", y = "log10(FC)", title = "Human logFC Distribution") + geom_violin(trim=FALSE)
Violin_1

############################################################################################################################

library(ggbreak)

### Total
## Circ type

Circ_type <- Circ_SD_Ns %>% 
  mutate(gene_type = replace_na(gene_type, "Intergenic")) %>%
  mutate(gene_type = strsplit(gene_type, ","), gene_name = strsplit(gene_name, ",") ) %>% 
  unnest_longer(c(gene_name, gene_type)) %>% filter(grepl("protein_coding|TR_C_gene|TR_J_gene|TEC|Intergenic", gene_type)) %>%
  group_by(CircRNA_ID, Timepoint, circ_type) %>% distinct(CircRNA_ID, Timepoint, gene_name, gene_type, circ_type) %>%
  summarise(gene_name = toString(gene_name), gene_type = toString(gene_type)) %>%
  ungroup() %>% select(circ_type, Timepoint) %>% as_tibble()


## Barplot
Circ_type_1 <- Circ_type %>%
  group_by(Timepoint, circ_type) %>%
  summarise(n = n()) %>%
  mutate(freq = round((n / sum(n)), 3)) %>% 
  ggplot(aes(x = Timepoint, y = freq, fill = circ_type, 
             label = paste0(n, "\n(", scales::percent(freq), ")"))) + 
  geom_col(width = 0.9,
           position=position_dodge(0.9)) + 
  geom_text(position = position_dodge(width = .9),    # move to center of bars
            vjust = -0.5,    # nudge above top of bar
            size = 3.7,
            fontface=2) + 
  scale_y_continuous(labels = scales::percent, 
                     expand = c(0,0), 
                     limits = c(0, 1.05)) +
  labs(x = "Timepoint", y = "Percentage Frequency", 
       title = "Genomic Origin of Total CircRNAs Across Timepoints",
       fill = "Circular type") +
  theme_classic() +
  theme(axis.text = element_text(size = 15), 
        axis.title = element_text(size = 15, face = "bold"),
        legend.title = element_text(size = 12, face="bold"),
        legend.text = element_text(size=11),
        plot.title=element_text(size = 15, face = "bold", hjust=0.5))



Circ_type_1

ggp <- (Annotation)/(Circ_type_1)
ggp

### DE

Circ_type <- Circ_SD_DE %>%   
  mutate(gene_type = replace_na(gene_type, "Intergenic")) %>%
  mutate(gene_type = strsplit(gene_type, ","), gene_name = strsplit(gene_name, ",") ) %>% 
  unnest_longer(c(gene_name, gene_type)) %>% filter(grepl("protein_coding|TR_C_gene|TR_J_gene|TEC|Intergenic", gene_type)) %>%
  group_by(CircRNA_ID, Timepoint, circ_type) %>% distinct(CircRNA_ID, Timepoint, gene_name, gene_type, circ_type) %>%
  summarise(gene_name = toString(gene_name), gene_type = toString(gene_type)) %>%
  ungroup() %>% select(circ_type, Timepoint) %>% as_tibble()

## Barplot
Circ_type <- Circ_type %>%
  group_by(Timepoint, circ_type) %>%
  summarise(n = n()) %>%
  mutate(freq = (n / sum(n))) %>% 
  ggplot(aes(x = Timepoint, y = freq, fill = circ_type, 
             label = paste0(n, "\n(", scales::percent(freq), ")"))) + 
  geom_col(position = 'dodge') + 
  geom_text(position = position_dodge(width = 10),    # move to center of bars
            vjust = -0.5,    # nudge above top of bar
            size = 6) + 
  scale_y_continuous(labels = scales::percent, 
                     expand = c(0,0), 
                     limits = c(0, 1.05)) +
  labs(x = "Timepoint", y = "Percentage Frequency", 
       title = "Genomic Origin of DE CircRNAs Across Timepoints",
       fill = "Circular type") +
  theme_classic() +
  theme(plot.title=element_text(hjust=0.5, face = "bold"))



Circ_type





###############################################################################################################
library(tidyr)

## number of genes Total

nog <- Circ_SD_Ns %>% mutate(gene_type = strsplit(gene_type, ","), gene_name = strsplit(gene_name, ",") ) %>% 
  unnest_longer(c(gene_name, gene_type)) %>% filter(grepl("protein_coding|TR_C_gene|TR_J_gene|TEC", gene_type)) %>%
  group_by(CircRNA_ID, Timepoint, circ_type) %>% distinct(CircRNA_ID, Timepoint, gene_name, gene_type, circ_type) %>%
  summarise(gene_name = toString(gene_name), gene_type = toString(gene_type)) %>%
  ungroup()

nog1 <- nog %>% dplyr::mutate(Gene_number = sapply(strsplit(gene_name, ","), length)) %>%
  as.data.frame() %>%
  dplyr::mutate(category = cut(Gene_number, breaks = c(0,1, 2, 3, Inf), 
                               labels=c("1", "2", "3", ">3"))) %>%
  group_by(Timepoint, category) %>%
  dplyr::summarise(n = n()) %>%
  dplyr::mutate(freq = (n / sum(n))) %>% 
  ggplot(aes(x = Timepoint, y = freq, fill = category, label = paste0(n, "\n(",scales::percent(freq %>% round(3)), ")"))) + 
  geom_col(position = 'dodge') + 
  geom_text(position = position_dodge(width = .9),    # move to center of bars
            vjust = -0.5,    # nudge above top of bar
            size = 3.7,
            fontface = 2) + 
  scale_y_continuous(labels = scales::percent, 
                     expand = c(0,0), 
                     limits = c(0, 1.1)) +
  labs(x = "Timepoint", y = "Percentage Frequency", 
       title = "Number of Parental Genes Across Timepoints \nAccording to Total CircRNA Type",
       fill = "Number of genes") +
  theme_classic() +
  theme(axis.text = element_text(size = 15), 
        axis.title = element_text(size = 15, face = "bold"),
        legend.title = element_text(size = 12, face="bold"),
        legend.text = element_text(size=11),
        plot.title=element_text(size = 15, face = "bold", hjust=0.5))

nog1



## number of genes DE

nog <- Circ_SD_DE %>% mutate(gene_type = strsplit(gene_type, ","), gene_name = strsplit(gene_name, ",") ) %>% 
  unnest_longer(c(gene_name, gene_type)) %>% filter(grepl("protein_coding|TR_C_gene|TR_J_gene|TEC", gene_type)) %>%
  group_by(CircRNA_ID, Timepoint, circ_type) %>% distinct(CircRNA_ID, Timepoint, gene_name, gene_type, circ_type) %>%
  summarise(gene_name = toString(gene_name), gene_type = toString(gene_type)) %>%
  ungroup()

nog1 <- nog %>% mutate(Gene_number = sapply(strsplit(gene_name, ","), length)) %>%
  as.data.frame() %>%
  mutate(category = cut(Gene_number, breaks = c(0,1, 2, 3, Inf), 
                        labels=c("1", "2", "3", ">3"))) %>%
  group_by(Timepoint, circ_type, category) %>%
  summarise(n = n()) %>%
  mutate(freq = (n / sum(n))) %>% 
  ggplot(aes(x = Timepoint, y = freq, fill = category, label = paste0(n, "\n(",scales::percent(freq %>% round(3)), ")"))) + 
  geom_col(position = 'dodge') + 
  geom_text(position = position_dodge(width = .9),    # move to center of bars
            vjust = -0.5,    # nudge above top of bar
            size = 3) + 
  scale_y_continuous(labels = scales::percent, 
                     expand = c(0,0), 
                     limits = c(0, 1.1)) +
  facet_wrap(~ circ_type) +
  labs(x = "Timepoint", y = "Percentage Frequency", 
       title = "Number of Parental Genes Across Timepoints \nAccording to DE CircRNA Type") +
  theme(axis.text = element_text(size = 15), 
        axis.title = element_text(size = 20, face = "bold"),
        legend.title = element_text(size = 15, face="bold"),
        legend.text = element_text(size=12),
        plot.title=element_text(size = 20, face = "bold", hjust=0.5))

nog1

#############################################################################################################

## Volcano plots
Circ_SD_Ns$diffexpressed <- "Not Significant"
Circ_SD_Ns$diffexpressed[Circ_SD_Ns$logFC > 2 & Circ_SD_Ns$FDR < 0.05] <- "Significant upregulated"
Circ_SD_Ns$diffexpressed[Circ_SD_Ns$logFC < -2 & Circ_SD_Ns$FDR < 0.05] <- "Significant downregulated"

VC <- list()

for (i in c(3, 6, 12)) {
  
  VC[[i]] <- ggplot(Circ_SD_Ns[(Circ_SD_Ns$Timepoint %in% paste0(i,"hrs")),], 
                    aes(x=logFC, y=-log10(FDR), col=diffexpressed), size = 2) + 
    geom_point() + scale_color_manual(values = c( "grey22", "red2", "green3")) + 
    geom_hline(yintercept = -log10(0.05),
               linetype = "dashed") + 
    geom_vline(xintercept = c(2, -2),
               linetype = "dashed",
               colour = c("red2", "green3")) +
    labs(x = "logFC", y = "-log10 FDR" ,colour = "Expression", title = paste0("CircRNA ", i, "hrs Volcano Plot")) +
    scale_alpha_manual(values=c(2,2,2)) +
    theme(axis.text = element_text(size = 20), 
          axis.title = element_text(size = 15, face = "bold"),
          plot.title = element_text(size = 15, face = "bold", hjust = 0.5), 
          legend.title = element_text(size=12, face="bold") )
  
}

VC[[12]] <- ggplot(Circ_SD_Ns[(Circ_SD_Ns$Timepoint %in% paste0("12hrs")),], 
                   aes(x=logFC, y=-log10(FDR), col=diffexpressed), size = 2) + 
  geom_point() + scale_color_manual(values = c( "grey22", "green3", "red2")) + 
  geom_hline(yintercept = -log10(0.05),
             linetype = "dashed") + 
  geom_vline(xintercept = c(2, -2),
             linetype = "dashed",
             colour = c("red2", "green3")) +
  labs(x = "logFC", y = "-log10 FDR" ,colour = "Expression", title = paste0("CircRNA 18hrs Volcano Plot")) +
  scale_alpha_manual(values=c(2,2,2)) +
  theme(axis.text = element_text(size = 20), 
        axis.title = element_text(size = 15, face = "bold"),
        plot.title = element_text(size = 15, face = "bold", hjust = 0.5), 
        legend.title = element_text(size=12, face = "bold") )




VC[[18]] <- ggplot(Circ_SD_Ns[(Circ_SD_Ns$Timepoint %in% paste0("18hrs")),], 
                   aes(x=logFC, y=-log10(FDR), col=diffexpressed), size = 2) + 
  geom_point() + scale_color_manual(values = c( "grey22", "green3", "red2")) + 
  geom_hline(yintercept = -log10(0.05),
             linetype = "dashed") + 
  geom_vline(xintercept = c(2, -2),
             linetype = "dashed",
             colour = c("red2", "green3")) +
  labs(x = "logFC", y = "-log10 FDR" ,colour = "Expression", title = paste0("CircRNA 18hrs Volcano Plot")) +
  scale_alpha_manual(values=c(2,2,2)) +
  theme(axis.text = element_text(size = 20), 
        axis.title = element_text(size = 15, face = "bold"),
        plot.title = element_text(size = 15, face = "bold", hjust = 0.5), 
        legend.title = element_text(size=12, face = "bold") )


library(patchwork)
ggp<- (VC[[3]] + VC[[6]])/(VC[[12]] + VC[[18]]) & theme(legend.position = "bottom")  
ggp + plot_layout(guides = "collect")

ggp<- (VC[[3]] + VC[[6]])/(VC[[12]] + VC[[18]])   
ggp



# ggplot(Circ_SD_Ns[(Circ_SD_Ns$Timepoint %in% "3hrs"),], aes(x=logFC, y=-log10(FDR), col=diffexpressed), size = 3) + 
#   geom_point() + scale_color_manual(values = c("grey22", "red2", "grey22", "grey22")) + 
#   geom_hline(yintercept = -log10(0.05),
#              linetype = "dashed") + 
#   geom_vline(xintercept = c(2, -2),
#              linetype = "dashed",
#              colour = c("red2", "green3")) +
#   labs(x = "logFC", y = "-log10 FDR" ,colour = "Expression", title =" CircRNA 3hrs Volcano Plot")

##############################################################################################################

## Heatmap

library(tidyverse)
library(pheatmap)
library(gplots)
library(dendextend)
library(colorspace)
library(ggplot2)
library(readxl)
library(stringr)


setwd('C:/Users/matha/Documents/CircRNA_Analysis_Compile/CIRI2_CIRIquant_Ind_NoStringent')

circ_3hrs_cpm <- read.table("3hrs/circRNA_CPM.csv", row.names=1, sep=",", header = T)
circ_6hrs_cpm <- read.table("6hrs/circRNA_CPM.csv", row.names=1, sep=",", header = T)
circ_12hrs_cpm <- read.table("12hrs/circRNA_CPM.csv", row.names=1, sep=",", header = T)
circ_18hrs_cpm <- read.table("18hrs/circRNA_CPM.csv", row.names=1, sep=",", header = T)

Circ_3hrs_merge <- Circ_3hrs %>% merge(.,circ_3hrs_cpm, by.x= "X", by.y = 0 , all=T) %>% select(c(1,9:12))
Circ_6hrs_merge <- Circ_6hrs %>% merge(.,circ_6hrs_cpm, by.x= "X", by.y = 0 , all=T) %>% select(c(1,9:12))
Circ_12hrs_merge <- Circ_12hrs %>% merge(.,circ_12hrs_cpm, by.x= "X", by.y = 0 , all=T) %>% select(c(1,9:12))
Circ_18hrs_merge <- Circ_18hrs %>% merge(.,circ_18hrs_cpm, by.x= "X", by.y = 0 , all=T) %>% select(c(1,9:12))


CPM_merge <- merge(merge(merge(Circ_3hrs_merge, Circ_6hrs_merge, by="X", all=T), 
                         Circ_12hrs_merge,by="X", all=T), Circ_18hrs_merge,by="X", all=T)

CPM_merge <- CPM_merge %>% 
  replace(is.na(.), 0) %>% 
  column_to_rownames("X")

colnames(CPM_merge) <- c("T3H_0", "T3H_1", "C3H_0", "C3H_1", "T6H_0", "T6H_1", "C6H_0", "C6H_1", 
                         "T12H_0", "T12H_1", "C12H_0", "C12H_1", "T18H_0", "T18H_1", "C18H_0", "C18H_1")

logCPM_merge <- CPM_merge%>%
  mutate(log10(across(where(is.numeric), ~ .x + 1)))

logCPM_merged<- array(0, dim=dim(logCPM_merge))

for (i in 1:ncol(logCPM_merge)) {logCPM_merged[,i] <- logCPM_merge[,i]}
logCPM_merged[, 1:16] <- sapply(logCPM_merged[, 1:16], as.numeric)

colnames(logCPM_merged) <-colnames(logCPM_merge)
rownames(logCPM_merged) <- rownames(logCPM_merge)


logCPM_merged <- as.data.frame(logCPM_merged)

logCPM_merged <- logCPM_merged[rowSums(logCPM_merged[])>0,] %>% as.matrix()


# Arrange col order
col.order <- c("C3H_0", "C3H_1", "C6H_0", "C6H_1", "C12H_0", "C12H_1", "C18H_0", "C18H_1", 
               "T3H_0", "T3H_1", "T6H_0", "T6H_1", "T12H_0", "T12H_1", "T18H_0", "T18H_1")
logCPM_merged <- logCPM_merged[, col.order]

#  Define colours
hmcols <- rev(redgreen(2750));


# plot heatmap
my_heatmap <- pheatmap(logCPM_merged,
                       cluster_rows = TRUE,
                       cluster_cols = FALSE, clustering_distance_rows = "euclidean",
                       clustering_distance_cols = "euclidean", clustering_method = "complete",
                       cutree_rows = 1, cutree_cols = 1, fontsize = 12, 
                       fontsize_row = 8, color=hmcols,scale="row", show_rownames = F,
                       main = "Human H1N1 CircRNA Heatmap")




## Plot DE heatmap


circDE_list <- Circ_SD_DE$CircRNA_ID %>% unique()

logCPM_DE_merged <- logCPM_merged[rownames(logCPM_merged) %in% circDE_list,] 


my_heatmap <- pheatmap(logCPM_DE_merged,
                       cluster_rows = TRUE,
                       cluster_cols = FALSE, clustering_distance_rows = "euclidean",
                       clustering_distance_cols = "euclidean", clustering_method = "complete",
                       cutree_rows = 1, cutree_cols = 1, fontsize = 12, 
                       fontsize_row = 8, color=hmcols,scale="row", show_rownames = F,
                       main = "Human H1N1 DE CircRNA Heatmap")



# common_DE <- generics::intersect(Circ_SD_DE$CircRNA_ID[(Circ_SD_DE$Timepoint %in% "12hrs")], 
#                                  Circ_SD_DE$CircRNA_ID[(Circ_SD_DE$Timepoint %in% "18hrs")]
#                                  )
# 
# common_DE <- Reduce(intersect, list(Circ_SD_Ns$CircRNA_ID[(Circ_SD_Ns$Timepoint %in% "12hrs")], 
#                                     Circ_SD_Ns$CircRNA_ID[(Circ_SD_Ns$Timepoint %in% "18hrs")],
#                                     Circ_SD_Ns$CircRNA_ID[(Circ_SD_Ns$Timepoint %in% "3hrs")],
#                                     Circ_SD_Ns$CircRNA_ID[(Circ_SD_Ns$Timepoint %in% "6hrs")]))
#              
# 
# Circ_heatmap <- Circ_SD_DE[(Circ_SD_DE$CircRNA_ID %in% common_DE),]
# 
# Circ_heatmap <- Circ_SD_Ns %>% distinct(CircRNA_ID, Timepoint, gene_name, .keep_all = T) %>% 
#   dplyr::select(logFC, Timepoint) %>% group_by(Timepoint) %>%
#   mutate(row = row_number()) %>%
#   # mutate(row = row_number()) %>%
#   # mutate(logCPM = exp(logCPM)) %>%
#   # mutate(logCPM = log10(logCPM + 1)) %>%
#   tidyr::pivot_wider(names_from = Timepoint, values_from = logFC) %>%
#   dplyr::select(-row)
# 
# 
# 
# Circ_heatmap <- Circ_heatmap %>% replace(is.na(.), 0)
# 
# #  Define colours
# hmcols <- rev(redgreen(2750));
# 
# # plot heatmap
# my_heatmap <- pheatmap(Circ_heatmap,
#                        cluster_rows = T,
#                        cluster_cols = F, clustering_distance_rows = "euclidean",
#                        clustering_distance_cols = "euclidean", clustering_method = "complete",
#                        cutree_rows = 1, cutree_cols = 1, fontsize = 12, 
#                        fontsize_row = 8, color=hmcols,scale="row", show_rownames = F, na_col = "grey",
#                        main = "Human DE Circular Heatmap")

###############################################################################################################

## Strand type

Circ_strand <- Circ_SD_Ns %>% distinct(CircRNA_ID, logFC, FDR, Timepoint, .keep_all = T) %>% 
  select(CircRNA_ID, strand, Timepoint) %>% as_tibble()

Circ_strand_1 <- Circ_strand %>%
  separate(CircRNA_ID, 
           into = c("Chr", "Position"), 
           sep = ":"
  )  %>% filter(!(grepl("GL000219.1|KI270711.1", Chr))) %>%
  group_by(Chr, Timepoint, strand) %>%
  summarise(n = n()) %>%
  mutate(freq = (n / sum(n)))


Level_order <- factor(Circ_strand$Chr, levels = c('chr1', 'chr2', 'chr3', 'chr4',
                                                  'chr5', 'chr6', 'chr7', 'chr8', 'chr9', 'chr10',
                                                  'chr11', 'chr12', 'chr13', 'chr14', 'chr15', 
                                                  'chr16', 'chr17', 'chr18', 'chr19', 'chr20', 
                                                  'chr21', 'chr22', 'chrX', 'chrY'))

Circ_strand <- Circ_strand %>%
  ggplot(aes(x = Level_order, y = n, fill = strand, label = paste0(n, "\n (", scales::percent(round(freq, 3)), ")"))) + 
  geom_col(position = 'stack') + 
  geom_text(position = position_stack(vjust = 0.5),    # move to center of bars
            size = 3) + 
  labs(x = "Chromosomal Location", y = "Frequency", 
       title = "CircRNA Strand Type Across Timepoints") +
  facet_wrap_paginate(~Timepoint, ncol = 1, nrow = 1, page = 4, scales = "free")

Circ_strand 


##################################################################################

## Chr Distribution

chr_length <- read.table("chromosome_length.txt", sep = "\t", header = T)
chr_length$Chr <- paste0("chr", chr_length$Chr)

Circ_chr <- Circ_SD_Ns %>% distinct(CircRNA_ID, .keep_all = T) %>% 
  select(CircRNA_ID, strand, Timepoint) %>% as_tibble()


Circ_chr <- Circ_strand %>%
  separate(CircRNA_ID, 
           into = c("Chr", "Position"), 
           sep = ":"
  )  %>% filter(!(grepl("GL000219.1|KI270711.1", Chr))) %>%
  group_by(Chr, Timepoint) %>%
  summarise(n = n())

Circ_chr <- merge(Circ_chr, chr_length, by = "Chr")

Circ_chr <- Circ_chr %>% mutate(Normalized_counts = as.numeric(n)/(as.numeric(gsub(",", "",Length.bp.))/1000))



ggp <- ggplot(Circ_chr %>% filter(Timepoint == "18hrs") %>% select(c(1, 3, 5)))  + 
  geom_bar(aes(x=Chr, y=n),stat="identity", fill="cyan",colour="#006000")+
  geom_line(aes(x=Chr, y=100000*Normalized_counts),stat="identity",color="red",linewidth=1, group = 1)+
  labs(title= "Distribution of circRNAs per chromosome",
       x="Chr number",y="Number of circRNAs")+
  scale_y_continuous(sec.axis=sec_axis(~.*0.00001,name="circRNA/kb"))
ggp






Level_order <- factor(Circ_strand$Chr, levels = c('chr1', 'chr2', 'chr3', 'chr4',
                                                  'chr5', 'chr6', 'chr7', 'chr8', 'chr9', 'chr10',
                                                  'chr11', 'chr12', 'chr13', 'chr14', 'chr15', 
                                                  'chr16', 'chr17', 'chr18', 'chr19', 'chr20', 
                                                  'chr21', 'chr22', 'chrX', 'chrY'))

Circ_strand <- Circ_strand %>%
  ggplot(aes(x = Level_order, y = n, fill = strand, label = paste0(n, "\n (", scales::percent(round(freq, 3)), ")"))) + 
  geom_col(position = 'stack') + 
  geom_text(position = position_stack(vjust = 0.5),    # move to center of bars
            size = 3) + 
  labs(x = "Chromosomal Location", y = "Frequency", 
       title = "CircRNA Strand Type Across Timepoints") +
  facet_wrap_paginate(~Timepoint, ncol = 1, nrow = 1, page = 4, scales = "free")

Circ_strand 


##################################################################################

library(ggVennDiagram)

Circ_DE_12hrs <- Circ_SD_Ns %>% filter(Timepoint == "12hrs", round(abs(logFC),0) >= 2, round(FDR,2) <= 0.05)

Circ_DE_18hrs <- Circ_SD_Ns %>% filter(Timepoint == "18hrs", round(abs(logFC),0) >= 2, round(FDR,2) <= 0.05)

S = list(Circ_DE_12hrs$CircRNA_ID, Circ_DE_18hrs$CircRNA_ID)

ggVennDiagram(S, category.names = c("DE 12hrs", "DE 18hrs"))
venn <- Venn(S)
data <- process_data(venn)
L = as.data.frame(venn_region(data))
L_1 = as.data.frame(L[,5])
L_2 = L_1[2,]/(L_1[1,] + L_1[3,])
L_3 = as.data.frame(paste0(t(L_1), "\n  ", L_2, "","%"))
colnames(L_3) <- "Both"
venn_region(data) %>%
  mutate(Algorithms = c("12 hrs", "18 hrs", "Common")) %>% 
  mutate(count = paste0(count, "\n  ", round(prop.table(count)*100),"%")) -> R
R %>% ggplot() +
  geom_sf(aes(fill=Algorithms), colour = alpha(c("dodgerblue4", "#B30000", "#006D2C"), 0.8), show.legend = F) +
  #geom_sf(aes(fill=Algorithms, colour = Algorithms), show.legend = T) +
  #geom_sf(aes(color=id), data = venn_setedge(data), show.legend = F) +
  geom_sf_text(aes(label = c("12 hrs", "18 hrs")), data = venn_setlabel(data), nudge_y = -10, nudge_x = 2, size = 6, fontface= "bold") +
  geom_sf_label(aes(label=L_1$`L[, 5]`), fontface = "bold", data = venn_region(data), size = 7, colour = "black", label.padding = unit(0.55, "lines")) +
  scale_fill_manual(values = alpha(c("dodgerblue4", "#B30000", "#00441B"), 0.58)) + labs(title = " Overlap of DE circRNAs between \n12 hrs and 18 hrs") +
  theme_void() + theme(plot.title=element_text(family='', face='bold', size=16, hjust=0.5, vjust=1))







# Writing common circRNAs from 12hrs and 18hrs
common <- Circ_DE_12hrs[Circ_DE_12hrs$CircRNA_ID %in% Circ_DE_18hrs$CircRNA_ID,] %>% dplyr::select(CircRNA_ID, strand, circ_type, gene_id, gene_name, Annotation_status,
                                                                                                   gene_type, logFC)

circ_type <- common %>% 
  dplyr::select(circ_type) %>% 
  group_by(circ_type) %>% 
  summarise(n = n()) %>%
  mutate(freq = (n / sum(n)))

mycols <- c("#0073C2FF", "#CD534CFF", "#868686FF", "#EFC000FF")

pie.chart <- ggplot(circ_type, aes(x = "", y = n, fill = circ_type)) +
  geom_bar(width = 1, stat = "identity", color = "white") +
  coord_polar("y", start = 0)+
  geom_text(aes(label = paste0(n, "\n", scales::percent(freq))), 
            color = "white", position = position_stack(vjust = 0.5))+
  scale_fill_manual(values = mycols) +
  theme_void()


pie.chart 


# Extract exons and remove lincRNA
common_exon <- common %>% filter(circ_type == "exon", !grepl("lincRNA", gene_type))

## Remove those gene within genes that overlap with each other
common_exon_pc <- common_exon %>%
  dplyr::select(-gene_id) %>%
  mutate(gene_type = strsplit(gene_type, ","), gene_name = strsplit(gene_name, ",") ) %>%
  unnest(c(gene_name, gene_type)) %>%
  filter(gene_type == "protein_coding", !grepl("TAF9", gene_name)) %>%
  group_by(CircRNA_ID) %>%
  mutate(gene_name = str_c(gene_name, collapse=","), gene_type = str_c(gene_type, collapse=",")) %>%
  distinct(CircRNA_ID, .keep_all = T) %>%
  mutate(Gene_number = sapply(strsplit(gene_name, ","), length)) %>%
  filter(!Gene_number > 1) %>% as.data.frame()

write.xlsx2(common_exon_pc, "C:/Users/matha/Documents/CircRNA_Analysis_Compile/CIRI2_CIRIquant_Ind_NoStringent/Common_DE_12&18hrs.xlsx", col.names = T,
            row.names = F, sep="\t", sheetName ="Common_DE_12hrs&18hrs")




##################################################################################

### Number of circRNAs

CircRNA_DE_table <- Circ_SD_DE %>% group_by(Timepoint) %>% summarise(n = n())

CircRNA_table <- Circ_SD_Ns %>% group_by(Timepoint) %>%
  summarise(n = n())

#################################################################################
library(clusterProfiler)
library(enrichplot)
library(org.Hs.eg.db)
library(tidyr)

##### GO and KEGG

### 12hrs 

Circ_SD_DE_12hrs <- Circ_SD_DE %>% mutate(gene_type = strsplit(gene_type, ","), gene_name = strsplit(gene_name, ",") ) %>% 
  unnest(c(gene_type, gene_name)) %>% 
  filter(Timepoint == "12hrs", circ_type == "exon", grepl("protein_coding|TR_C_gene|TR_J_gene|TEC", gene_type))


Circ_SD_DE_12hrs <- clusterProfiler::bitr(Circ_SD_DE_12hrs$gene_name, fromType = "SYMBOL", toType = "ENTREZID", OrgDb = org.Hs.eg.db)


# Enrich GO terms
GOenrich <- enrichGO(Circ_SD_DE_12hrs$ENTREZID, keyType = "ENTREZID", org.Hs.eg.db, ont = "All", 
                     pvalueCutoff = 0.05, pAdjustMethod = "BH", qvalueCutoff = 0.05)


dotplot(GOenrich,showCategory = 10, size = NULL,font.size = 8, 
        title = "Circular 12hrs GO  Enrichment", split = "ONTOLOGY") + facet_grid_sc(ONTOLOGY ~ ., scales = "free")



# Enrich KEGG terms
keg_res <- enrichKEGG(Circ_SD_DE_12hrs$ENTREZID, organism="hsa", keyType="kegg", pvalueCutoff = 0.05, 
                      pAdjustMethod = "BH", qvalueCutoff = 0.05)

# Dotplot for KEGG
dotplot(keg_res,showCategory = 10, size = NULL,font.size = 8, title = "Circular 12hrs KEGG enrichment") 

######################################################################################################

### 18hrs 

Circ_SD_DE_18hrs <- Circ_SD_DE %>% mutate(gene_type = strsplit(gene_type, ","), gene_name = strsplit(gene_name, ",") ) %>% 
  unnest(c(gene_type, gene_name)) %>% 
  filter(Timepoint == "18hrs", circ_type == "exon", grepl("protein_coding|TR_C_gene|TR_J_gene|TEC", gene_type))


Circ_SD_DE_18hrs <- clusterProfiler::bitr(Circ_SD_DE_18hrs$gene_name, fromType = "SYMBOL", toType = "ENTREZID", OrgDb = org.Hs.eg.db)


# Enrich GO terms
GOenrich <- enrichGO(Circ_SD_DE_18hrs$ENTREZID, keyType = "ENTREZID", org.Hs.eg.db, ont = "All", 
                     pvalueCutoff = 0.05, pAdjustMethod = "BH", qvalueCutoff = 0.05)


dotplot(GOenrich,showCategory = 10, size = NULL,font.size = 8, 
        title = "Circular 18hrs GO  Enrichment", split = "ONTOLOGY") + facet_grid(ONTOLOGY ~ ., scales = "free")



# Enrich KEGG terms
keg_res <- enrichKEGG(Circ_SD_DE_18hrs$ENTREZID, organism="hsa", keyType="kegg", pvalueCutoff = 0.05, 
                      pAdjustMethod = "BH", qvalueCutoff = 0.05)

# Dotplot for KEGG
dotplot(keg_res,showCategory = 10, size = NULL,font.size = 8, title = "Circular 18hrs KEGG enrichment") 


############################################################################################################

### cOMMON 12hrs and 18hrs

Circ_SD_DE_12_18hrs <- Circ_SD_DE_12hrs[Circ_SD_DE_12hrs$SYMBOL %in% Circ_SD_DE_18hrs$SYMBOL,]


# Enrich GO terms
GOenrich <- enrichGO(Circ_SD_DE_12_18hrs$ENTREZID, keyType = "ENTREZID", org.Hs.eg.db, ont = "All", 
                     pvalueCutoff = 0.05, pAdjustMethod = "BH", qvalueCutoff = 0.05)


dotplot(GOenrich,showCategory = 10, size = NULL,font.size = 8, 
        title = "Circular 12hrs and 18hrs GO  Enrichment", split = "ONTOLOGY") + facet_grid_sc(ONTOLOGY ~ ., scales = "free")



# Enrich KEGG terms
keg_res <- enrichKEGG(Circ_SD_DE_12_18hrs$ENTREZID, organism="hsa", keyType="kegg", pvalueCutoff = 0.05, 
                      pAdjustMethod = "BH", qvalueCutoff = 0.05)

# Dotplot for KEGG
dotplot(keg_res,showCategory = 10, size = NULL,font.size = 8, title = "Circular 12hrs and 18hrs KEGG enrichment") 


############################################################################################################

### CeRNA 12hrs and 18hrs

library(readxl)
library(pathview)

circRNA_ceRNA <- paste0("hsa_DE_", c(3, 4, 6, 11, 13, 14, 15, 20, 24, 31, 33, 34, 38, 39, 42, 43, 45, 47, 48, 50, 
                                     51, 53, 58, 59, 60, 62))

ceRNA_enrichment <- read_excel("C:/Users/matha/Documents/ceRNA_GO_KEGG/CircRNAs_overlap12_18.xlsx", sheet = 1)

ceRNA_enrichment_circs <- ceRNA_enrichment[ceRNA_enrichment$customName %in% circRNA_ceRNA,] %>% dplyr::select(gene_name) %>% 
  mutate(gene_name = strsplit(gene_name, ",") ) %>% 
  unnest(gene_name)

ceRNA_enrichment_circs <- clusterProfiler::bitr(ceRNA_enrichment_circs$gene_name, fromType = "SYMBOL", toType = "ENTREZID", OrgDb = org.Hs.eg.db)


# Enrich GO terms
GOenrich <- enrichGO(ceRNA_enrichment_circs$ENTREZID, keyType = "ENTREZID", org.Hs.eg.db, ont = "All", 
                     pvalueCutoff = 0.01, pAdjustMethod = "BH", qvalueCutoff = 0.05)


dotplot(GOenrich,showCategory = 10, size = NULL,font.size = 8, 
        title = "CeRNA Circular 12hrs and 18hrs GO  Enrichment", split = "ONTOLOGY") + facet_grid_sc(ONTOLOGY ~ ., scales = "free")



# Enrich KEGG terms
keg_res <- enrichKEGG(ceRNA_enrichment_circs$ENTREZID, organism="hsa", keyType="kegg", pvalueCutoff = 0.05, 
                      pAdjustMethod = "BH", qvalueCutoff = 0.05)
keg_res_symbol <- setReadable(keg_res, OrgDb = org.Hs.eg.db, keyType="ENTREZID")

# Pathway view
hsa05164 <- pathview(gene.data  = ceRNA_enrichment_circs$ENTREZID,
                     pathway.id = "hsa05164",
                     species    = "hsa")



write.xlsx2(keg_res_symbol@result, file = "C:/Users/matha/Desktop/Master/Progress/7th Nov 2022/H1N1/kegg_ceRNA_DEcircRNA.xlsx", 
            col.names = T)

# Dotplot for KEGG
dotplot(keg_res,showCategory = 10, size = NULL,font.size = 8, title = "CeRNA Circular 12hrs and 18hrs KEGG enrichment") 

##################################################################################################################################


ggp<- (Annotation) /(Circ_type_1)   
ggp


