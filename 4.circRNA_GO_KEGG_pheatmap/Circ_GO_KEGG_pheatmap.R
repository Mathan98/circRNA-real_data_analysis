## GO KEGG pheatmap CircRNA


## Circular
library(ggplot2)
library(dplyr)
library(rtracklayer)
library(pheatmap)
library(org.Hs.eg.db)
library(stringr)
library(tidyr)
library(clusterProfiler)
library(tibble)
library(RColorBrewer)
library(org.Hs.eg.db)

## Read and combine data from different timepoints

# setwd('C:/Users/matha/Documents/CircRNA_Analysis_Compile/CIRIquant_Superdepth_noStringency/Allhrs_Superdepth_NoStringency/')
setwd('D:/circRNA/CIRIquant_Superdepth_noStringency/')


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
  GTF[[i]] <- rtracklayer::import(paste0( i, ".gtf", sep="")) %>% as.data.frame() %>%
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



#####################################################################################################################################################

# 6hrs
Circ_SD_DE_6hrs <- Circ_SD_DE %>% 
                          mutate(gene_type = strsplit(gene_type, ","), gene_name = strsplit(gene_name, ",") ) %>% 
                          unnest(c(gene_type, gene_name)) %>% 
                          filter(Timepoint == "6hrs", grepl("protein_coding|TR_C_gene|TR_J_gene|TEC", gene_type))
Circ_SD_DE_6hrs$gene_id <- sub("\\..*", "",Circ_SD_DE_6hrs$gene_id)
Circ_SD_DE_6hrs <- clusterProfiler::bitr(Circ_SD_DE_6hrs$gene_name, fromType = "ENSEMBL", toType = c("SYMBOL","ENTREZID"), OrgDb = org.Hs.eg.db)



# 12hrs 
Circ_SD_DE_12hrs <- Circ_SD_DE %>% 
                          mutate(gene_type = strsplit(gene_type, ","), gene_name = strsplit(gene_name, ",") ) %>% 
                          unnest(c(gene_type, gene_name)) %>% 
                          filter(Timepoint == "12hrs",  grepl("protein_coding|TR_C_gene|TR_J_gene|TEC", gene_type))
Circ_SD_DE_12hrs$gene_id <- sub("\\..*", "",Circ_SD_DE_12hrs$gene_id)
Circ_SD_DE_12hrs <- clusterProfiler::bitr(Circ_SD_DE_12hrs$gene_id, fromType = "ENSEMBL", toType = "ENTREZID", OrgDb = org.Hs.eg.db)


# 18hrs 
Circ_SD_DE_18hrs <- Circ_SD_DE %>% 
                          mutate(gene_type = strsplit(gene_type, ","), gene_name = strsplit(gene_name, ",") ) %>% 
                          unnest(c(gene_type, gene_name)) %>% 
                          filter(Timepoint == "18hrs",  grepl("protein_coding|TR_C_gene|TR_J_gene|TEC", gene_type))
Circ_SD_DE_18hrs$gene_id <- sub("\\..*", "",Circ_SD_DE_18hrs$gene_id)
Circ_SD_DE_18hrs <- clusterProfiler::bitr(Circ_SD_DE_18hrs$gene_id, fromType = "ENSEMBL", toType = "ENTREZID", OrgDb = org.Hs.eg.db)



#####################################################################################################################################################

# Enrich GO terms

# 12hrs
go_res_12 <- enrichGO(Circ_SD_DE_12hrs$ENTREZID, keyType = "ENTREZID", org.Hs.eg.db, ont = "All", 
                     pvalueCutoff = 0.05, pAdjustMethod = "BH", qvalueCutoff = 0.05, readable = T)

# 18hrs
go_res_18 <- enrichGO(Circ_SD_DE_18hrs$ENTREZID, keyType = "ENTREZID", org.Hs.eg.db, ont = "All", 
                     pvalueCutoff = 0.05, pAdjustMethod = "BH", qvalueCutoff = 0.05, readable = T)

#####################################################################################################################################################

# Bind the GO results from different timepoints
GO_rbind <- rbind(go_res_12@result %>% mutate(Timepoint = "12hrs"),
                  go_res_18@result %>% mutate(Timepoint = "18hrs")) %>%
  dplyr::select(c("ONTOLOGY", "Description", "p.adjust", "Timepoint"))


# Change from long to wide
GO_wide <- GO_rbind %>%
  group_by(Description) %>%
  pivot_wider(names_from = Timepoint, values_from = p.adjust) %>%
  column_to_rownames("Description")


# Split according to Ontology
GO_wide <- split(GO_wide, GO_wide$ONTOLOGY)


# Get All GO terms
GO_wide_ALL <- GO_rbind %>%
  group_by(Description) %>%
  pivot_wider(names_from = Timepoint, values_from = p.adjust) %>%
  column_to_rownames("Description") %>% 
  dplyr::select(c(2:3)) %>% 
  as.matrix()


# Make a matrix out of the split Ontology
GO_BP <- GO_wide[["BP"]] %>% dplyr::select(c(2:3)) %>% as.matrix()
GO_MF <- GO_wide[["MF"]] %>% dplyr::select(c(2:3)) %>% as.matrix()
GO_CC <- GO_wide[["CC"]] %>% dplyr::select(c(2:3)) %>% as.matrix() ## NO CC for circRNA

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


# GO_BP
# GO_10 <- GO_rbind %>% 
#   arrange(p.adjust) %>% 
#   group_by(ONTOLOGY) %>%
#   filter(ONTOLOGY == "BP") %>%
#   distinct(Description) %>%
#   slice(1:10) %>%
#   as.list()

GO_BP_1 <- GO_BP %>% as.data.frame()
GO_BP_1$`3hrs` <- NA
GO_BP_1$`6hrs` <- NA

GO_BP <- GO_BP_1 %>%  dplyr::select(3, 4, 1, 2)
GO_BP_1 <- GO_BP %>% slice(1:7)
GO_BP_2 <- GO_BP %>% slice(9:11)
GO_BP <- rbind(GO_BP_1, GO_BP_2)

# GO_BP_10 <- GO_BP[(rownames(GO_BP) %in% GO_10[["Description"]]),]

GO_BP_HM <- pheatmap(GO_BP, cellwidth = 40, cellheight = 40, 
                     cluster_cols = F, cluster_rows = F,
                     #treeheight_row = 40, treeheight_col = 40,
                     clustering_distance_rows = "euclidean",
                     clustering_distance_cols = "euclidean", clustering_method = "complete",
                     fontsize = 14,
                     fontsize_col = 20, fontsize_row = 18,
                     border_color='grey60', show_rownames = T,
                     labels_row = make_bold_names(GO_BP, rownames, c('protein polyubiquitination', 'defense response to virus', 'proteasome-mediated ubiquitin-dependent protein catabolic process')),
                     labels_col = make_bold_names(GO_BP, colnames, colnames(GO_BP)),
         main = "Top 10 GO BP Enriched Terms")
         #filename = "GO_BP_10_pheatmap.png")
         

# GO_MF
# GO_10 <- GO_rbind %>% 
#   arrange(p.adjust) %>% 
#   group_by(ONTOLOGY) %>%
#   filter(ONTOLOGY == "MF") %>%
#   distinct(Description) %>%
#   slice(1:10) %>%
#   as.list()


GO_MF_1 <- GO_MF %>% as.data.frame()
GO_MF_1$`3hrs` <- NA
GO_MF_1$`6hrs` <- NA

GO_MF <- GO_MF_1 %>%  dplyr::select(3, 4, 1, 2)
GO_MF_1 <- GO_MF %>% slice(1:6)
GO_MF_2 <- GO_MF %>% slice(8:11)
GO_MF <- rbind(GO_MF_1, GO_MF_2)

# GO_MF_10 <- GO_MF[(rownames(GO_MF) %in% GO_10[["Description"]]),]


GO_MF_HM <-pheatmap(GO_MF, cellwidth = 40, cellheight = 40, 
                    cluster_cols = F, cluster_rows = F,
         #treeheight_row = 40, treeheight_col = 40,
         clustering_distance_rows = "euclidean",
         clustering_distance_cols = "euclidean", clustering_method = "complete",
         fontsize = 14,
         fontsize_col = 20, fontsize_row = 18,
         border_color='grey60', show_rownames = T,
         labels_row = make_bold_names(GO_MF, rownames, c('ubiquitin-like protein ligase activity', 'aminoacyltransferase activity', 'ubiquitin-like protein transferase activity')),
         labels_col = make_bold_names(GO_MF, colnames, colnames(GO_MF)),
         main = "Top 10 GO MF Enriched Terms")
         #filename = "GO_MF_10_pheatmap.png")


## Save plots
library(ggpubr)
library(ggplotify)
p1 <- as.ggplot(GO_BP_HM)
p2 <- as.ggplot(GO_MF_HM)


# use patchwork to arrange them together
figure <- ggarrange(p1, p2,
                    #labels = c("A", "B"),
                    ncol = 2, nrow = 1)
figure

#####################################

## CC Ontology

####################################
GO_CC_1 <- GO_CC %>% as.data.frame()
GO_CC_1$`3hrs` <- NA
GO_CC_1$`6hrs` <- NA

GO_CC <- GO_CC_1 %>%  dplyr::select(3, 4, 1, 2)
GO_CC_1 <- GO_CC %>% slice(1:8)
GO_CC_2 <- GO_CC %>% slice(10:11)
GO_CC <- rbind(GO_CC_1, GO_CC_2)

# GO_BP_10 <- GO_BP[(rownames(GO_BP) %in% GO_10[["Description"]]),]

GO_CC_HM <- pheatmap(GO_CC, cellwidth = 40, cellheight = 40, 
                     cluster_cols = F, cluster_rows = F,
                     #treeheight_row = 40, treeheight_col = 40,
                     clustering_distance_rows = "euclidean",
                     clustering_distance_cols = "euclidean", clustering_method = "complete",
                     fontsize = 14,
                     fontsize_col = 20, fontsize_row = 18,
                     border_color='grey60', show_rownames = T,
                     labels_row = make_bold_names(GO_CC, rownames, c('cytoplasmic stress granule', 'lysosomal membrane', 'early endosome')),
                     labels_col = make_bold_names(GO_CC, colnames, colnames(GO_CC)),
                     main = "Top 10 GO CC Enriched Terms")



#######################################################################################################################


## Create a long df
logFC_circ <- Circ_SD_Ns %>% 
                  dplyr::select(CircRNA_ID, logFC, Timepoint, gene_name)
                  # transform(gene_name = strsplit(gene_name, ",")) %>% 
                  # unnest(gene_name)

LogFC_gene <- Circ_SD_Ns %>% select(CircRNA_ID, gene_name) %>% unique()
logFC_circ <- Circ_SD_Ns %>% select(CircRNA_ID, logFC, Timepoint)



## Convert to wider df
LogFC_Circmx <- logFC_circ %>% 
                    pivot_wider(names_from = Timepoint, values_from = logFC) %>%
                    as.data.frame() %>%
                    mutate_all( ~replace(., lengths(.)==0, 0))
LogFC_Circmx2 <- merge(LogFC_Circmx, LogFC_gene, by="CircRNA_ID", all.x = T)
LogFC_Circmx <- LogFC_Circmx2 %>% na.omit()
LogFC_Circmx$circ_name <- paste0("circ", ave(LogFC_Circmx$gene_name, LogFC_Circmx$gene_name, FUN = function(x) paste0(x, "_", seq_along(x))))



#### Extract DE genes from matrix
## DE gene list
DE <- Circ_SD_DE %>% 
  dplyr::select(CircRNA_ID)



DE_geneLogFC_Circmx <- LogFC_Circmx %>%
                              filter(CircRNA_ID %in% DE$CircRNA_ID) %>%
  transform(gene_name = strsplit(gene_name, ",")) %>% 
  unnest(gene_name)



####################################################################################################################
############## GO_BP

# Get the gene list and the counts for each timepoint
GO_geneList <- rbind(go_res_12@result %>% mutate(Timepoint = "12hrs"),
                     go_res_18@result %>% mutate(Timepoint = "18hrs")) %>%
  dplyr::select(c("ONTOLOGY", "Description", "geneID", "Count","Timepoint"))


# Get the top3 terms for BP
Terms <- c('defense response to virus', 'protein polyubiquitination', 'proteasome-mediated ubiquitin-dependent protein catabolic process')



# Filter the top3 terms in GO_geneList from all timepoints
GO_BP_3_geneList <- GO_geneList[(GO_geneList$Description %in% Terms),]


# Filter distinct terms each with the highest count
GO_BP_3_geneList <- GO_BP_3_geneList %>%
  arrange(desc(Count)) %>%
  distinct(Description, .keep_all = T)


# Get the Genes
Genes <- GO_BP_3_geneList %>% 
  transform(geneID = strsplit(geneID, "/")) %>% 
  unnest(geneID)


Genes_terms <- data.frame(geneID = Genes$geneID,
                          Description = Genes$Description)


Genes_terms <- Genes_terms %>%
  pivot_wider(names_from = Description, values_from = Description) %>%
  as.data.frame()


Genes_terms_1 <- merge(as.data.frame(Genes_terms), as.data.frame(DE_geneLogFC_Circmx), by.x = "geneID", by.y = "gene_name", all.x =T)
Genes_terms_1 <- Genes_terms_1 %>% filter(!is.na(CircRNA_ID))

## Set rownames
# rownames(Genes_terms_1) <- paste0(Genes_terms_1$CircRNA_ID, "/", Genes_terms_1$geneID)
# Genes_terms_1$circ_name <- paste0("circ", make.unique(Genes_terms_1$geneID, sep = "_"))
Genes_terms_1 <- Genes_terms_1[,c(2,4,3, 10)] %>% unique()
rownames(Genes_terms_1) <- paste0(Genes_terms_1$circ_name)
colnames(Genes_terms_1) <- c("BP Term 3", "BP Term 2", "BP Term 1", "circ_name")



## Add labels
# colnames(FPKM_DE_logs) <-colnames(FPKM_DE_log)
# rownames(FPKM_DE_logs) <- rep_FPKM$gene_name




############


# Filter out DE from rep FPKM 
logFC_BP <- DE_geneLogFC_Circmx %>% 
  as.data.frame(.) %>% 
  filter(gene_name %in% Genes$geneID) %>% as.data.frame()


rownames(logFC_BP) <- as.character(paste0(logFC_BP$circ_name))
logFC_BP <- logFC_BP %>% dplyr::select(c(2:5)) %>% data.matrix()


# Arrange col order

# col.order <- c("C3H_0", "C3H_1", "C6H_0", "C6H_1", "C12H_0", "C12H_1", "C18H_0", "C18H_1", 
#                "T3H_0", "T3H_1", "T6H_0", "T6H_1", "T12H_0", "T12H_1", "T18H_0", "T18H_1")
# FPKM_BP <- FPKM_BP[, col.order]


logFC_BP_HM <- pheatmap(logFC_BP %>% as.matrix(.), annotation_row = Genes_terms_1[,c(1,2,3)], 
                        cluster_rows = T,
                        labels_row = Genes_terms_1$circ_name,
               cluster_cols = FALSE, clustering_distance_rows = "euclidean",
               clustering_distance_cols = "euclidean", clustering_method = "complete",
               fontsize = 20, annotation_names_row = T,
               fontsize_row = 10, fontsize_col = 24, 
               #color=hmcols,
               scale="none", show_rownames = T,
               main = "Gene Expression of Interested GO BP terms", border_color=NA)
               # filename = "Gene Expression of Top 3 GO BP terms.png")


logFC_BP_HM

##########

# MF
# Get the top2 terms for MF
Terms <- c('ubiquitin-like protein ligase activity', 'aminoacyltransferase activity', 'ubiquitin-like protein transferase activity')



# Filter the top3 terms in GO_geneList from all timepoints
GO_MF_3_geneList <- GO_geneList[(GO_geneList$Description %in% Terms),]


# Filter distinct terms each with the highest count
GO_MF_3_geneList <- GO_MF_3_geneList %>%
  arrange(desc(Count)) %>%
  distinct(Description, .keep_all = T)


# Get the Genes
Genes <- GO_MF_3_geneList %>% 
  transform(geneID = strsplit(geneID, "/")) %>% 
  unnest(geneID)


Genes_terms <- data.frame(geneID = Genes$geneID,
                          Description = Genes$Description)


Genes_terms <- Genes_terms %>%
  pivot_wider(names_from = Description, values_from = Description) %>%
  as.data.frame()

Genes_terms_1 <- merge(as.data.frame(Genes_terms), as.data.frame(DE_geneLogFC_Circmx), by.x = "geneID", by.y = "gene_name", all.x =T)
Genes_terms_1 <- Genes_terms_1 %>% filter(!is.na(CircRNA_ID))

## Set rownames
# rownames(Genes_terms_1) <- paste0(Genes_terms_1$CircRNA_ID, "/", Genes_terms_1$geneID)
# Genes_terms_1$circ_name <- paste0("circ", make.unique(Genes_terms_1$geneID, sep = "_"))
Genes_terms_1 <- Genes_terms_1[,c(2,4,3, 10)] %>% unique()
rownames(Genes_terms_1) <- paste0(Genes_terms_1$circ_name)
colnames(Genes_terms_1) <- c("MF Term 2", "MF Term 3", "MF Term 1", "circ_name")





## Add labels
# colnames(FPKM_DE_logs) <-colnames(FPKM_DE_log)
# rownames(FPKM_DE_logs) <- rep_FPKM$gene_name




############


# Filter out DE from logFC 
logFC_MF <- DE_geneLogFC_Circmx %>% 
  as.data.frame(.) %>% 
  filter(gene_name %in% Genes$geneID) %>% as.data.frame()


rownames(logFC_MF) <- as.character(paste0(logFC_MF$circ_name))
logFC_MF <- logFC_MF %>% dplyr::select(c(2:5)) %>% data.matrix()



# Arrange col order

# col.order <- c("C3H_0", "C3H_1", "C6H_0", "C6H_1", "C12H_0", "C12H_1", "C18H_0", "C18H_1", 
#                "T3H_0", "T3H_1", "T6H_0", "T6H_1", "T12H_0", "T12H_1", "T18H_0", "T18H_1")
# FPKM_BP <- FPKM_BP[, col.order]


logFC_MF_HM <- pheatmap(logFC_MF, annotation_row = Genes_terms_1[,c(2, 1,3)], cluster_rows = T,
         labels_row = Genes_terms_1$circ_name,
         cluster_cols = FALSE, clustering_distance_rows = "euclidean",
         clustering_distance_cols = "euclidean", clustering_method = "complete",
         fontsize = 20, annotation_names_row = T,
         fontsize_row = 10, fontsize_col = 24, 
         #color=hmcols,
         scale="none", show_rownames = T,
         main = "Gene Expression of Interested GO MF terms", border_color=NA)
         #filename = "Gene Expression of Top 3 GO MF terms.png")

logFC_MF_HM
###########

### CC


# Get the top3 terms for BP
Terms <- c('cytoplasmic stress granule', 'lysosomal membrane', 'early endosome')



# Filter the top3 terms in GO_geneList from all timepoints
GO_CC_3_geneList <- GO_geneList[(GO_geneList$Description %in% Terms),]


# Filter distinct terms each with the highest count
GO_CC_3_geneList <- GO_CC_3_geneList %>%
  arrange(desc(Count)) %>%
  distinct(Description, .keep_all = T)


# Get the Genes
Genes <- GO_CC_3_geneList %>% 
  transform(geneID = strsplit(geneID, "/")) %>% 
  unnest(geneID)


Genes_terms <- data.frame(geneID = Genes$geneID,
                          Description = Genes$Description)


Genes_terms <- Genes_terms %>%
  pivot_wider(names_from = Description, values_from = Description) %>%
  as.data.frame()

Genes_terms_1 <- merge(as.data.frame(Genes_terms), as.data.frame(DE_geneLogFC_Circmx), by.x = "geneID", by.y = "gene_name", all.x =T)
Genes_terms_1 <- Genes_terms_1 %>% filter(!is.na(CircRNA_ID))

## Set rownames
# rownames(Genes_terms_1) <- paste0(Genes_terms_1$CircRNA_ID, "/", Genes_terms_1$geneID)
# Genes_terms_1$circ_name <- paste0("circ", make.unique(Genes_terms_1$geneID, sep = "_"))
Genes_terms_1 <- Genes_terms_1[,c(2,4,3, 10)] %>% unique()
rownames(Genes_terms_1) <- paste0(Genes_terms_1$circ_name)
colnames(Genes_terms_1) <- c("CC Term 2", "CC Term 1", "CC Term 3", "circ_name")



## Add labels
# colnames(FPKM_DE_logs) <-colnames(FPKM_DE_log)
# rownames(FPKM_DE_logs) <- rep_FPKM$gene_name




############


# Filter out DE from rep FPKM 
logFC_CC <- DE_geneLogFC_Circmx %>% 
  as.data.frame(.) %>% 
  filter(gene_name %in% Genes$geneID) %>% as.data.frame()


rownames(logFC_CC) <- as.character(paste0(logFC_CC$circ_name))
logFC_CC <- logFC_CC %>% dplyr::select(c(2:5)) %>% data.matrix()


# Arrange col order

# col.order <- c("C3H_0", "C3H_1", "C6H_0", "C6H_1", "C12H_0", "C12H_1", "C18H_0", "C18H_1", 
#                "T3H_0", "T3H_1", "T6H_0", "T6H_1", "T12H_0", "T12H_1", "T18H_0", "T18H_1")
# FPKM_CC <- FPKM_CC[, col.order]


logFC_CC_HM <- pheatmap(logFC_CC, annotation_row = Genes_terms_1[,c(3, 1,2)], cluster_rows = T,
                        labels_row = Genes_terms_1$circ_name,
                        cluster_cols = FALSE, clustering_distance_rows = "euclidean",
                        clustering_distance_cols = "euclidean", clustering_method = "complete",
                        fontsize = 20, annotation_names_row = T,
                        fontsize_row = 10, fontsize_col = 24, 
                        #color=hmcols,
                        scale="none", show_rownames = T,
                        main = "Gene Expression of Interested GO CC terms", border_color=NA)
# filename = "Gene Expression of Top 3 GO CC terms.png")
logFC_CC_HM


#######################################################################################################################################################################################

## KEGG

# Enrich KEGG
library(R.utils)
R.utils::setOption("clusterProfiler.download.method","auto")

## 12hrs
keg_res_12 <- enrichKEGG(Circ_SD_DE_12hrs$ENTREZID, organism= "hsa", keyType= "kegg",
                         pvalueCutoff = 0.05, pAdjustMethod = "BH")
## 18hrs
keg_res_18 <- enrichKEGG(Circ_SD_DE_18hrs$ENTREZID, organism= "hsa", keyType= "kegg",
                         pvalueCutoff = 0.05, pAdjustMethod = "BH")


# Change to gene name
keg_res_12 <- setReadable(keg_res_12, OrgDb = org.Hs.eg.db, keyType="ENTREZID")
keg_res_18 <- setReadable(keg_res_18, OrgDb = org.Hs.eg.db, keyType="ENTREZID")

## Write kegg enrichment info
library(xlsx)
write.xlsx2(keg_res_12@result, file = "KEGG_enrichment_12hrs.xlsx") 

write.xlsx2(keg_res_18@result, file = "KEGG_enrichment_18hrs.xlsx")


## KEGG enrichment
# Bind the KEGG results from different timepoints
KEGG_rbind <- rbind(keg_res_12@result %>% mutate(Timepoint = "12hrs"), 
                    keg_res_18@result %>% mutate(Timepoint = "18hrs")) %>%
  dplyr::select(c("Description", "p.adjust", "Timepoint"))


# Change from long to wide
KEGG_wide_All <- KEGG_rbind %>% 
  group_by(Description) %>%
  pivot_wider(names_from = Timepoint, values_from = p.adjust) %>%
  column_to_rownames("Description")


# Replace NA with max values
# KEGG_wide_All$`3hrs` <- KEGG_wide_All$`3hrs` %>% replace(is.na(.), max(keg_res_3@result$p.adjust))
# KEGG_wide_All$`6hrs` <- KEGG_wide_All$`6hrs` %>% replace(is.na(.), max(keg_res_6@result$p.adjust))
# KEGG_wide_All$`12hrs` <- KEGG_wide_All$`12hrs` %>% replace(is.na(.), max(keg_res_12@result$p.adjust))
# KEGG_wide_All$`18hrs` <- KEGG_wide_All$`18hrs` %>% replace(is.na(.), max(keg_res_18@result$p.adjust))

Breaks <- c(0, 0.05, 0.1, 0.2, 0.3, 0.4, 1.0)
# Color <- colorRampPalette(rev(brewer.pal(n = 4, name = "RdYlBu")))(8)
# pheatmap(KEGG_wide_All, cellwidth = 15, cellheight = 12, color = Color,
#          legend_breaks = Breaks, legend_labels = c(0, 0.05, 0.1, 0.2, 0.6, 0.8, 1.0),
#          cluster_cols = F, cluster_rows = F, filename = "KEGG_wide_All_pheatmap.png")


## Top 10
# KEGG_All
# KEGG_10 <- KEGG_rbind %>% 
#   arrange(p.adjust) %>%
#   distinct(Description) %>%
#   slice(1:10) %>%
#   as.list()


KEGG_wide_All_1 <- KEGG_wide_All %>% as.data.frame()
KEGG_wide_All_1$`3hrs` <- NA
KEGG_wide_All_1$`6hrs` <- NA

KEGG_wide_All <- KEGG_wide_All_1 %>%  dplyr::select(3, 4, 1, 2) %>% slice(1:10)


# KEGG_All_10 <- KEGG_wide_All[(rownames(KEGG_wide_All) %in% KEGG_10[["Description"]]),]

KEGG_PW <- pheatmap(KEGG_wide_All, cellwidth = 40, cellheight = 40, 
         legend_breaks = seq(min(KEGG_rbind$p.adjust), max(KEGG_rbind$p.adjust), length.out = 7), 
         legend_labels = c(0, 0.05, 0.06, 0.07, 0.075, 0.8, 0.9),
         cluster_cols = F, cluster_rows = F,
         #treeheight_row = 40, treeheight_col = 40,
         clustering_distance_rows = "euclidean",
         clustering_distance_cols = "euclidean", clustering_method = "complete",
         fontsize = 14,
         fontsize_col = 20, fontsize_row = 18,
         border_color='grey60', show_rownames = T,
         labels_row = make_bold_names(KEGG_wide_All, rownames, rownames(KEGG_wide_All)[c(1,2,4)]),
         labels_col = make_bold_names(KEGG_wide_All, colnames, colnames(KEGG_wide_All)),
         main = "Top 10 KEGG Enriched Pathways")
         #filename = "KEGG_All_10_pheatmap.png")





################################################################
## KEGG geneList logFC expression
# Get the gene list and the counts for each timepoint

KEGG_geneList <- rbind(keg_res_12@result %>% mutate(Timepoint = "12hrs"), 
                       keg_res_18@result %>% mutate(Timepoint = "18hrs")) %>%
  dplyr::select(c("Description", "geneID", "Count","Timepoint"))


# Get the top3 terms for BP
Terms <- c("TNF signaling pathway", "Endocytosis", "Chemokine signaling pathway")


# Filter the top3 terms in GO_geneList from all timepoints

KEGG_3_geneList <- KEGG_geneList[(KEGG_geneList$Description %in% Terms),]


# Filter distinct terms each with the highest count

KEGG_3_geneList <- KEGG_3_geneList %>%
  arrange(desc(Count)) %>%
  distinct(Description, .keep_all = T)

# Get the Genes

Genes <- KEGG_3_geneList %>% 
  transform(geneID = strsplit(geneID, "/")) %>% 
  unnest(geneID)


Genes_terms <- data.frame(geneID = Genes$geneID,
                          Description = Genes$Description)


Genes_terms <- Genes_terms %>%
  pivot_wider(names_from = Description, values_from = Description) %>%
  as.data.frame()

Genes_terms_1 <- merge(as.data.frame(Genes_terms), as.data.frame(DE_geneLogFC_Circmx), by.x = "geneID", by.y = "gene_name", all.x =T)
Genes_terms_1 <- Genes_terms_1 %>% filter(!is.na(CircRNA_ID))

## Set rownames
# rownames(Genes_terms_1) <- paste0(Genes_terms_1$CircRNA_ID, "/", Genes_terms_1$geneID)
# Genes_terms_1$circ_name <- paste0("circ", make.unique(Genes_terms_1$geneID, sep = "_"))
Genes_terms_1 <- Genes_terms_1[,c(2,4,3, 10)] %>% unique()
rownames(Genes_terms_1) <- paste0(Genes_terms_1$circ_name)
colnames(Genes_terms_1) <- c("KEGG Pathway 3", "KEGG Pathway 2", "KEGG Pathway 1", "circ_name")






## Add labels
# colnames(FPKM_DE_logs) <-colnames(FPKM_DE_log)
# rownames(FPKM_DE_logs) <- rep_FPKM$gene_name




############
# Filter out DE from rep FPKM 
logFC_KEGG <- DE_geneLogFC_Circmx %>% 
  as.data.frame(.) %>% 
  filter(gene_name %in% Genes$geneID) 


rownames(logFC_KEGG) <- as.character(paste0(logFC_KEGG$circ_name))
logFC_KEGG <- logFC_KEGG %>% dplyr::select(c(2:5)) %>% data.matrix()


# Arrange col order

col.order <- c("3hrs", "6hrs", "12hrs", "18hrs")
logFC_KEGG <- logFC_KEGG[, col.order]


# col.order <- c("C3H_0", "C3H_1", "C6H_0", "C6H_1", "C12H_0", "C12H_1", "C18H_0", "C18H_1", 
#                "T3H_0", "T3H_1", "T6H_0", "T6H_1", "T12H_0", "T12H_1", "T18H_0", "T18H_1")
# FPKM_BP <- FPKM_BP[, col.order]


KEGG_logFC_expression <- pheatmap(logFC_KEGG, annotation_row = Genes_terms_1[,c(1,2,3)], cluster_rows = T,
         labels_row = Genes_terms_1$circ_name,
         cluster_cols = FALSE, clustering_distance_rows = "euclidean",
         clustering_distance_cols = "euclidean", clustering_method = "complete",
         fontsize = 20, annotation_names_row = T,
         fontsize_row = 10, fontsize_col = 24, 
         #color=hmcols,
         scale="none", show_rownames = T,
         main = "Gene Expression of Interested KEGG Pathways", border_color=NA)
         # filename = "Gene Expression of Top 3 KEGG terms.png")

KEGG_logFC_expression


## Save plots
library(ggpubr)
library(ggplotify)
g1 <- as.ggplot(KEGG_PW)
g2 <- as.ggplot(KEGG_logFC_expression)

# use patchwork to arrange them together
figure <- ggarrange(g1, g2,
                    #labels = c("A", "B"), 
                    #heights = c(3,1), 
                    #label.y = 0.99, label.x = 0.13,
                    ncol = 2, nrow = 1)
figure
