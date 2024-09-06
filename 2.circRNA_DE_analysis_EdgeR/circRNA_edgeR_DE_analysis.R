library(tidyr)
library(dplyr)
library(ggplot2)
library(ggpubr)
library(stringr)
library(pheatmap)
library(colorspace)
library(dendextend)
library(gplots)
library(edgeR)
library(limma)
library(tidyverse)
library(magrittr)
library(statmod)
library(reshape2)
library(affy)

# setwd('C:/Users/matha/Documents/CircRNA_Analysis_Compile/CIRIquant_Superdepth_noStringency/Allhrs_Superdepth_NoStringency/')

# setwd('D:/circRNA/CIRIquant_Superdepth_noStringency')
# 
# #load data
# lib_mtx_03 <- read.csv("3hrs_Superdepth_noStringency/library_info.csv", row.names = 1) %>% mutate(Timepoint = "03hrs")
# lib_mtx_06 <- read.csv("6hrs_Superdepth_noStringency/library_info.csv", row.names = 1) %>% mutate(Timepoint = "06hrs")
# lib_mtx_12 <- read.csv("12hrs_Superdepth_noStringency/library_info.csv", row.names = 1) %>% mutate(Timepoint = "12hrs")
# lib_mtx_18 <- read.csv("18hrs_Superdepth_noStringency/library_info.csv", row.names = 1) %>% mutate(Timepoint = "18hrs")
# 
# #put all data frames into list
# lib_mtx <- rbind(lib_mtx_03, lib_mtx_06, lib_mtx_12, lib_mtx_18)
# rownames(lib_mtx) <- c("03hrs_CASE1", "03hrs_CASE2", "03hrs_CONTROL1", "03hrs_CONTROL2",
#                        "06hrs_CASE1", "06hrs_CASE2", "06hrs_CONTROL1", "06hrs_CONTROL2",
#                        "12hrs_CASE1", "12hrs_CASE2", "12hrs_CONTROL1", "12hrs_CONTROL2",
#                        "18hrs_CASE1", "18hrs_CASE2", "18hrs_CONTROL1", "18hrs_CONTROL2")
# 
# 
# 
# #gene_mtx
# gene_mtx_03 <- read.csv("3hrs_Superdepth_noStringency/gene_count_matrix.csv", row.names = 1)
# gene_mtx_06 <- read.csv("6hrs_Superdepth_noStringency/gene_count_matrix.csv", row.names = 1)
# gene_mtx_12 <- read.csv("12hrs_Superdepth_noStringency/gene_count_matrix.csv", row.names = 1)
# gene_mtx_18 <- read.csv("18hrs_Superdepth_noStringency/gene_count_matrix.csv", row.names = 1)
# 
# #put all data frames into list
# gene_mtx_36 <- merge(gene_mtx_03, gene_mtx_06, by = "row.names", all = T)
# gene_mtx_3612 <- merge(gene_mtx_36, gene_mtx_12, by.x = "Row.names", by.y = "row.names", all = T)
# gene_mtx <- merge(gene_mtx_3612, gene_mtx_18, by.x = "Row.names", by.y = "row.names", all = T) %>%
#   tibble::column_to_rownames(var = "Row.names")
# gene_mtx[is.na(gene_mtx)] <- 0
# 
# colnames(gene_mtx) <- c(rownames(lib_mtx))
# #gene_mtx <- gene_mtx[,rownames(lib_mtx)]
# 
# 
# 
# ##bsj_mtx
# bsj_mtx_03 <- read.csv("3hrs_Superdepth_noStringency/circRNA_bsj.csv", row.names=1, sep = ",")
# bsj_mtx_06 <- read.csv("6hrs_Superdepth_noStringency/circRNA_bsj.csv", row.names=1, sep = ",")
# bsj_mtx_12 <- read.csv("12hrs_Superdepth_noStringency/circRNA_bsj.csv", row.names=1, sep = ",")
# bsj_mtx_18 <- read.csv("18hrs_Superdepth_noStringency/circRNA_bsj.csv", row.names=1, sep = ",")
# 
# #put all data frames into list
# bsj_mtx_36 <- merge(bsj_mtx_03, bsj_mtx_06, by = "row.names", all = T)
# bsj_mtx_3612 <- merge(bsj_mtx_36, bsj_mtx_12, by.x = "Row.names", by.y = "row.names", all = T)
# bsj_mtx <- merge(bsj_mtx_3612, bsj_mtx_18, by.x = "Row.names", by.y = "row.names", all = T) %>%
#   tibble::column_to_rownames(var = "Row.names")
# bsj_mtx[is.na(bsj_mtx)] <- 0
# 
# colnames(bsj_mtx) <- c(rownames(lib_mtx))

### Allhrs compiled from step2 and step3 CIRIquant
#load data

samples <- c("CASE1", "CASE2", "CONTROL1", "CONTROL2",
             "CASE3", "CASE4", "CONTROL3", "CONTROL4",
             "CASE5", "CASE6", "CONTROL5", "CONTROL6",
             "CASE7", "CASE8", "CONTROL7", "CONTROL8")

lib_mtx <- read.csv('library_info.csv', row.names = 1)
lib_mtx <- lib_mtx[samples,]
lib_mtx$Timepoint <- c(rep('03hrs',4), rep('06hrs',4), rep('12hrs',4), rep('18hrs',4) )
rownames(lib_mtx) <- c("03hrs_CASE1", "03hrs_CASE2", "03hrs_CONTROL1", "03hrs_CONTROL2",
                       "06hrs_CASE1", "06hrs_CASE2", "06hrs_CONTROL1", "06hrs_CONTROL2",
                       "12hrs_CASE1", "12hrs_CASE2", "12hrs_CONTROL1", "12hrs_CONTROL2",
                       "18hrs_CASE1", "18hrs_CASE2", "18hrs_CONTROL1", "18hrs_CONTROL2")

gene_mtx <- read.csv('gene_count_matrix.csv', row.names = 1)
gene_mtx <- gene_mtx[,samples]
colnames(gene_mtx) <- c(rownames(lib_mtx))

bsj_mtx <- read.csv('circRNA_bsj.csv', row.names=1, sep = ",")
bsj_mtx <- bsj_mtx[,samples]
colnames(bsj_mtx) <- c(rownames(lib_mtx))

# Filtering and normalization
gene_DGE <- DGEList(counts = gene_mtx, group = lib_mtx$Group, samples = lib_mtx$Timepoint) ## Create one for each timepoint
# gene_DGE <- edgeR::cbind.DGEList(gene_DGE_3hrs, gene_DGE_6hrs, gene_DGE_12hrs, gene_DGE_18hrs) # combine all

# gene_idx <- filterByExpr(gene_DGE, design=design, group = )
# gene_DGE <- gene_DGE[gene_idx, keep.lib.sizes=FALSE]
# 
# 
# counts <- gene_DGE$counts %>% cpm() %>% as.data.frame() %>% rownames_to_column()

genes2keep  <- gene_DGE %>% cpm() %>%
            is_greater_than(1) %>%
            rowSums() %>%
            is_greater_than(4)

gene_DGE <- gene_DGE[genes2keep, keep.lib.sizes=FALSE]
gene_DGE <- calcNormFactors(gene_DGE) # apply TMM normalization to account for the compositional biases



###################################################################################################################3
## Boxplot 

par(mar=c(9,2,1,1))

boxplot(gene_DGE %>% cpm(log=F) %>% +1 %>% log10(), xlab="", ylab="Log2 counts per million",las=2)
# Let's add a blue horizontal line that corresponds to the median logCPM
abline(h=median(gene_DGE %>% cpm(log=F) %>% +1 %>% log10()),col="blue")
title("Boxplots of log(CPM + 1) (Normalised)")

#######
log10cpm <- gene_DGE %>% cpm(log=F) %>% log10() %>%  as.data.frame %>% gather(key = sample, value = val)

ggboxplot(log10cpm, x = "sample", y = "val", width = 0.8)

# ## Boxplot
Boxplot_1 = ggplot(log10cpm , aes(x=sample, y= as.numeric(val), fill=sample)) +
  labs(x = "Samples", y = "log10(CPM)", title = "Boxplot of Human Linear Genes (H3N2)", fill = "Sample_names") +
  geom_boxplot() + 
  theme(axis.text.x = element_blank(),
        axis.title.x = element_blank(),
        axis.text.y = element_text(size = 20),
        axis.title = element_text(size = 20),
        plot.title = element_text(size = 25, hjust = 0.5),
        legend.position = "right")
Boxplot_1
###########################################################################################################

## Reads mapped per sample
par(mar=c(9,5,5, 5))

colSums(as.matrix(gene_DGE)) %>% barplot(., ylab="Reads mapped per sample", las=2)

############################################################################################################

## Must be normal distributed graph
gene_DGE %>%
  cpm(log = TRUE) %>%
  plotDensities(legend = FALSE, main = "B. After Filtering")

plotDensity(gene_DGE %>% cpm(log = TRUE), lwd=2,lty=1,xlab("Density"),ylab("Counts"), main = "Density plot")



#############################################################################################

## Barplot on the number of genes detected

## Sample vs genes
par(mar=c(9,5,5,5))
{barplot(colSums(gene_DGE %>% cpm(log=F) >2),ylab="Number of detected genes",las=2)
  abline(h=median(colSums(gene_DGE %>% cpm(log=F) >2)))}

## Genes vs samples
{barplot(rowSums(gene_DGE %>% cpm(log=F)>2),xlab="Genes",ylab="Number of samples",las=2,names.arg="")
  abline(h=median(rowSums(gene_DGE %>% cpm(log=F)>2)))}

######################################################################################################

## MDS plot
samples <- c("red", "red", "orange", "orange",
             "darkgreen", "darkgreen", "green", "green",
             "darkblue", "darkblue","lightblue", "lightblue",
             "purple", "purple", "yellow", "yellow")
plotMDS(gene_DGE %>% cpm(log=F) %>% +1 %>% log10(), col=samples)

# plotMDS(gene_DGE,dim=c(1,2))


## MD plot
plotMD(gene_DGE, column=10)
abline(h=0, col="red", lty=2, lwd=2)

######################################################################################################

## Mean difference plots to check whether composition bias were addressed

plotMA(gene_DGE, column=8)
abline(h=0, col="red", lty=2, lwd=2)

############################################################################################################################
## PCA

library(factoextra)

par(mar=c(5,5,1,1))

cpm_log <- cpm(gene_DGE, log = F) %>% +1 %>% log10()

pca <- princomp(cpm_log)
fviz_eig(pca, addlabels = TRUE)
fviz_pca_var(pca, col.var = "black")


# pca <- prcomp(t(cpm_log), scale. = TRUE)
# plot(pca$x[, 1], pca$x[, 2], pch = ".", xlab = "PC1", ylab = "PC2")
# text(pca$x[, 1], pca$x[, 2], labels = colnames(cpm_log))
# summary(pca)

############################################################################################################################
############################################################################################################################
### Heatmap
hmcols <- rev(redgreen(2750));



# Arrange col order
counts_log <- gene_DGE %>% cpm(log = F) %>% +1 %>% log10()

col.order <- c("03hrs_CONTROL1", "03hrs_CONTROL2", "06hrs_CONTROL1", "06hrs_CONTROL2", "12hrs_CONTROL1", "12hrs_CONTROL2", "18hrs_CONTROL1", "18hrs_CONTROL2", 
               "03hrs_CASE1", "03hrs_CASE2", "06hrs_CASE1", "06hrs_CASE2", "12hrs_CASE1", "12hrs_CASE2", "18hrs_CASE1", "18hrs_CASE2")
counts_log <- counts_log[, col.order]


pheatmap(counts_log,
         cluster_rows = TRUE,
         cluster_cols = FALSE, clustering_distance_rows = "euclidean",
         clustering_distance_cols = "euclidean", clustering_method = "complete",
         cutree_rows = 1, cutree_cols = 1, fontsize = 12, 
         fontsize_row = 8, color=hmcols,scale="row", show_rownames = F,
         main = "Human H3N2 Linear Gene Heatmap")


############################################################################################################################
############################################################################################################################

## dendogram

dist <- dist(t(log10(gene_df + 1)) , diag=TRUE)

# Hierarchical Clustering with hclust
hc <- hclust(dist)

# Plot the result
plot(hc)



######################################################################################################

## Design matrix

Time <- factor(lib_mtx$Timepoint, levels = c("03hrs", "06hrs", "12hrs", "18hrs"))
treat <- factor(lib_mtx$Group, levels=c("C", "T"))

design <- model.matrix(~Time + Time:treat)


######################################################################################################3
gene_DGE <- estimateDisp(gene_DGE, design, robust = TRUE)
gene_fit <- glmFit(gene_DGE, design)
gene_lrt <- glmLRT(gene_fit, coef=7)

gene_df <- gene_lrt$table
gene_order <- order(gene_lrt$table$PValue)
gene_df$DE <- decideTestsDGE(gene_lrt)
gene_df <- gene_df[gene_order, ]
gene_df$FDR <- p.adjust(gene_df$PValue, method="fdr")
write.csv(gene_df, file="New_design_linear_mRNA_18hrs.tsv", quote = FALSE)

with(gene_df, plot(logCPM, logFC, pch=20, cex=.8,  main="logFC vs Abundance"))
with(subset(gene_df, FDR<0.05), points(logCPM, logFC, cex=.8 ,pch=20, col="red"))


######################################################################################################
# # design matrix
# if ("Subject" %in% colnames(lib_mtx)) {
#   subject <- factor(lib_mtx$Subject)
#   treat <- factor(lib_mtx$Group, levels=c("C", "T"))
# 
#   design <- model.matrix(~subject + treat)
# } else {
#   treat <- factor(lib_mtx$Group, levels=c("C", "T"))
#   design <- model.matrix(~treat)
# }

#######################################################################################################

#linear gene DE analysis
# gene_df <- gene_lrt$table
# gene_order <- order(gene_lrt$table$PValue)
# gene_df$DE <- decideTestsDGE(gene_lrt)
# gene_df <- gene_df[gene_order, ]
# gene_df$FDR <- p.adjust(gene_df$PValue, method="fdr")


######################################################################################################

# create DGEList object from the count table
circ_DGE <- DGEList(counts = bsj_mtx,
                    group = lib_mtx$Group,
                    lib.size = gene_DGE$samples[, "lib.size"],
                    norm.factors = gene_DGE$samples[, "norm.factors"])


circ_DGE_cpm <- circ_DGE %>% edgeR::cpm() %>%
  is_greater_than(0.05) %>%
  rowSums() %>%
  is_greater_than(2)


## keep circRNA of interest
circ_DGE <- circ_DGE[circ_DGE_cpm, keep.lib.sizes=TRUE]


# ## write cpm to file
# library(xlsx)
# counts <- cpm(circ_DGE, normalized.lib.sizes = T)
# write.xlsx(counts, file="circ_cpm_new.xlsx")


# circ_DGE[["samples"]]$replicates <- rep(c(1,2), 8)

# circRNAs to keep
# circ2keep  <- circ_DGE %>% cpm() %>%
#   rowSums() %>%
#   is_greater_than(0.1)

# circ_DGE <- circ_DGE[circ2keep, keep.lib.sizes=T]
# 
# circ_counts <- circ_DGE$counts %>% as.data.frame() %>% rownames_to_column()

# circ_idx <- filterByExpr(circ_DGE, min.count=)
# circ_DGE <- circ_DGE[circ_idx, , keep.lib.sizes=TRUE]
# head(circ_df[rowSums(bsj_mtx >= 2) >= nrow(lib_mtx) / 2,])



############################################################################################################
## Boxplot


colnames(circ_DGE) <- c("03hrs_Trt_1", "03hrs_Trt_2", "03hrs_Ctrl_1", "03hrs_Ctrl_2",
                             "06hrs_Trt_1", "06hrs_Trt_2", "06hrs_Ctrl_1", "06hrs_Ctrl_2",
                             "12hrs_Trt_1", "12hrs_Trt_2", "12hrs_Ctrl_1", "12hrs_Ctrl_2",
                             "18hrs_Trt_1", "18hrs_Trt_2", "18hrs_Ctrl_1", "18hrs_Ctrl_2")

counts_log <- cpm(circ_DGE) %>% +1 %>% log10() %>% as.data.frame()
counts_log$circRNA <- rownames(counts_log)
counts_log_boxplot <- melt(counts_log, id.vars=c("circRNA"))
counts_log_boxplot$variable <- factor(counts_log_boxplot$variable, levels = c("03hrs_Ctrl_1", "03hrs_Ctrl_2", "06hrs_Ctrl_1", "06hrs_Ctrl_2",
                                                                              "12hrs_Ctrl_1", "12hrs_Ctrl_2", "18hrs_Ctrl_1", "18hrs_Ctrl_2",
                                                                              "03hrs_Trt_1", "03hrs_Trt_2", "06hrs_Trt_1", "06hrs_Trt_2", 
                                                                              "12hrs_Trt_1", "12hrs_Trt_2", "18hrs_Trt_1", "18hrs_Trt_2"))

ggplot(counts_log_boxplot, aes(x = variable, y = value, fill = variable)) +
  geom_boxplot() +
  # theme(aspect.ratio = 1) +
  labs(x = "Samples", y = "log10(CPM)") +
  theme(axis.text.x = element_text(angle = 65, vjust = 0.5, hjust=0.5,
                                   face = "bold", size="14"),
        axis.text.y = element_text(size = 14, face = "bold", vjust = 1, hjust = 1),
        axis.title.x = element_text(size = 15, face= "bold"),
        axis.title.y = element_text(size = 15, face = "bold"),
        legend.position  = "none")


############################################################################################################
# # 
# boxplot(counts_log, xlab="", ylab="Log2 counts per million",las=2)
# # Let's add a blue horizontal line that corresponds to the median logCPM
# abline(h=median(counts_log),col="blue")
# title("Boxplots of logCPMs (normalised)")

############################################################################################################
## Heatmap

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


#  Define colours
hmcols <- rev(redgreen(2750));


# Arrange col order
col.order <- c("03hrs_Ctrl_1", "03hrs_Ctrl_2", "06hrs_Ctrl_1", "06hrs_Ctrl_2",
               "12hrs_Ctrl_1", "12hrs_Ctrl_2", "18hrs_Ctrl_1", "18hrs_Ctrl_2",
               "03hrs_Trt_1", "03hrs_Trt_2", "06hrs_Trt_1", "06hrs_Trt_2", 
               "12hrs_Trt_1", "12hrs_Trt_2", "18hrs_Trt_1", "18hrs_Trt_2")
## counts_log <- log(cpm + 1)
counts_log <- counts_log[, col.order]


pheatmap(counts_log,
         cluster_rows = TRUE,
         cluster_cols = FALSE, clustering_distance_rows = "euclidean",
         clustering_distance_cols = "euclidean", clustering_method = "complete",
         cutree_rows = 1, cutree_cols = 1, fontsize = 18, face="bold",
         fontsize_row = 8, color=hmcols,scale="row", show_rownames = F,
         labels_col = make_bold_names(counts_log, colnames, col.order))




############################################################################################################

# ## Must be normal distributed graph
# circ_DGE %>%
#   cpm(log = TRUE) %>%
#   plotDensities(legend = FALSE, main = "B. After Filtering")

############################################################################################################

## MDS plot
samples <- c("red", "red", "orange", "orange",
             "darkgreen", "darkgreen", "green", "green",
             "darkblue", "darkblue","lightblue", "lightblue",
             "purple", "purple", "yellow", "yellow")
plotMDS(circ_DGE, col=samples, dim = c(1,2))

############################################################################################################

#### MA plot
par(mfrow=c(1,2))
plotMD(circ_DGE,column = 7)
abline(h=0,col="grey")


############################################################################################################
#### PCA plot
library(factoextra)

par(mar=c(5,5,1,1))

cpm_log <- cpm(circ_DGE, log = F) %>% +1 %>% log10()

pca <- princomp(cpm_log)
fviz_eig(pca, addlabels = TRUE)
fviz_pca_var(pca, col.var = "black")

fviz_pca_biplot(pca,
                label="var", labelsize = 6, arrowsize = 0.7,
                repel = F,
                col.var = factor(c("C3H_0", "C3H_1", "C6H_0", "C6H_1", "C12H_0", "C12H_1", "C18H_0", "C18H_1", 
                                   "T3H_0", "T3H_1", "T6H_0", "T6H_1", "T12H_0", "T12H_1", "T18H_0", "T18H_1"))) +
  xlim(0,0.4) +
  ylim(-0.1, 0.1) +
  theme(axis.text = element_text(size = 12),
        text = element_text(size = 12),
        axis.title = element_text(size = 15, face = "bold"))

fviz_pca_biplot(pca,
                #individuals
                geom.ind="point",
                col.ind=Iris$Species,
                pointshape=Iris$Site,
                
                #variables
                col.var="steelblue")


############################################################################################################


circ_DGE <- estimateDisp(circ_DGE, design, robust = TRUE) # estimate the genewise dispersion estimates over all genes, allowing for a possible abundance trend
#circ_DGE$common.dispersion # display common dispersion value
#plotBCV(circ_DGE) # plot BCV plot


circ_fit <- glmFit(circ_DGE, design) # Fit a negative binomial generalized log-linear model to the read counts for each gene. 

## coef = 5 = 3hrs, coef = 6 = 6hrs ...
circ_lrt <- glmLRT(circ_fit, coef=5) # Conduct genewise statistical tests for a given coefficient or coefficient contrast.
circ_lrt # display circ_lrt datasets
circ_df <- circ_lrt$table
circ_order <- order(circ_lrt$table$PValue)
circ_df$DE <- as.vector(decideTests(circ_lrt))
circ_df <- circ_df[circ_order, ]
circ_df$FDR <- p.adjust(circ_df$PValue, method="fdr")
write.csv(cpm(circ_DGE), file="03hrs_cpm.csv", quote = F)
write.csv(circ_df, file="New_design_circRNA_03hrs.tsv", quote = FALSE)


# ## circ boxplot
# 
# cpm_long <- cpm %>% as.data.frame() %>% pivot_longer(cols = c("CASE1", "CASE2", "CONTROL1", "CONTROL2"),
#                                         names_to = "Dataset",
#                                         values_to = "CPM")
# 
# 
# 
# 
# ## To output a table
# circ_df <- circ_lrt$table
# circ_order <- order(circ_lrt$table$PValue)
# circ_df$DE <- as.vector(decideTestsDGE(circ_lrt))
# circ_df <- circ_df[circ_order, ]
# circ_df$FDR <- p.adjust(circ_df$PValue, method="fdr")
# write.csv(circ_df, file=opt$out, quote = FALSE)


