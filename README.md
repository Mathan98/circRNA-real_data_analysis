# CircRNA real_data_analysis
CircRNA analysis (Identification, Quantification, Characterisation, Functional Annotation)



# MRNA quality control analysis, GO and KEGG enrichment
Primarily, aligned RNA-seq reads by Cufflinks suite coupled with TopHat2 aligner are inspected through quality control metrics to ensure good alignment rate to the genome and their distribution across samples. 
DE mRNAs are determined when the |log2 fold-change| >= 2 and q-value <= 0.05
Moreover, mRNAs are also inspected for an innate immune response reaction due to IAV infection through GO and KEGG analysis across 4 timepoints. 

Scripts can be found in **1.mRNA_QC_GO_KEGG_Heatmap** folder



# CircRNA DE analysis through edgeR
CircRNAs are detected through CIRI2 and CIRIquant, where the differential expression (DE) analysis script was customised from the one provided by CIRIquant.
The DE of circRNAs are inspected through statistical tests across control and IAV-infected samples in a time-series analysis to retrieve confident DE circRNAs. In this case , we established that DE circRNAs have |log2 fold-change| >= 2 and FDR <= 0.05

Scripts can be found in **2.circRNA_DE_analysis_EdgeR** folder



# CircRNA characterization
CircRNAs are annotated based on two circular databases circAtlas and circBase. Both these databases represent the first established and current updated databases. 
CircRNAs are characterized based on circular host gene origin type and number of host genes.

Scripts can be found in **3.circRNA_characterization** folder



# CircRNA GO and KEGG with pheatmap of significance and log2FC expression across timepoints
By using the host genes of DE circRNAs, GO and KEGG enrichment were conducted using the clusterProfiler R package. This is to determine the functionally enriched terms and pathways of the DE circRNAs.
Then, a heatmap across the timepoints were plotted using the significance value of GO and KEGG output followed by picking the top 3 relevant terms or pathways and using the log2FC expression of circRNAs to plot another heatmap.

Scripts can be found in **4.circRNA_GO_KEGG_pheatmap** folder



# CircRNA functions as microRNA sponges that inadvertently regulates mRNA expression
## CircRNA-miRNA prediction

CircR is a tool used to predict the circRNA-miRNA interactions given user input circRNA coordinates in a bed file and supplementary files provided by CircR developers ([drive folder](https://drive.google.com/drive/folders/1zJVyzEFAMtvZTTueWRocxXs63jUxsl-U?usp=sharing))

Example run:

```
python Circr.py --threads 8 -i CircRNA_Overlap12_18hrs.bed --gtf ./support_files/human/hg38/hg38.ensGene.gtf \
--genome ./support_files/human/hg38/hg38.fa --rRNA ./support_files/human/hg38/hg38.rRNA.bed \
--miRNA ./support_files/miRNA/hsa_mature.fa --AGO ./support_files/human/hg38/hg38.AGO.bed \
--validated_interactions ./support_files/human/hg38/hg38.INT.bed -o overlap_12_18_miRanda_Nostringent.csv

```
Scripts can be found in **5.circRNA-miRNA_analysis** folder



## MiRNA-mRNA prediction

To mimic the previous circRNA-miRNA prediction, we have **CUSTOM MODIFIED** the CircR.py script to enable the input to be mRNA bed file instead of circRNA. In this case, we use the 3' UTR mRNA bed12 file obtained from UCSC. This is because miRNA predominantly binds to the 3' UTR region of human mRNAs.

Example run:

```
python miRNA-mRNA.py --threads 8 -i miRNA-mRNA_input12.bed --coord --gtf ./support_files/human/hg38/hg38.ensGene.gtf \
--genome ./support_files/human/hg38/hg38.fa --rRNA ./support_files/human/hg38/hg38.rRNA.bed \
--miRNA miRNA_miRanda.fa --AGO ./support_files/human/hg38/hg38.AGO.bed \
--validated_interactions ./support_files/human/hg38/hg38.INT.bed --TS_miRNA_lib miRNA_TS.txt -o miRNA_mRNA_results.csv

```

Scripts can be found in **6.miRNA-mRNA_analysis** folder

