# CircRNA real_data_analysis
CircRNA analysis (Identification, Quantification, Characterisation, Functional Annotation)

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

## MiRNA-mRNA prediction

To mimic the previous circRNA-miRNA prediction, we have **CUSTOM MODIFIED** the CircR.py script to enable the input to be mRNA bed file instead of circRNA. In this case, we use the 3' UTR mRNA bed12 file obtained from UCSC. This is because miRNA predominantly binds to the 3' UTR region of human mRNAs.

Example run:

```
python miRNA-mRNA.py --threads 8 -i miRNA-mRNA_input12.bed --coord --gtf ./support_files/human/hg38/hg38.ensGene.gtf \
--genome ./support_files/human/hg38/hg38.fa --rRNA ./support_files/human/hg38/hg38.rRNA.bed \
--miRNA miRNA_miRanda.fa --AGO ./support_files/human/hg38/hg38.AGO.bed \
--validated_interactions ./support_files/human/hg38/hg38.INT.bed --TS_miRNA_lib miRNA_TS.txt -o miRNA_mRNA_results.csv

```
