
Acute lymphoblastic leukemia (ALL) is a type of cancer in which the bone marrow makes too many immature lymphocytes caused by Somatic mutations which are genetic variants that only accumulate in the tumor or the affected tissues, and they are not transmitted from generation to generation. 
# Aim of mini-project
In this project we ll be analysing DNA sequences,corresponding of regions of Chromosome (9,22) from 2 patients that were diagnosed with ALL-B inorder to identify mutations caused by somatic variation that may affect proteins functions or specific alterations that may be susceptible to tailored therapy.
# Material
Bone marrow samples obtained at diagnosis were used to prepare DNA from blastic cells, and this DNA was further enriched for sequences corresponding 2 regions of chromosome 9 and chromosome 22, respectively, using a capture probes.
Captured DNA was then subjected to illumina paired-end sequencing.
# Method 
##  QC analysis of the reads and Cleaning the fastq files using FastQC
Analyse the quality of reads inorder to make the decision whether to continue the analysis of the request resequencing due to low quality reads. I used AfterQC software (https://github.com/OpenGene/AfterQC.git)[1],that makes Automatic Filtering, Trimming, Error Removing and Quality Control statistics for fastq data,3 folders are generated(good,bad,QC), the final fastq files  that ll be used for the coming analysis are in the folder called 'good'.
The quality of the reads has been further investigeted using Fastq and multqc softwares in the Galaxy platforms revealing a good quality of long reads with low % of duplicates.
The workflow of genomics primary analysis will be the following:

## Galaxy platform workflow

The Jupyter,uranus data sets and the  human reference genome GRCh38/hg38 from the UCSC genome browser (http://hgdownload.soe.ucsc.edu/goldenPath/hg38/bigZips/hg38.fa.gz) has been uploaded to Galaxy from a local machine. After QC'ing we move on to map the reads using BWAmem.To check the genomic coverage of both mapped reads data sets, we used the Genome coverage tool from Bedtools software. Next using Picard's MarkDuplicates tool we process output of bwa mem. This step produces two files, a collection of deduplicated BAMs and a collection of duplicate metrics data produced by MarkDuplicates tool. Reads reallignment around Indels using Realigner target creator and Indel realigner tools from GATK is a crucial step to minimize alignment errors of reads ends.Finally the output Bam file from Indel realigner were used as input to Filter SAM or BAM tool to retain only properly mapped reads with mapping quality above 20 and mapping only to chr9 and chr22. Finally output of the filtering step is merged with MergeSAM tool and displayed in the UCSC Genome Browser.
To end up we need to visualize genomic coverage and the corresponding genes using the UCSC genome browser using the table browser tool in the UCSC,by selection one isoforme for each gene from GENCODE v29 database.

## 3- Somatic Structural aberation
Structural aberation is the alteration of the structure of the chromosome(sequence of genes or kind of genes in chromosome or no. of genes)
FusorSV (https://github.com/TheJacksonLaboratory/SVE),  is a tool that use different SV caller to obtain more accurate results[2]

