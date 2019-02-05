# ALL_B-Project
Genomic analysis of two DNA samples from acute lymphoblastic leukemia patients
Acute lymphoblastic leukemia (ALL) is a type of cancer in which the bone marrow makes too many immature lymphocytes caused by Somatic mutations which are genetic variants that only accumulate in the tumor or the affected tissues, and they are not transmitted from generation to generation.

# Aim of mini-project
In this project we ll be analysing DNA sequences,corresponding of regions of Chromosome (9,22) from 2 patients that were diagnosed with ALL-B inorder to identify mutations genetic aberations.

# Material
Bone marrow samples obtained at diagnosis were used to prepare DNA from blastic cells, and this DNA was further enriched for sequences corresponding 2 regions of chromosome 9 and chromosome 22, respectively, using a capture probes. Captured DNA was then subjected to illumina paired-end sequencing.

# Method
##QC analysis of the reads and Cleaning the fastq files using FastQC
Analyse the quality of reads inorder to make the decision whether to continue the analysis of the request resequencing due to low quality reads. I used AfterQC software (https://github.com/OpenGene/AfterQC.git),that makes Automatic Filtering, Trimming, Error Removing and Quality Control statistics for fastq data,3 folders are generated(good,bad,QC), the final fastq files that ll be used for the coming analysis are in the folder called 'good'. The quality of the reads has been further investigeted using Fastq and multqc softwares in the Galaxy platforms revealing a good quality of long reads with low % of duplicates. We have medium length reads that are less than 200 pb with almost good sequencing quality with small duplication percentage. 

The workflow of genomics primary analysis will be the following:

## Galaxy platform workflow
The Jupyter,uranus data sets and the human reference genome GRCh38/hg38 from the UCSC genome browser (http://hgdownload.soe.ucsc.edu/goldenPath/hg38/bigZips/hg38.fa.gz) has been uploaded to Galaxy from a local machine. After QC'ing we move on to map the reads using BWAmem.To check the genomic coverage of both mapped reads data sets, we used the Genome coverage tool from Bedtools software. Next using Picard's MarkDuplicates tool we process output of bwa mem. This step produces two files, a collection of deduplicated BAMs and a collection of duplicate metrics data produced by MarkDuplicates tool. Reads reallignment around Indels using Realigner target creator and Indel realigner tools from GATK is a crucial step to minimize alignment errors of reads ends.Finally the output Bam file from Indel realigner were used as input to Filter SAM or BAM tool to retain only properly mapped reads with mapping quality above 20 and mapping only to chr9 and chr22. Finally output of the filtering step is merged with MergeSAM tool and displayed in the UCSC Genome Browser. To end up we need to visualize genomic coverage and the corresponding genes using the UCSC genome browser using the table browser tool in the UCSC,by selection one isoforme for each gene from GENCODE v29 database obtaining list of canonical genes. 
For joint variant call, i used Freebayes a variant detector designed to find small polymorphisms,SNPs, indels ,MNPs (multi-nucleotide polymorphisms), and complex events (composite insertion and substitution events) smaller than the length of a short-read sequencing alignment, the output is a vcf file withreference and alteration
## Genomic covered region and corresponding genes 
chr9:12133-138286430 with 6945 one isoforme genes
chr22:10736170-50783662 with 4615 one isoforme genes
List of genes are stored in two files corresponding to each chromosome (please check canonical genes folder for the list of genes for details.
## Structural aberation
Structural aberation is the alteration of the structure of the chromosome(sequence of genes or kind of genes in chromosome or no. of genes),it Can be divided into insertions, deletions, inversions, and translocations (either inter or intra-chromosomal). From the literature, chr22 or Ph1 is the most frequent abberation in ALLB tumor 
From the Genome Browser we visualized Database of Genomic Variants: Structural Variation parameters and CNVs and Indels where displayed either in bleu there is a gain in size relative to the reference or in red if there is a loss in size relative to the reference.

FusorSV (https://github.com/TheJacksonLaboratory/SVE), is a tool that use different SV caller to obtain more accurate results


## Comments

If we want to extend the analysis protocole for bigger data set of patients, it might be agood idea to have longer sequencing reads so we can obtain better mapping matches and this ll influence complex rearrengements and oncogenes amplifications search so i might suggest to use the pacbio platform for sequencing than benchmark the results of mapping and SV call.
Structural aberation at the level of  chr9 and 22(gene translocation from chr9 to chr22) is a marker of tumor will certainly affect the protein since this mutant gene is differentially expressed in ALLB patients, this might be a confirmation of the clinical diagnosis but we can go larger and sequence the whole genome of patients so we can creat kind of map of all the differentially expressed genes to have a clearer idea about specific genetic markers and early genetic diagnosis of the tumor. This might open new doors toward targeted genetic treatments and help the patients to avoid putative side effects of non necessary therapy.
What i suggest also is providing normal or reference samples with the tumor samples, this is very important and facilitate the analysis since the major part of SV tools are based on comparing normal vs tumor. I have tried to use different tools (Dellu,Pindel,Mutec2 and others) but couldnt go further. 
It very important the filter out the germline variants to identify somatic one so using public database. 
Filtering out the germeline variation can provide better identification of somatic ones. 
