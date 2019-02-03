
Childhood acute lymphoblastic leukemia (ALL) is a type of cancer in which the bone marrow makes too many immature lymphocytes caused by Somatic mutations which are genetic variants that only accumulate in the tumor or the affected tissues, and they are not transmitted from generation to generation. 
# Aim of mini-project
In this project we ll be analysing DNA sequences,corresponding of regions of Chromosome (9,22) from 2 patients that were diagnosed with ALL-B inorder to identify mutations caused by somatic variation that may affect proteins functions or specific alterations that may be susceptible to tailored therapy.
# Material
Bone marrow samples obtained at diagnosis were used to prepare DNA from blastic cells, and this DNA was further enriched for sequences corresponding 2 regions of chromosome 9 and chromosome 22, respectively, using a capture probes.
Captured DNA was then subjected to illumina paired-end sequencing.
# Method 
##  QC analysis of the reads and Cleaning the fastq files
Analyse the quality of reads inorder to make the decision whether to continue the analysis of the request resequencing due to low quality reads. I used AfterQC software (https://github.com/OpenGene/AfterQC.git)[1],that makes Automatic Filtering, Trimming, Error Removing and Quality Control statistics for fastq data,3 folders are generated(good,bad,QC), the final fastq files  that ll be used for the coming analysis are in the folder called 'good'.
The quality of the reads has been further investigeted using Fastq and multqc softwares in the Galaxy platforms revealing a good quality of long reads with low % of duplicates.
The workflow of genomics primary analysis will be the following:

# Galaxy platform workflow
The Jupyter,uranus data sets and the  human reference genome GRCh38/hg38 from the UCSC genome browser (http://hgdownload.soe.ucsc.edu/goldenPath/hg38/bigZips/hg38.fa.gz) has been uploaded to Galaxy from a local machine, inorder to speed up the analysis so we run once all the further analysis for both data sets, the data set has been collected into a paired collection (a list) that groups the paired reads of each patients under one file. After QC'ing we move on to map the reads using BWAmem, process the resulting BAM datasets, and visualize coverage in a genome browser.
Mapping this collection to human genome with bwa mem produces a flat collection of BAM datasets.it s very important to set when readgroups parameter in this step. This allows us to merge individual BAM datasets into one at the end of this analysis. Next using Picard's MarkDuplicates tool we process output of bwa mem. This step produces two files, a collection of deduplicated BAMs and a collection of duplicate metrics data produced by MarkDuplicates tool. We then filter BAM collection produced by MarkDuplicates using Filter SAM or BAM tool to retain only properly mapped reads with mapping quality above 20 and mapping only to chr9 and chr22. Finally output of the filtering step is merged with MergeSAM tool and displayed in the UCSC Genome Browser. Again, merging is only possible because we have set the readgroups during the mapping step.


The Jupyter,uranus data sets and the  human reference genome GRCh38/hg38 from the UCSC genome browser (http://hgdownload.soe.ucsc.edu/goldenPath/hg38/bigZips/hg38.fa.gz) has been uploaded to Galaxy from a local machine. After QC'ing we move on to map the reads using BWAmem. Next using Picard's MarkDuplicates tool we process output of bwa mem. This step produces two files, a collection of deduplicated BAMs and a collection of duplicate metrics data produced by MarkDuplicates tool. Reads reallignment around Indels using Realigner target creator and Indel realigner tools from GATK is a crucial step to minimize alignment errors of reads ends
Base quality score recalibration/recalculation based on the alignment data using GATK.Finally,we filter the recalibrated BAM file


using Filter SAM or BAM tool to retain only properly mapped reads with mapping quality above 20 and mapping only to chr9 and chr22. Finally output of the filtering step is merged with MergeSAM tool and displayed in the UCSC Genome Browser. Again, merging is only possible because we have set the readgroups during the mapping step.


We then filter BAM collection produced by MarkDuplicates using Filter SAM or BAM tool to retain only properly mapped reads with mapping quality above 20 and mapping only to chr9 and chr22.




Finally output of the filtering step is merged with MergeSAM tool and displayed in the UCSC Genome Browser. Again, merging is only possible because we have set the readgroups during the mapping step.

To end up we need to visualize genomic coverage and the corresponding genes using the UCSC genome browser.







# Python software workflow


## 2- Reads mapping to the reference human genome
The reads are now aligned to the human reference genome, using BWA-mem, inorder to identify the covered genomic region and the corresponding genes, further steps are needed to remove any putative errors in the mapping results. The steps are the following:
        1- Read pairs with identical start and end positions were assumed to be PCR duplicates and were removed using Picard
        2- Reads reallignment around Indels to minimize alignment errors of reads ends using GATK
        3- Base quality score recalibration/recalculation based on the alignment data using GATK
## 3- Somatic mutations call (single nucleotide variant ,indels,structural variation )
SNP and genotype identification
Comparison of the sequenced reads to their point of alignment on the human genome providing a list of genomic variations that is organized according to their genomic location (chromosome and position) and the variant allele. These variants are knwon as somatic variants


Single nucleotide variants (SNVs) were called with SomaticSniper.
Structural variations were called by Pindel
Short indels where called by Delly
Copy number variation called by VarScan2
After identifying the variants, a further step is needed in order to filter false positive variants,using GATK, to investigate the ability of a read containing a variant allele to map uniquely to a single location on the human reference genome.
Variants were filtered based on sequence coverage and quality scores
## 4- Annotation
All variants were annotated against the Ensembldatabase.
# Structural aberation
Structural aberation is the alteration of the structure of the chromosome(sequence of genes or kind of genes in chromosome or no. of genes)
FusorSV (https://github.com/TheJacksonLaboratory/SVE),  is a tool that use different SV caller to obtain more accurate results[2]
Genomic region
chr9:1-138,394,717
chr22:1-50,818,468
