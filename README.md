Somatic mutations: variants only accumulate in the tumor or the affected tissues, and they are not transmitted from generation to generation. 
In this project we ll be analysing DNA sequences,corresponding of regions of Chromosome (9,22) from 2 patients that were diagnosed with ALL-B. 
Aim of work :the search for specific alterations that may be susceptible to tailored therapy
Identify mutations caused by somatic variation that may affect proteins functions
The workflow of genomics analysis will be the following:
## 1- Cleaning the fastq files.
Trimming the adapters and cleaning low quality bases in the reads using Trimmomatic is necessary to have better alignment in the further steps
Tool: Trimmomatic
## 2- QC analysis of the reads
Analyse the quality of reads inorder to make the decision whether to continue the analysis of the request resequencing due to low quality reads.
## 3- Reads mapping to the reference human genome
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
MetaSV is a tool that use different SV caller to obtain more accurate results
