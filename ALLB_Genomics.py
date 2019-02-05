
import os
import sys
import csv
from itertools import product
import fileinput
from datetime import datetime
start_time = datetime.now()
new=open("/Users/sondeskalboussi/Desktop/Time1.txt","w")
class BSMTBGan(object):
    def __init__(self,global_Dir,user_ref,user_ref_fasta,File_table,mappers):
      
      #Settings parameters and paths
        self.global_Dir=global_Dir# path to home directory in the cluster or local machine
        self.results=self.global_Dir+"/Results/"#results folder
        self.tools="/Users/sondeskalboussi/Tools"#Tools binaries
        self.ref_gen_Dir="/Users/sondeskalboussi/Desktop/ALL-b_pipeline/Reference/"#Hg38 reference
        self.genome=user_ref
        self.fasta=user_ref_fasta
        self.ref_genome="/Users/sondeskalboussi/Desktop/ALL-b_pipeline/Reference/"+self.genome+"/FASTA/"#fasta file of reference genome
        self.ref="/Users/sondeskalboussi/Desktop/ALL-b_pipeline/Reference/"+self.genome+"/FASTA/"+self.fasta#fasta file of reference genome
        self.bwa=self.tools+"/bwa-0.7.15/"#path to bwa mem software
        self.trimmomatic=self.tools+"/Trimmomatic-0.36/"
        self.GATK=self.tools+"/GenomeAnalysisTK.jar"
        self.samtools="/samtools/"
        self.picard=self.tools+"/picard.jar"
        self.bedtools="bedtools/"
        self.bcftools=self.tools+"/bcftools-1.4/"
        self.mappers=mappers
        self.dbSNP="dbSNP/dbSNP.vcf"
        self.illumina_adapters="illumina_adapters.fna.fasta"
        self.fastQFileData=File_table
        self.Sample=[]
        self.libraryPE={}
        self.librarySE={}
        self.valide_bam=[]#bam files validated by statistics for variant call
        self.Delly=self.tools+"/delly/src/delly"
        self.Normal='path to controle sample'

    print 
    print
    print
    
    #Process the reference genome for the later use for mapping and GATK_
    def process_ref_genome(self):
           print ("""
                ============================================================================
                     Reference genomes processing for the annotation is starting"!
                ============================================================================
            """)
           for fl in os.listdir(self.ref_genome):
               if fl.endswith(".fa") or fl.endswith(".fasta"):
                        os.chdir(self.ref_genome)
                        ID=os.path.splitext(fl)[0]
                        os.system("{}bwa index {}""".format(self.bwa,fl))#index the reference fasta for bwa
                        os.system("{}novoindex {}.nix {} ".format(self.novocraft,ID,fl))#index the reference fasta for novoalign
                        os.system("samtools faidx {}""".format(fl))
                        os.system("java -jar {} CreateSequenceDictionary R={} O={}.dict""".format(self.picard,fl,ID))#GATK
           return "Processed genome files are ready !"
   
        print ("""
                  ====================================
                       The trimming starts now!
                  ====================================
        """)

        if not os.path.isdir(self.results+"Trimming"):
               os.makedirs(self.results+"Trimming")
        for row in file:
                SM=row[3]
                if SM not in self.Sample:
                    self.Sample.append(SM)
                ID=row[1]
                PL=row[2]
                LB=row[4]
                if "R1" or "R2" in row[0] and SM in row[0] :
                    read=row[0]
                    if "R1" in read:
                        read1=read
                        LB1=row[4]
                    if "R2"in read:
                        read2=read
                        LB2=row[4]
                        if LB1==LB2:
                            output1=self.results+"Trimming/"+SM+"_"+ID+"_"+LB+"_"+PL+"_R1.PE_paired_trimed.fq"
                            output2=self.results+"Trimming/"+SM+"_"+ID+"_"+LB+"_"+PL+"_R2.PE_paired_trimed.fq"
                            os.system("""java -jar {}trimmomatic-0.36.jar PE -threads 4 -phred33 {} {} {} {} {} {} ILLUMINACLIP:TruSeq3-PE:2:30:10 LEADING:3 TRAILING:3 MINLEN:36""".format(self.trimmomatic,read1,read2,output1,output2))#change extension of trimmomatic for the cluster
    # Create a dictionary where the keys are samples ID and value is the reads ID when they are coming from more than sequencing library so later you can merge all the output for the same sample 
        for fl in os.listdir(self.results+"Trimming/"):
                os.chdir(self.results+"Trimming/")
                if fl.endswith("R1.PE_paired_trimed.fq"):
                                    IDPE=fl.split("R1.PE_paired_trimed.fq")[0]
                                    sm=IDPE.split("_")[0]
                                    lib=IDPE.split("_")[2]
                                    if sm not in self.libraryPE.keys():
                                        self.libraryPE[sm]=[lib]
                                    else:
                                        self.libraryPE[sm].append(lib)
               
        return "Trimming is done!"
   
   # Mapping reads using BWA-mem, mappers is a list of software that we want to use 
    def Mapping(self):
        print ("""
              ========================
                Mapping starts now!
              ========================
            """)
        for fl in os.listdir(self.results+"Trimming/"):
          os.chdir(self.results+"Trimming/")
          if ".fq" in fl:
            SM=fl.split("_")[0]
            ID=fl.split("_")[1]
            PL=fl.split("_")[3]
            LB=fl.split("_")[2]
            readGroup="@RG\\tID:"+ID+"\\tSM:"+SM+"\\tLB:"+LB+"\\tPL:"+PL
            input=os.path.join(self.results+"Trimming/",fl)
            for map in self.mappers:# mappers is a list of mapping software name (BWA-mem)
                    if not os.path.isdir(self.results+map+"/Alignment/Final_Bam"):
                        os.makedirs(self.results+map+"/Alignment/Final_Bam")         
                    #Mapping the PE reads
                    if "PE_paired_trimed" in input:
                        if "R1" in input:
                            R1=input
                        if "R2" in input:
                            R2=input
                            output_map=self.results+map+"/Alignment/"+SM+"_"+LB+"_PE_"+map+".sam"
                            if map=="BWA":
                                os.system("""{}bwa mem -M -t 4 -R "{}" {} {} {} > {}""".format(self.bwa,readGroup,self.ref,R1,R2,output_map))
                                os.system("""samtools view -Sb -o {} {}""".format(self.results+map+"/Alignment/"+SM+"_"+LB+"_PE_"+map+".bam",output_map))
        for map in self.mappers:
                for fl in os.listdir(self.results+map+"/Alignment"):
                    os.chdir(self.results+map+"/Alignment")
                    if fl.endswith(".sam"):
                            os.system("""JAVA -jar {} SortSam INPUT={} OUTPUT={} SORT_ORDER=coordinate""".format(self.picard,fl,fl.replace(".sam","_sorted.bam")))
                            os.system("""JAVA -jar {} BuildBamIndex INPUT={}""".format(self.picard,fl.replace(".sam","_sorted.bam")))
        return "Mapping is done"
    #Mark and remove the PCR duplicates
    def PCR_dup_mark(self):
           
            print ("""
                            ==========================================
                                PCR duplication marking begins!
                            ==========================================
                          """)
            for map in self.mappers: 
                for fl in os.listdir(self.results+map+"/Alignment/"):
                    os.chdir(self.results+map+"/Alignment/")
                    if fl.endswith(map+"_sorted.bam"):
                            input=os.path.join(self.results+map+"/Alignment/",fl)
                            output=os.path.join(self.results+map+"/Alignment/",fl.replace(map+"_sorted.bam",map+"_dedup.bam"))
                            metric=os.path.join(self.results+map+"/Alignment/",fl.replace(map+"_sorted.bam",map+"_dedup.bam.metrics"))
                            os.system("""java -jar {} MarkDuplicates INPUT={} OUTPUT={} METRICS_FILE={}""".format(self.picard,input,output,metric))#add{}
                            os.system("""java -jar {} BuildBamIndex INPUT={}""".format(self.picard,output))
            return "PCR duplicates removal is done!"
                        
   # Realign around the indels using GATK the final output is a sorted indexed bam file
    def realignment(self):                 
        print ("""
                =====================================================
                      Realignment around the indels using GATK
                =====================================================
            """)
        for map in self.mappers:           
            for fl in os.listdir(self.results+map+"/Alignment/"):
                  os.chdir(self.results+map+"/Alignment/")
                  if fl.endswith("_dedup.bam"):
                        input=os.path.join(self.results+map+"/Alignment/",fl)
                        output1=os.path.join(self.results+map+"/Alignment/",fl.replace("_dedup.bam",".intervals"))
                        output2=os.path.join(self.results+map+"/Alignment/",fl.replace("_dedup.bam","_realg.bam"))
                        output3=os.path.join(self.results+map+"/Alignment/",fl.replace("_dedup.bam","_realg_sorted.bam"))
                        os.system("""java -jar {} -T RealignerTargetCreator -nt 4 -R {} -I {} -o {}""".format(self.GATK,self.ref,input,output1))#change thread
                        os.system("""java -jar {} -T IndelRealigner -R {} -I {} -targetIntervals {} -o {}""".format(self.GATK,self.ref,input,output1,output2))
                        os.system("""JAVA -jar {} SortSam INPUT={} OUTPUT={} SORT_ORDER=coordinate""".format(self.picard,output2,output3))
                        #os.system("""samtools sort -o {} {}""".format(output3,output2))
                        os.system("""java -jar {} BuildBamIndex INPUT={}""".format(self.picard,output3))
        return "realignment around indels is done!"
    #Base quality score recalibration
    def base_qual_recal(self):                   
        for map in self.mappers:
             print ("""
                =================================================================
                    {} alignment base quality score recalibration starts now!
                =================================================================
             """.format(map))                     
             for fl in os.listdir(self.results+map+"/Alignment/"):
                  os.chdir(self.results+map+"/Alignment/")
                  if fl.endswith(map+"_realg_sorted.bam"):
                        input=os.path.join(self.results+map+"/Alignment/",fl)
                        output1=os.path.join(self.results+map+"/Alignment/",fl.replace("_realg_sorted.bam","_recal_data.table"))
                        output2=os.path.join(self.results+map+"/Alignment/",fl.replace("_realg_sorted.bam","_recal.bam"))
                        output3=os.path.join(self.results+map+"/Alignment/",fl.replace("_realg_sorted.bam","_recal_sorted.bam"))
                        os.system("""java -jar {} -T BaseRecalibrator -nct 4 -R {} -I {} -knownSites {} -o {}""".format(self.GATK,self.ref,input,self.dbSNP,output1))#add path to GATK cluster
                        os.system("""java -jar {} -T PrintReads -nct 4 -R {} -I {} -BQSR {}  -o {}""".format(self.GATK,self.ref,input,output1,output2))
                        #os.system("""samtools sort -o {} {}""".format(output3,output2))
                        os.system("""JAVA -jar {} SortSam INPUT={} OUTPUT={} SORT_ORDER=coordinate""".format(self.picard,output2,output3))
                        os.system("""java -jar {} BuildBamIndex INPUT={}""".format(self.picard,output3))

        return "Base recalibration is finished!"
       #Check if there is samples from more than one sequencing library if yes merge them and move them to Final_bam folder otherwise use bam file in Alignment directory with the extension final_mapper.bam
   
    def mapping_stat(self):
      
            print ("""
                 ==================================================================
                   Genome coverage and reads mappability statistics are starting !
                 ==================================================================
                """)
            for map in self.mappers:
                MAP={}
                COV={}
                if not os.path.isdir(self.results+map+"/statistics"):
                    os.makedirs(self.results+map+"/statistics")
                if len(os.listdir(self.results+map+"/Alignment/Final_Bam/"))== 0:
                    Final_Bam=self.results+map+"/Alignment/"
                else:
                    Final_Bam=self.results+map+"/Alignment/Final_Bam/"
                for fl in os.listdir(Final_Bam):
                        os.chdir(Final_Bam)
                        if fl.endswith("_recal_sorted.bam") or "_merged_" in fl:
                            input=os.path.join(Final_Bam,fl)
                            output1=self.results+map+"/statistics/"+fl.replace(".bam",".coverage")
                            output2=self.results+map+"/statistics/"+fl.replace(".bam",".flagstat")
                            output3=self.results+map+"/statistics/"+fl.replace(".bam",".stats")
                            os.system("""java -jar {} -T DepthOfCoverage -R {} -I {} -o {} --omitDepthOutputAtEachBase --omitIntervalStatistics --omitLocusTable""".format(self.GATK,self.ref,input,output1))
                            os.system("""samtools flagstat {} > {}""".format(input,output2))
                            os.system("""samtools stats {} > {}""".format(input,output3))
                for fl in os.listdir(self.results+map+"/statistics/"):
                   os.chdir(self.results+map+"/statistics/")
                   ID=fl.split("_")[0]
                   if fl.endswith(".flagstat"):
                       for line in open(fl).readlines():
                            line=line.strip()
                            if "mapped (" in line:
                                x=line.split("mapped")[1].split(":")[0].split("(")[1].split("%")[0]
                                y=float(x)
                                MAP[ID]=y
                   if fl.endswith("_summary"):
                       for line in open(fl).readlines()[1:2]:
                                line=line.strip()
                                x=line.split("\t")[2]
                                y=float(x)
                                COV[ID]=y
                new=open(self.results+map+"/statistics/mapping_stat.txt","w")
                new1=open(self.results+map+"/statistics/Failed_mapping_stat.txt","w")
                new.write("""{}{}{}{}{}{}""".format("Sample","\t","%Mapped_reads","\t","Genome_coverage mean", "\n"))
                new1.write("""{}{}{}{}{}{}""".format("Sample","\t","%Mapped_reads","\t","Genome_coverage mean", "\n"))
                for key in MAP.keys():
                    if key in COV.keys():
                        if MAP[key]>=90 and COV[key]>=40 :
                            if key not in self.valide_bam:
                                self.valide_bam.append(key)
                            new=open(self.results+map+"/statistics/mapping_stat.txt","a")
                            new.write("""{} {}  {}{}""".format(key,MAP[key],COV[key],"\n"))
                    
                        else:
                            new1=open(self.results+map+"/statistics/Failed_mapping_stat.txt","a")
                            new1.write("""{} {}  {}{}""".format(key,MAP[key],COV[key],"\n"))
                new.close()
                new1.close()
            return " Statistics are done!"
    # Joint variant call SNP and indels using GATK
    def joint_variant_calling(self):
            print ("""
                  ==================================================================
                    Joint variant calling (SNP and indels) using GATK is starting !
                  ==================================================================
            """)                        
            for map in self.mappers:
                if not os.path.isdir(self.results+map+"/Joint_Variants"):
                    os.makedirs(self.results+map+"/Joint_Variants")
                if len(os.listdir(self.results+map+"/Alignment/Final_Bam/"))== 0:
                    Final_Bam=self.results+map+"/Alignment/"
                else:
                    Final_Bam=self.results+map+"/Alignment/Final_Bam/"
                for i in self.valide_bam:
                    for fl in os.listdir(Final_Bam):
                        if i in fl:
                            if fl.endswith("_recal_sorted.bam") or in fl:
                                os.chdir(Final_Bam)
                                input=os.path.join(Final_Bam,fl)
                                output=os.path.join(self.results+map+"/Joint_variants/",fl.replace(".bam","_GATK_snps_indels.vcf"))
                                os.system("""java -jar {} -T HaplotypeCaller -R {} -I {} -stand_call_conf 30 -o {}""".format(self.GATK,self.ref,input,output))#add path gatk cluster
            return "Joint variant calling (SNP and indels) using GATK is done!"
    #non-model organism, so variant re calibration isn't possible
    def joint_variant_calling_hard_filtering(self):
         print ("""
            ==================================================================
                  Joint_variant_calling_hard_filtering is starting !
            ==================================================================
            """)
         for map in self.mappers:
             for fl in os.listdir(self.results+map+"/Joint_Variants"):
                 os.chdir(self.results+map+"/Joint_Variants")
                 id=fl.split("_")[0]
                 if fl.endswith(".vcf"):
                     input=os.path.join(self.results+map+"/Joint_Variants",fl)
                     output1=self.results+map+"/Joint_variants/"+id+"_raw_snps.vcf"
                     output11=self.results+map+"/Joint_variants/"+id+"_filtered_snps.vcf"
                     output2=self.results+map+"/Joint_variants/"+id+"_raw_indels.vcf"
                     output22=self.results+map+"/Joint_variants/"+id+"_filtered_indels.vcf"
                     #variant selection
                     os.system("""java -jar {} -T SelectVariants -R {} -V {} -selectType SNP -o {}""".format(self.GATK,self.ref,input,output1))
                     os.system("""java -jar {} -T SelectVariants -R {} -V {} -selectType INDEL -o {}""".format(self.GATK,self.ref,input,output2))
                     #variant filtration
                     os.system("""java -jar {} -T VariantFiltration -R {} -V {} --filterExpression "QD < 2.0 || FS > 60.0 || MQ < 40.0 || MQRankSum < -12.5 || ReadPosRankSum < -8.0" --filterName "my_snp_filter" -o {}""".format(self.GATK,self.ref,output1,output11))
                     os.system("""java -jar {} -T VariantFiltration -R {} -V {} --filterExpression "QD < 2.0 || FS > 200.0 || ReadPosRankSum < -20.0" --filterName "my_indel_filter" -o {}""".format(self.GATK,self.ref,output2,output22))
         return " variant_calling_hard_filtering is done!"
    #number of SNPs identified with average quality scores and average mapping quality for filtered snp.csv and raw.vcf, parse the MQ mapping quality, QUAl field: quality score, wc - l nbre of snp then compute average
    def SNP_statistics(self):
        print ("""
            ==================================================================
                  SNP statistics is starting !
            ==================================================================
            """)
        for map in self.mappers:
            new=open(self.results+map+"/statistics/VCF_stat.txt","w")
            new.write("""{}{}{}{}{}{}{}{}""".format("Sample","\t","Nbre_SNP","\t","AV_MQ","\t","AV_QUAL","\n"))
            for fl in os.listdir(self.results+map+"/Joint_variants/"):
                os.chdir(self.results+map+"/Joint_variants/")
                f=os.path.join(self.results+map+"/Joint_variants/",fl)
                if f.endswith("_filtered_snps.vcf"):
                    Nbre_SNP=0
                    AV_MQ=0
                    AV_QUAL=0
                    for lines in filter(None, open(f).readlines()):
                        if not lines.startswith("##") and not lines.startswith("#") and "my_snp_filter" not in lines:
                            Nbre_SNP+=1
                            line=lines.split("\t")
                            AV_QUAL+=float(line[5])/Nbre_SNP
                            L=line[7].split(";")[8]
                            if "MQ" in L:
                                AV_MQ+=float(line[7].split(";")[8].split("=")[1])/Nbre_SNP
                    new.write("""{}{}{}{}{}{}{}{}""".format(fl.split("_")[0],"\t",Nbre_SNP,"\t",AV_MQ,"\t",AV_QUAL,"\n"))
            new.close()
        return "VCF statistics is done"
    def Genotype_structural_variation_calling(self):
        print ("""
                  ======================================================
                    Somatic structural_variation_calling is starting !
                  ======================================================
        """)                    
        for map in self.mappers:
            if not os.path.exists(self.results+map+"/Struc_Variants"):
                    os.makedirs(self.results+map+"/Struc_Variants")
            if len(os.listdir(self.results+map+"/Alignment/Final_Bam/"))== 0:
                Final_Bam=self.results+map+"/Alignment/"
            else:
                Final_Bam=self.results+map+"/Alignment/Final_Bam/"
            id=fl.split("_")[0]
                 if fl.endswith(".vcf"):
                     input=os.path.join(self.results+map+"/Joint_Variants",fl)
                     output1=self.results+map+"/Joint_variants/"+id+"_raw_snps.vcf"    
            for id in self.valide_bam:
                for fl in os.listdir(Final_Bam):
                    if id in fl:
                        if fl.endswith("_recal_sorted.bam") in fl and id in os.listdir(self.Normal):
                            os.chdir(Final_Bam)
                            input1=os.path.join(Final_Bam,fl)
                            input2=os.path.join(self.Normal,fl)
                            output=self.results+map+"/Struc_Variants/"+fl.replace("_realg_sorted.bam","_SV"                          os.system("""delly call -x hg38.excl -o {}.bcf -g {} {} {}""".format(id,self.ref,output,input1,input2))#add delly
                            os.system("""delly call -x hg38.excl -o {}.bcf -g {} {} {}""".format(id,self.ref,output,input1,input2))#somatic sv call
                            os.system("""delly filter -f somatic -o {}.pre.bcf -s {}.tsv {}.bcf""".format(id,output,input1,id))#somatic sv prefiltering
                            os.system(""" delly call -g hg38.fa -v {}.pre.bcf -o geno.bcf -x hg38.excl {}.bam {}.bam ... controlN.bam""".format(id,input1,input2,id))#Genotype pre-filtered somatic sites
                            os.system("""delly filter -f somatic -o {}.somatic.bcf -s {}.tsv geno.bcf""".format(id,input1,id))#Post-filter for somatic SVs using all control samples.
                            os.system("""bcftools view {}.bcf > {}.vcf """.format(output,output))#convert bcf to vcf
        return "GSV is done!"        
                 delly call -x hg19.excl -o t1.bcf -g hg19.fa tumor1.bam control1.bam


def main():
        pipeline=BSMTBGan("path to Pipeline","path to HG38.fa","path to Fastq.csv", mappers=["BWA"])
        try:
            #print pipeline.process_ref_genome()
            #print pipeline.trimming()
            #print pipeline.Mapping()
            #print pipeline.PCR_dup_mark()
            #print pipeline.realignment()
            #print pipeline.base_qual_recal()
            #print pipeline.mapping_stat()
            #print pipeline.joint_variant_calling()
            #print pipeline.joint_variant_calling_hard_filtering()
            #print pipeline.SNP_statistics()
            #print pipeline.Genotype_structural_variation_calling()

        except IOError as e:
            print("I/O error: {0}".format(e))    
    
if __name__ == '__main__':main()
end_time = datetime.now()
new.write('BWA pipeline execution time: {}'.format(end_time - start_time))
new.close()

