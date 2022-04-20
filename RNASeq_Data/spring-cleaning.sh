#!/bin/bash

#Stephen T. and Ally S. 
#S4B
#Spring 2022

######################################################################################
############################# RNA-seq Data Cleanup ###################################
######################################################################################

#RNA-seq data analysis to compare gene expression of two treatment groups of mice

#In this script we will:
	#Check the quality of the reads
	#Clean reads
	#Index the genome
	#Map the reads to a reference genome
	#Get read counts per gene
	#Use the read counts to find differentially expressed genes (DEGs)

#Packages: FastQC, TrimGalore, STAR, 

######################################################################################
######################################################################################
######################################################################################

#1) Check the quality of the reads

function quality_check {

	#input files: *.fastq.gz in ~/RNASeq_Data/Case and ../Control - hardcoded here, but this can be customized by the user.
	#output files: *.html to be viewed in a web browser
		#This file will be located in ~/RNASeq_Data/FastQC as specified in the code. This can be modified by the user if they choose.
	#Packages: FastQC

	###############################################################################

	module load fastqc/0.11.9 #loading the fastqc module from ASC

	mkdir FastQC #make a new directory for output files
	cd Case #move into the directory with the fastq files
	fastqc *.fastq.gz -o ~/s4b-project/RNASeq_Data/FastQC #Use FastQC to perform a quality check of sequences
	cd ../Control #move into directory with the control fastq files
	fastqc *.fastq.gz -o ~/s4b-project/RNASeq_Data/FastQC #Use FastQC to perform a quality check of sequences

}

#2) Clean the reads by trimming

function trim_reads {

	#input files: *.fastq.gz in ~/s4b-project/RNASeq_Data/Case and ../Control - hardcoded here, but this can be customized by the user.
        #output files: *val_1.fq.gz (for R1) or *val_2.fq.gz (for R2) files 
		#These files will be located in ~/s4b-project/RNASeq_Data/TrimmedReads
	#packages: trimgalore
		#trimgalore will automatically detect and cut sequences at illumina adapters

	###############################################################################

	module load trimgalore/0.6.6 #loading trimgalore module from ASC

	#making new directories in ~/s4b-project/RNASeq_Data for the trimmed data to be stored
	mkdir TrimmedReads #create new
	mkdir TrimmedReads/Case
	mkdir TrimmedReads/Control

	cd Case #move into the directory with the fastq files
		#TrimGalore was struggling to pair files, so this had to be hard coded for each pair
        trim_galore --paired --output_dir ~/s4b-project/RNASeq_Data/TrimmedReads/Case 4040-KH-14.4040-KH-14_0_filtered_R1.fastq.gz 4040-KH-14.4040-KH-14_0_filtered_R2.fastq.gz
	trim_galore --paired --output_dir ~/s4b-project/RNASeq_Data/TrimmedReads/Case 4040-KH-16.4040-KH-16_0_filtered_R1.fastq.gz 4040-KH-16.4040-KH-16_0_filtered_R2.fastq.gz
	trim_galore --paired --output_dir ~/s4b-project/RNASeq_Data/TrimmedReads/Case 4040-KH-21.4040-KH-21_0_filtered_R1.fastq.gz 4040-KH-21.4040-KH-21_0_filtered_R2.fastq.gz
	trim_galore --paired --output_dir ~/s4b-project/RNASeq_Data/TrimmedReads/Case 4040-KH-22.4040-KH-22_0_filtered_R1.fastq.gz 4040-KH-22.4040-KH-22_0_filtered_R2.fastq.gz
	trim_galore --paired --output_dir ~/s4b-project/RNASeq_Data/TrimmedReads/Case 4040-KH-23.4040-KH-23_0_filtered_R1.fastq.gz 4040-KH-23.4040-KH-23_0_filtered_R2.fastq.gz
	trim_galore --paired --output_dir ~/s4b-project/RNASeq_Data/TrimmedReads/Case 4040-KH-24.4040-KH-24_0_filtered_R1.fastq.gz 4040-KH-24.4040-KH-24_0_filtered_R2.fastq.gz
	trim_galore --paired --output_dir ~/s4b-project/RNASeq_Data/TrimmedReads/Case 4040-KH-25.4040-KH-25_0_filtered_R1.fastq.gz 4040-KH-25.4040-KH-25_0_filtered_R2.fastq.gz

        cd ../Control #move into directory with the control fastq files
	trim_galore --paired --output_dir ~/s4b-project/RNASeq_Data/TrimmedReads/Control 4040-KH-1.4040-KH-1_0_filtered_R1.fastq.gz 4040-KH-1.4040-KH-1_0_filtered_R2.fastq.gz
	trim_galore --paired --output_dir ~/s4b-project/RNASeq_Data/TrimmedReads/Control 4040-KH-4.4040-KH-4_0_filtered_R1.fastq.gz 4040-KH-4.4040-KH-4_0_filtered_R2.fastq.gz
	trim_galore --paired --output_dir ~/s4b-project/RNASeq_Data/TrimmedReads/Control 4040-KH-5.4040-KH-5_0_filtered_R1.fastq.gz 4040-KH-5.4040-KH-5_0_filtered_R2.fastq.gz
	trim_galore --paired --output_dir ~/s4b-project/RNASeq_Data/TrimmedReads/Control 4040-KH-6.4040-KH-6_0_filtered_R1.fastq.gz 4040-KH-6.4040-KH-6_0_filtered_R2.fastq.gz
	trim_galore --paired --output_dir ~/s4b-project/RNASeq_Data/TrimmedReads/Control 4040-KH-17.4040-KH-17_0_filtered_R1.fastq.gz 4040-KH-17.4040-KH-17_0_filtered_R2.fastq.gz
	trim_galore --paired --output_dir ~/s4b-project/RNASeq_Data/TrimmedReads/Control 4040-KH-18.4040-KH-18_0_filtered_R1.fastq.gz 4040-KH-18.4040-KH-18_0_filtered_R2.fastq.gz

}

#3) Check the quality of the trimmed reads

function qc_trimmed {

        #input files: *.fastq.gz in ~/s4b-project/RNASeq_Data/TrimmedReads/Case and ../Control - hardcoded here, but this can be customized by the user.
        #output files: *.html to be viewed in a web browser
                #This file will be located in ~/s4b-project/RNASeq_Data/TrimmedReads/trimmed-FastQC as specified in the code. This can be modified by the user if they need. 
        #Packages: FastQC

        ###############################################################################

        module load fastqc/0.11.9 #loading the fastqc module from ASC

        mkdir TrimmedReads/trimmed-FastQC #make a new directory for output files
        
	cd TrimmedReads/Case #move into the directory with the fastq files
        fastqc *.fq.gz -o ~/s4b-project/RNASeq_Data/TrimmedReads/trimmed-FastQC #Use FastQC to perform a quality check of sequences
        cd ../Control #move into directory with the control fastq files
        fastqc *.fq.gz -o ~/s4b-project/RNASeq_Data/TrimmedReads/trimmed-FastQC #Use FastQC to perform a quality check of sequences

}


#4) Index the genome
#5) Align the reads to the genome
#6) Count reads

function mapping {

	#input files: sequences to be mapped - *.fq.gz in ~/s4b-project/RNASeq_Data/TrimmedReads/Case and ../Control - hardcoded here, but this can be customized by the user of this script. 
		# genome to be indexed and mapped to - ./ncbi-genomes-2022-04-19/GCF_000001635.27_GRCm39_genomic.fna
        #output files: 
		#Indexing Output: Log files that are useful for quality checking and debugging, SJ.out.tab (splice junctions), Genome, SAindexes, chrLengths
		#Alignment Output: Aligned.out.sam file with mapped sequences will be in RNASeq_Data directory along with Log files to document the quality of mapping (Log files are very useful for quality control and debugging) 
        #Packages: STAR (Spliced Transcripts Alignment to a Reference)
		#requires 8 CPU cores and 30 gb memory to run

        ###############################################################################

	source /opt/asn/etc/asn-bash-profiles-special/modules.sh
	module load star/2.7.0e 

	###############################INDEXING GENOME#################################

	#STAR --runThreadN ___ \\ number of cores
	#--runMode genomeGenerate \\ genome mode
	#--genomeDir ___ \\ path to the directory for output
	#--genomeFastaFiles ___ \\ path to the FASTA files of the genome

	###############################################################################
	
	#STAR --runThreadN 8 --runMode genomeGenerate --genomeDir /home/aubars001/s4b-project/RNASeq_Data/ncbi-genomes-2022-04-19 --genomeFastaFiles /home/aubars001/s4b-project/RNASeq_Data/ncbi-genomes-2022-04-19/GCF_000001635.27_GRCm39_genomic.fna
	
	###########################MAPPING RNA-SEQ TO INDEX#############################

	#STAR --runThreadN ___ \\ number of cores
        #--genomeDir ___ \\ path to the directory containing indexed genome
        #--readFilesCommand zcat \\ tells STAR that files are gunzipped (.gz)
	#--readFilesIn ___,___,___ \\ files to be mapped to genome separated by commas
	#--quantMode GeneCounts \\ tells STAR to count number reads per gene while mapping (a read is counted if it overlaps [1nt or more] one an donly one gene)

        ###############################################################################

	STAR --runThreadN 8 --quantMode GeneCounts --genomeDir /home/aubars001/s4b-project/RNASeq_Data/ncbi-genomes-2022-04-19 --readFilesCommand zcat --readFilesIn /home/aubars001/s4b-project/RNASeq_Data/TrimmedReads/Case/4040-KH-21.4040-KH-21_0_filtered_R1_val_1.fq.gz,/home/aubars001/s4b-project/RNASeq_Data/TrimmedReads/Case/4040-KH-21.4040-KH-21_0_filtered_R2_val_2.fq.gz,/home/aubars001/s4b-project/RNASeq_Data/TrimmedReads/Control/4040-KH-18.4040-KH-18_0_filtered_R1_val_1.fq.gz,/home/aubars001/s4b-project/RNASeq_Data/TrimmedReads/Control/4040-KH-18.4040-KH-18_0_filtered_R2_val_2.fq.gz

}


function main {

	#quality_check /home/aubars001/s4b-project/RNASeq_Data
	#trim_reads /home/aubars001/s4b-project/RNASeq_Data
	#qc_trimmed /home/aubars001/s4b-project/RNASeq_Data
	mapping /home/aubars001/s4b-project/RNASeq_Data

}

#############################################################################################################################
################################################## CODING REGION BELOW ######################################################
#############################################################################################################################

echo "Time to get our seq on!"

main

echo "We just got seq-y. CongRATS!"

#That's a wrap!
