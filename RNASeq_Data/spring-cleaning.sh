#!/bin/bash

#Stephen T. and Ally S. 
#S4B
#Spring 2022

######################################################################################
##############################RNA-seq Data Cleanup####################################
######################################################################################

#In this script we will:
	#Check the quality of the reads
	#Clean reads
	#Map the reads to a genome
	#Get read counts per gene
	#Use the read counts to find differentially expressed genes (DEGs)

#Packages: FastQC, TrimGalore, 

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
        #output files: *_trimmed.fq files 
		#These files will be located in ~/s4b-project/RNASeq_Data/TrimmedReads
	#packages: trimgalore
		#trimgalore will automatically detect and cut sequences at illumina adapters

	###############################################################################

	module load trimgalore/0.6.6 #loading trimgalore module from ASC

	#making new directories in ~/s4b-project/RNASeq_Data for the trimmed data to be stored
	mkdir TrimmedReads
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
                #This file will be located in ~/s4b-project/RNASeq_Data/TrimmedReads/trimmed-FastQC as specified in the code. This can be modified by the user if th$
        #Packages: FastQC

        ###############################################################################

        module load fastqc/0.11.9 #loading the fastqc module from ASC

        mkdir trimmed-FastQC #make a new directory for output files
        
	cd Case #move into the directory with the fastq files
        fastqc *.fq.gz -o ~/s4b-project/RNASeq_Data/TrimmedReads/trimmed-FastQC #Use FastQC to perform a quality check of sequences
        cd ../Control #move into directory with the control fastq files
        fastqc *.fq.gz -o ~/s4b-project/RNASeq_Data/TrimmedReads/trimmed-FastQC #Use FastQC to perform a quality check of sequences

}

function main {
#	quality_check /home/RNASeq_Data/
#	trim_reads /home/RNASeq_Data/
	qc_trimmed /home/RNASeq_Data/TrimmedReads
}

main
