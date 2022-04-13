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

######################################################################################
######################################################################################
######################################################################################

#1) Check the quality of the reads

function check_quality {

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

	#input files: *.fastq.gz in ~/RNASeq_Data/Case and ../Control - hardcoded here, but this can be customized by the user.
        #output files: *_trimmed.fq files 
	#packages: trimgalore
		#trimgalore should cut sequences at illumina adapters

	###############################################################################

	module load trimgalore/0.6.6 #loading trimgalore module from ASC

	cd Case #move into the directory with the fastq files
        trim_galore --paired *R1.fastq.gz *R2.fastq.gz
        cd ../Control #move into directory with the control fastq files
	trim_galore --paired *R1.fastq.gz *R2.fastq.gz

}

function main {
#	check_quality /home/RNASeq_Data/
	trim_reads /home/RNASeq_Data/
}

main
