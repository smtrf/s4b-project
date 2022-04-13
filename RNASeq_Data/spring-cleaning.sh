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

	module load fastqc/0.11.9 #loading the fastqc module from ASC

	mkdir FastQC #make a new directory for output files
	cd Case #move into the directory with the fastq files
	fastqc *.fastq.gz -o ~/s4b-project/RNASeq_Data/FastQC #Use FastQC to perform a quality check of sequences
	cd ../Control #move into directory with the control fastq files
	fastqc *.fastq.gz -o ~/s4b-project/RNASeq_Data/FastQC #Use FastQC to perform a quality check of sequences

}

#2) Clean the reads by trimming


function main {
	check_quality /home/RNASeq_Data/
}

main
