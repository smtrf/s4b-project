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

	module load fastqc/0.11.9 #loading the fastqc module

	mkdir FastQC #make a new directory for output files
	cd Case #move into the directory with the fastq files
	fastqc *.fastqc -o ~/RNASeq_Data/FastQC #Use FastQC to perform a quality check of sequences
	cd ../Control #move into directory with the control fastq files
	fastqc *.fastqc -o ~/RNASeq_Data/FastQC #Use FastQC to perform a quality check of sequences

}

#2) Clean the reads by trimming


function main {
	check_quality /home/RNASeq_Data/
}

main
