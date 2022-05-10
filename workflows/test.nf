#!/bin/bash nextflow

nextflow.enable.dsl=2

    // Step 0 NF way
    chr_ch = Channel.from( 1..22 )

    arrayId_ch = Channel.from( 1..2 ) // 1-573 phenotypes

    slidingWindowPos_ch = Channel.from( 1..2 ) // for loop slidingWindow 1-200

    agdsFiles_ch = Channel.fromPath(params.agdsFiles, checkIfExists:true)
    aGDSdir_ch = Channel.fromPath(params.aGDSdir, checkIfExists:true).view()
    jobNum_ch = Channel.fromPath(params.jobNum, checkIfExists:true)
    nullModel_ch = Channel.fromPath(params.nullModel, checkIfExists:true)
    nameCatalog_ch = Channel.fromPath(params.nameCatalog, checkIfExists:true)
/*
========================================================================================
    IMPORT LOCAL MODULES/SUBWORKFLOWS
========================================================================================
*/

    // Step 0: Preparation for association analysis of whole-genome/whole-exome sequencing studies
        // Input: aGDS files of all 22 chromosomes. For more details, please see the R script.
        // Output: agds_dir.Rdata, Annotation_name_catalog.Rdata, jobs_num.Rdata.
        // Script: Association_Analysis_PreStep.r
        // path new script : /nfs/team151/software/STAARpipeline_INTERVAL/final

    process analysisPreStep {     
        input:
            path aGDS

        output:
            path "*_dir.Rdata", emit: agds_dir

        script:
        """
        #!/usr/bin/env Rscript

        ## load required packages
        library(gdsfmt)
        library(SeqArray)
        library(SeqVarTools)

        ###############################
        #           Input
        ###############################

## file directory of aGDS file (genotype and annotation data) 
        # dir.geno <- "/lustre/scratch119/realdata/mdt2/projects/interval_wgs/final_release_freeze_GDS/gt_phased_GDS/"
dir_geno <- "${params.aGDSdir}"

## file name of aGDS, separate by chr number 
adgs_file_name_1 <- "interval_wgs.chr"
agds_file_name_2 <- ".gt_phased.gds"

## channel name of the QC label in the GDS/aGDS file
        #QC_label <- "annotation/info/QC_label"
QC_label <- "${params.qcLabel}"

## file directory for the output files
        #output_path <- "/lustre/scratch119/realdata/mdt2/projects/interval_wgs/analysis/STAARpipeline/data/input/"
output_path <- "${params.output}"

        ###############################
        #        Main Function
        ###############################

#### aGDS directory
agds_dir <- paste0(dir_geno,adgs_file_name_1,seq(1,22),agds_file_name_2) 
        #save(agds_dir,file=paste0(output_path,"agds_dir.Rdata",sep=""))
save(agds_dir,file=paste0(".","agds_dir.Rdata",sep=""))

#### Annotation dir -> SEEMS ITS NOT NEEDED IN THIS STEP
        #Annotation_name_catalog <- "/lustre/scratch119/realdata/mdt2/projects/interval_wgs/analysis/STAARpipeline/data/input/Annotation_name_catalog.txt"
#Annotation_name_catalog <- ${params.annotationNameCatalog}

#### jobs_num
jobs_num <- matrix(rep(0,66),nrow=22)
for(chr in 1:22)
{
    print(chr)
    gds.path <- agds_dir[chr]
    genofile <- seqOpen(gds.path)
    
    filter <- seqGetData(genofile, QC_label)
    SNVlist <- filter == "PASS" 
    
    position <- as.numeric(seqGetData(genofile, "position"))
    position_SNV <- position[SNVlist]
    
    jobs_num[chr,1] <- chr
    jobs_num[chr,2] <- min(position[SNVlist])
    jobs_num[chr,3] <- max(position[SNVlist])
    
    seqClose(genofile)
}

# Individual Analysis
jobs_num <- cbind(jobs_num,ceiling((jobs_num[,3]-jobs_num[,2])/10e6))

# Sliding Window
jobs_num <- cbind(jobs_num,ceiling((jobs_num[,3]-jobs_num[,2])/5e6))

# SCANG
jobs_num <- cbind(jobs_num,ceiling((jobs_num[,3]-jobs_num[,2])/1.5e6))

colnames(jobs_num) <- c("chr","start_loc","end_loc","individual_analysis_num","sliding_window_num","scang_num")
jobs_num <- as.data.frame(jobs_num)

        # save(jobs_num,file=paste0(output_path,"jobs_num.Rdata",sep=""))
save(jobs_num,file=paste0(".","jobs_num.Rdata",sep=""))

        """
    }
workflow TEST {
    //take:
    //    arrayId
    //    slidingWindowPos_ch

    //main:
    //Step 0: Preparation for association analysis of whole-genome/whole-exome sequencing studies
    analysisPreStep(aGDSdir_ch)

}