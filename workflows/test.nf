#!/bin/bash nextflow

nextflow.enable.dsl=2

    // Step 0 NF way
    chr_ch = Channel.from( 1..22 )

    arrayId_ch = Channel.from( 1..2 ) // 1-573 phenotypes

    slidingWindowPos_ch = Channel.from( 1..2 ) // for loop slidingWindow 1-200

    // agdsFiles_ch = Channel.fromPath(params.agdsFiles, checkIfExists:true)
    aGDSdir_ch = Channel.fromPath(params.aGDSdir, checkIfExists:true)

    // Step 1 NF way
    phenoDir = "/lustre/scratch119/humgen/projects/interval_wgs/analysis/phenotypes/all_phenotypes/fam/by_trait/wgs_cov_adj"
    pNum= '2'

    // Step 3.2
    //jobNum = params.jobNum
    //variantType = params.variantType
    //gene_missing_imputation = params.gene_missing_imputation

    /*
========================================================================================
    IMPORT LOCAL MODULES/SUBWORKFLOWS
========================================================================================
*/

    // Step 0: Preparation for association analysis of whole-genome/whole-exome sequencing studies
        // Input: aGDS files of all 22 chromosomes. For more details, please see the R script.
        // Output:  agds_dir.Rdata, 
        //          ~~Annotation_name_catalog.Rdata, 
        //          jobs_num.Rdata.
        // Script: Association_Analysis_PreStep.r
        // path new script : /nfs/team151/software/STAARpipeline_INTERVAL/final

    process analysisPreStep {    
        publishDir "${params.outdir}/step00", mode: 'copy', overwrite: false, pattern: "*_dir.Rdata"
        publishDir "${params.outdir}/step00", mode: 'copy', overwrite: false, pattern: "*_num.Rdata"
        
        input:
            path aGDS
        output:
            path "*_dir.Rdata", emit: agds_dirFile
            path '*_num.Rdata', emit: jobs_num
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
        output_path <- "./"

                ###############################
                #        Main Function
                ###############################

        #### aGDS directory
        agds_dir <- paste0(dir_geno,adgs_file_name_1,seq(1,22),agds_file_name_2) 
                #save(agds_dir,file=paste0(output_path,"agds_dir.Rdata",sep=""))
        save(agds_dir,file=paste0(output_path,"agds_dir.Rdata",sep=""))

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
        save(jobs_num,file=paste0(output_path,"jobs_num.Rdata",sep=""))

        """
    }
    process fitNullModel { //   SCRIPT ->>        /nfs/team151/software/STAARpipeline_INTERVAL/final/STAARpipeline_Null_Model.R
        publishDir "${params.outdir}/step01", mode: 'copy', overwrite: false, pattern: "obj.STAAR.*.Rdata"

        input:
            path phenoDir
            val pNum
        output:
            path "obj.STAAR.*.Rdata", emit: objNullModel
        script:
        """
        #!/usr/bin/env Rscript

        library(gdsfmt)
        library(SeqArray)
        library(SeqVarTools)
        library(STAAR)
        library(STAARpipeline)
        library(data.table)         ## added for fread function
        ###############################
        #           Input
        ###############################
        #       pnum <- as.numeric(commandArgs(TRUE)[1])
        pnum <- as.numeric(${pNum})

        ## Phenotype data
        ##      d <- dir("/lustre/scratch119/humgen/projects/interval_wgs/analysis/phenotypes/all_phenotypes/fam/by_trait/wgs_cov_adj", full.names=TRUE)
        d <- dir("${phenoDir}", full.names=TRUE)
        ##      pname <- dir("/lustre/scratch119/humgen/projects/interval_wgs/analysis/phenotypes/all_phenotypes/fam/by_trait/wgs_cov_adj")
        pname <- dir("${phenoDir}")
        pname <- gsub(".fam", "", pname)
        pheno <- fread(d[pnum], stringsAsFactors=F, data.table=F)
        print(pheno)
        ## file directory for the output file 
        ##      output_path <- "/lustre/scratch119/realdata/mdt2/projects/interval_wgs/analysis/STAARpipeline/results/Null_Model/"
        output_path <- "./"
        ## output file name
        output_name <- paste("obj.STAAR.",pname[pnum],".Rdata", sep="")


        ###########################################################
        #           Main Function 
        ###########################################################

        ### fit null model
        #obj.nullmodel <- fit_nullmodel(LDLadj.norm~age+age2+sex+PC1+PC2+PC3+PC4+PC5+PC6+PC7+PC8+PC9+PC10, data = phenotype, kins = sgrm, kins_cutoff = 0.022, id = "userId", use_sparse = TRUE,family = gaussian(link = "identity"), verbose=T)
        obj.nullmodel <- fit_nullmodel(V6~1, data=pheno, kins=NULL, id="V1", family=gaussian(link="identity"), verbose=T)

        save(obj.nullmodel,file = paste0(output_path,output_name))


        ## Edit file name - not needed anymore
        ## d <- dir("/lustre/scratch119/realdata/mdt2/projects/interval_wgs/analysis/STAARpipeline/results/Null_Model", full.names=TRUE)
        ## for (i in 1:length(d))  {
        ##     print(i)
        ##     load(d[i])   ## obj.nullmodel
        ##     outfile <- gsub(".fam", "", d[i])
        ##     save(obj.nullmodel, file=outfile)
        ## }

        """
    }
    // Step 2: Individual analysis
        // Input: aGDS files and the STAAR null model. For more details, please see the R script.
        // Output: Rdata files with the user-defined names.
        // Script: STAARpipeline_Individual_Analysis.r

        //TODO -> needed the jobs_num.Rdata TOO
        //        arrayid <- as.numeric(commandArgs(TRUE)[1])        

    // Step 3.1: Gene-centric coding analysis   
        // Input: aGDS files and the STAAR null model. For more details, please see the R scripts.
        // Output: 381 Rdata files with the user-defined names. For more details, please see the R scripts.
        // Script: STAARpipeline_Gene_Centric_Coding.r and STAARpipeline_Gene_Centric_Coding_Long_Masks.r

        // TODO -> gene_num_in_array <- 50  
        //         table(genes_info[,2])
        //  ### exclude large genes     <-- magic numbers

    // Step 3.2: Gene-centric noncoding analysis
        // Input: aGDS files and the STAAR null model. For more details, please see the R scripts.
        // Output: 387 Rdata files with the user-defined names for protein-coding genes and 223 Rdata files with the user-defined names for ncRNA genes. For more details, please see the R scripts.
        //          Script: STAARpipeline_Gene_Centric_Noncoding.r, STAARpipeline_Gene_Centric_Noncoding_Long_Masks.r, STAARpipeline_Gene_Centric_ncRNA.r and STAARpipeline_Gene_Centric_ncRNA_Long_Masks.r
        // SCRIPT --> /nfs/team151/software/STAARpipeline_INTERVAL/final/STAARpipeline_Gene_Centric_Noncoding.R

    process geneCentricNoCoding {   
            publishDir "${params.outdir}/step3_2", mode: 'copy', overwrite: false, pattern: "results_gene_centric_noncoding_*"
        
        input:
            path jobsNum
            file agds_dirFile
            file nullModel
            val variantType
            val gene_missing_imputation

        output:
            path "results_gene_centric_noncoding_*", emit: results

        script:
        """
        #!/usr/bin/env Rscript

        ## load required package
        library(gdsfmt)
        library(SeqArray)
        library(SeqVarTools)
        library(STAAR)
        library(STAARpipeline)
        library(GenomeInfoDb)
        library(TxDb.Hsapiens.UCSC.hg38.knownGene)

        ###############################
        #           Input
        ###############################
            ## job nums
        #   jobs_num <- get(load("/lustre/scratch119/realdata/mdt2/projects/interval_wgs/analysis/STAARpipeline/data/input/jobs_num.Rdata"))
        jobs_num <- get(load("${jobsNum}"))
            ## agds dir
        #   agds_dir <- get(load("/lustre/scratch119/realdata/mdt2/projects/interval_wgs/analysis/STAARpipeline/data/input/agds_dir.Rdata"))
        agds_dir <- get(load("${agds_dirFile}"))
            ## Null Model
        #   obj_nullmodel <- get(load("/lustre/scratch119/realdata/mdt2/projects/interval_wgs/analysis/STAARpipeline/results/Null_Model/obj.STAAR.fbc_neut.Rdata"))
        obj_nullmodel <- get(load("${nullModel}"))

        #trait <- "fbc_neut"
        trait <- "${params.trait}"

            ## QC_label
        #   QC_label <- "annotation/info/QC_label"
        QC_label <- "${params.qcLabel}"

            ## variant_type
        #   variant_type <- "SNV"
        variant_type <- "${variantType}"

            ## geno_missing_imputation
        #   geno_missing_imputation <- "mean"
        geno_missing_imputation <- "${gene_missing_imputation}"

            ## Annotation_dir
        #   Annotation_dir <- "annotation/info/FunctionalAnnotation/FunctionalAnnotation"
        Annotation_dir <- "${params.annotationDir}"
            ## Annotation channel
        #   Annotation_name_catalog <- read.delim("/lustre/scratch119/realdata/mdt2/projects/interval_wgs/analysis/STAARpipeline/data/input/Annotation_name_catalog.txt", stringsAsFactors=F)
        Annotation_name_catalog <- read.delim("${params.annotationNameCatalog}", stringsAsFactors=F)
            ## Use_annotation_weights
        #   Use_annotation_weights <- TRUE
        Use_annotation_weights <- ${params.AnnotationWeights}
            
            ## Annotation name
        Annotation_name <- c("CADD","LINSIGHT","FATHMM.XF","aPC.EpigeneticActive","aPC.EpigeneticRepressed","aPC.EpigeneticTranscription",
                            "aPC.Conservation","aPC.LocalDiversity","aPC.Mappability","aPC.TF","aPC.Protein","aPC.Liver")

            ## output path
        #   output_path <- paste("/lustre/scratch119/realdata/mdt2/projects/interval_wgs/analysis/STAARpipeline/results/Gene_Centric_NonCoding/", trait, sep="")
        output_path <- paste("./")

        #cmd <- paste("mkdir", output_path)
        #system(cmd)
            ## output file name
        output_file_name <- paste("results_gene_centric_noncoding_", trait, sep="")

            ## input array id from batch file (Harvard FAS cluster)
        #   arrayid <- as.numeric(commandArgs(TRUE)[1])                    ## from 1 to sum(group.num.allchr) which is 379
        arrayid <- as.numeric("${params.ArrayId}")

        ###########################################################
        #           Main Function 
        ###########################################################
        ## gene number in job
# FIX THE LENGHT
        gene_num_in_array <- 2 #50 
        group.num.allchr <- ceiling(table(genes_info[,2])/gene_num_in_array)
        sum(group.num.allchr)                                        

        chr <- which.max(arrayid <= cumsum(group.num.allchr))
        group.num <- group.num.allchr[chr]

        if (chr == 1){
        groupid <- arrayid
        }  else  {
        groupid <- arrayid - cumsum(group.num.allchr)[chr-1]
        }

        genes_info_chr <- genes_info[genes_info[,2]==chr,]
        sub_seq_num <- dim(genes_info_chr)[1]

        if (groupid < group.num)
        { 
            sub_seq_id <- ((groupid - 1)*gene_num_in_array + 1):(groupid*gene_num_in_array)
        }  else  {
            sub_seq_id <- ((groupid - 1)*gene_num_in_array + 1):sub_seq_num
        }	
# FIX THE LENGHT
        ### exclude large genes
        jobid_exclude <- c(21,39) # ,44,45,46,53,55,83,88,103,114,127,135,150,154,155,163,164,166,180,189,195,200,233,280,285,295,313,318,319,324,327,363,44,45,54)
        sub_seq_id_exclude <- c(1009,1929,182,214,270,626,741,894,83,51,611,385,771,493,671,702,238,297,388,352,13,303,600,170,554,207,724,755,1048,319,324,44,411,195,236,677)
            
        for (i in 1:length(jobid_exclude))
        {
            if (arrayid==jobid_exclude[i])
            {
                sub_seq_id <- setdiff(sub_seq_id,sub_seq_id_exclude[i])
            }
        }

        ### gds file
        gds.path <- agds_dir[chr]
        genofile <- seqOpen(gds.path)

        genes <- genes_info

        results_noncoding <- c()
        for(kk in sub_seq_id)
        {
            print(kk)
            gene_name <- genes_info_chr[kk,1]
            results <- Gene_Centric_Noncoding(chr=chr,gene_name=gene_name,genofile=genofile,obj_nullmodel=obj_nullmodel,
                                        rare_maf_cutoff=0.01,rv_num_cutoff=2,
                                        QC_label=QC_label,variant_type=variant_type,geno_missing_imputation=geno_missing_imputation,
                                        Annotation_dir=Annotation_dir,Annotation_name_catalog=Annotation_name_catalog,
                                        Use_annotation_weights=Use_annotation_weights,Annotation_name=Annotation_name)
            
            results_noncoding <- append(results_noncoding,results)
        }

        save(results_noncoding, file=paste0(output_path,"/",output_file_name,"_",arrayid,".Rdata"))

        seqClose(genofile)
        """
    }

// Step 4: Sliding window analysis
        // Input: aGDS files and the STAAR null model. For more details, please see the R script.
        // Output: Rdata files with the user-defined names.
        // Script: STAARpipeline_Sliding_Window.r
        //          /nfs/team151/software/STAARpipeline_INTERVAL/final/STAARpipeline_Sliding_Window.R

    process slidingWindow {  
        publishDir "${params.outdir}/step04", mode: 'copy', overwrite: false, pattern: "results_sliding_window_*.Rdata"
        tag "arrayId - $arrayId"

        input:
            val (arrayId) 
            file (aGDSdir) 
            path (nullModel) 
            path (jobNum) 

        output:
            path "*.Rdata", emit: slidingWindow_out

        script:
        """
        #!/usr/bin/env Rscript
    # modified the library path from lustre to docker
        library(gdsfmt)        
        library(SeqArray)
        library(SeqVarTools)
        library(STAAR)
        library(STAARpipeline)
        ###############################
        #           Input
        ###############################
        ##  ## LOAD R OBJECTS
            ## job nums
        jobs_num <- get(load("${jobNum}"))
            ## agds dir
        agds_dir <- get(load("${aGDSdir}"))
            ## Null Model
        obj_nullmodel <- get(load("${nullModel}"))

    ## defined in the bash 1-573
        ## from 1 to max(cumsum(jobs_num\$sliding_window_num)) which is 573
        arrayid <- as.numeric("${arrayId}")

        #### LABELS
        # trait <- "fbc_neut"  # used in #output_path <- paste( .... and not used
            ## QC_label                 --> used in the TRY
        QC_label <- "annotation/info/QC_label"
            ## variant_type             --> used in the TRY
        variant_type <- "SNV"
            ## geno_missing_imputation  --> used in the TRY
        geno_missing_imputation <- "mean"

        ##  ## ANNOTATION
    # WHY? are thet for input or for output?
            ## Annotation_dir
        Annotation_dir <- "annotation/info/FunctionalAnnotation/FunctionalAnnotation"
            ## Annotation channel
        Annotation_name_catalog <- read.delim("${params.annotationNameCatalog}")

    # boolean by default, maybe can be a param of pipeline ?
            ## Use_annotation_weights
        Use_annotation_weights <- TRUE
    # same, why? input or output? 
    #   is static?
            ## Annotation name      ## size = 11
        Annotation_name <- c("CADD","LINSIGHT","FATHMM.XF","aPC.EpigeneticActive","aPC.EpigeneticRepressed","aPC.EpigeneticTranscription",
                            "aPC.Conservation","aPC.LocalDiversity","aPC.Mappability","aPC.TF","aPC.Protein")
        
        ##  ## OUTPUT
            ## output path
        #   output_path <- paste("/lustre/scratch119/realdata/mdt2/projects/interval_wgs/analysis/STAARpipeline/results/Sliding_Window/", trait, sep="")
        output_path <- paste("./")
        #   cmd <- paste("mkdir", output_path)
        #   system(cmd)
            ## output file name
    # can be defined at the begining of the script
        output_file_name <- paste("results_sliding_window_", "", sep="")

    ## input array id from batch file               
    #    SBATCH --array=1-573 --mem=11000
        #arrayid <- as.numeric(commandArgs(TRUE)[1])     ## from 1 to max(cumsum(jobs_num\$sliding_window_num)) which is 573

        ###############################
        #        Main Function
        ###############################
        chr <- which.max(arrayid <= cumsum(jobs_num\$sliding_window_num))
        group.num <- jobs_num\$sliding_window_num[chr]

        if (chr==1){
            groupid <- arrayid
        }  else  {
            groupid <- arrayid - cumsum(jobs_num\$sliding_window_num)[chr-1]
        }

        ### gds file
        gds.path <- agds_dir[chr]
        genofile <- seqOpen(gds.path)

        start_loc <- (groupid-1)*5e6 + jobs_num\$start_loc[chr]
        end_loc <- start_loc + 1000*25 - 1

        results_sliding_window <- c()

    #>>  TODO  << This should be unrapped
    # Why it is 200 and not 2k ? -> this can be unwrapped with a ch.value(1..200)
        for(kk in 1:200)              
        {
            print(kk)

            start_loc_sub <- start_loc + 1000*25*(kk-1)
            end_loc_sub <- end_loc + 1000*25*(kk-1) + 1000
            
            end_loc_sub <- min(end_loc_sub,jobs_num\$end_loc[chr])

            # If unwrapped, all the files of results() will need to be merged after the process            
            results <- c()
            if(start_loc_sub < end_loc_sub)
            {
                results <- try(Sliding_Window(chr=chr, start_loc=start_loc_sub, end_loc=end_loc_sub, genofile=genofile, obj_nullmodel=obj_nullmodel, 
                                type="multiple",QC_label=QC_label,variant_type=variant_type,geno_missing_imputation=geno_missing_imputation,
                                Annotation_dir=Annotation_dir,Annotation_name_catalog=Annotation_name_catalog,
                                Use_annotation_weights=Use_annotation_weights,Annotation_name=Annotation_name))
                
                if(class(results)!="try-error")
                {
                    results_sliding_window <- rbind(results_sliding_window,results)
                }

            }
        }
    # saving path '.' and NF will handle the file
        save(results_sliding_window, file=paste0(".","/",output_file_name,"_",arrayid,".Rdata"))

        seqClose(genofile)
        """
    }

workflow TEST {
    //take:
    //    arrayId
    //    slidingWindowPos_ch

    //main:
    //Step 0: Preparation for association analysis of whole-genome/whole-exome sequencing studies
    analysisPreStep(aGDSdir_ch)

    //step 01
    fitNullModel(phenoDir, pNum)

    //step 3.2
    geneCentricNoCoding(params.jobNum,
                        analysisPreStep.out.agds_dirFile,
                        fitNullModel.out.objNullModel,
                        params.variantType,
                        params.gene_missing_imputation)
    // step 04
    slidingWindow(arrayId_ch, 
                analysisPreStep.out.agds_dirFile, 
                fitNullModel.out.objNullModel, 
                params.jobNum)
}