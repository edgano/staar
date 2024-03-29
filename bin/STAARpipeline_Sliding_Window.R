#####################################################################
# Sliding window analysis using STAARpipeline
# Xihao Li, Zilin Li
# 11/04/2021
#####################################################################

rm(list=ls())
gc()

### load required package
library(gdsfmt, lib="/nfs/team151/software/Rlibs_4.1")
library(SeqArray, lib="/nfs/team151/software/Rlibs_4.1")
library(SeqVarTools, lib="/nfs/team151/software/Rlibs_4.1")
library(STAAR, lib="/nfs/team151/software/Rlibs_4.1")
library(STAARpipeline, lib="/nfs/team151/software/Rlibs_4.1")

###########################################################
#           User Input
###########################################################
## job nums
jobs_num <- get(load("/lustre/scratch119/realdata/mdt2/projects/interval_wgs/analysis/STAARpipeline/data/input/jobs_num.Rdata"))
## agds dir
agds_dir <- get(load("/lustre/scratch119/realdata/mdt2/projects/interval_wgs/analysis/STAARpipeline/data/input/agds_dir.Rdata"))
## Null Model
obj_nullmodel <- get(load("/lustre/scratch119/realdata/mdt2/projects/interval_wgs/analysis/STAARpipeline/results/Null_Model/obj.SCANG_STAAR.fbc_neut.Rdata"))
trait <- "fbc_neut"

## QC_label
QC_label <- "annotation/info/QC_label"
## variant_type
variant_type <- "SNV"
## geno_missing_imputation
geno_missing_imputation <- "mean"

## Annotation_dir
Annotation_dir <- "annotation/info/FunctionalAnnotation/FunctionalAnnotation"
## Annotation channel
Annotation_name_catalog <- read.delim("/lustre/scratch119/realdata/mdt2/projects/interval_wgs/analysis/STAARpipeline/data/input/Annotation_name_catalog.txt")
## Use_annotation_weights
Use_annotation_weights <- TRUE
## Annotation name
Annotation_name <- c("CADD","LINSIGHT","FATHMM.XF","aPC.EpigeneticActive","aPC.EpigeneticRepressed","aPC.EpigeneticTranscription",
					"aPC.Conservation","aPC.LocalDiversity","aPC.Mappability","aPC.TF","aPC.Protein")
					
## output path
output_path <- paste("/lustre/scratch119/realdata/mdt2/projects/interval_wgs/analysis/STAARpipeline/results/Sliding_Window/", trait, sep="")
cmd <- paste("mkdir", output_path)
system(cmd)
## output file name
output_file_name <- paste("results_sliding_window_", trait, sep="")

## input array id from batch file 
arrayid <- as.numeric(commandArgs(TRUE)[1])     ## from 1 to max(cumsum(jobs_num$sliding_window_num)) which is 573


###########################################################
#           Main Function 
###########################################################
chr <- which.max(arrayid <= cumsum(jobs_num$sliding_window_num))
group.num <- jobs_num$sliding_window_num[chr]

if (chr==1){
   groupid <- arrayid
}  else  {
   groupid <- arrayid - cumsum(jobs_num$sliding_window_num)[chr-1]
}

### gds file
gds.path <- agds_dir[chr]
genofile <- seqOpen(gds.path)

start_loc <- (groupid-1)*5e6 + jobs_num$start_loc[chr]
end_loc <- start_loc + 1000*25 - 1

results_sliding_window <- c()
for(kk in 1:200)
{
	print(kk)

	start_loc_sub <- start_loc + 1000*25*(kk-1)
	end_loc_sub <- end_loc + 1000*25*(kk-1) + 1000
	
	end_loc_sub <- min(end_loc_sub,jobs_num$end_loc[chr])
	
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

save(results_sliding_window, file=paste0(output_path,"/",output_file_name,"_",arrayid,".Rdata"))

seqClose(genofile)



## r <- get(load("/lustre/scratch119/realdata/mdt2/projects/interval_wgs/analysis/STAARpipeline/results/Sliding_Window/fbc_neut/results_sliding_window_fbc_neut_1.Rdata"))