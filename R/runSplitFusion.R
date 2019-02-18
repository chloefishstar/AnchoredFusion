##!/bin/bash

## from BAM SA hits to structure variation (fusion, cDNA or large indel)
## input:
##	- bam/id_folder (current folder)
##	- bam/id_folder/consolidated.sam 
##
## output:
##	- strV_read_ID
##	- sa.cir0 (pre-filtered)
##	- sa.cir.panel (filtered)
##		## at least 25 bp exculsive between hit1 and hit2
##    		## at least 5 (or x$) split reads
##		## invovles targeting site (+/- 1kb)
##		## filtered within fusion < minBreakDistance


# enter your working directory
#. ../run.info.sh  ###Modified by Baifeng###
runSplitFusion <- function (step = c("bam-consolidate", "breakpoint-consolidate", "breakpoint-filter", "breakpoint-anno", "breakpoint-anno-post"), runInfo, output, sample.id){
	library(SplitFusion)
	options(width=204)
	source(runInfo)
#	subii='SUBID_'
#	dir.create(subii)
#	setwd(subii)
	dir.create(paste0(output,"/",sample.id), recursive = TRUE)
	setwd(paste0(output,"/",sample.id))
	SplitFusionPath <- system.file(package='SplitFusion')
	####==== 1: consolidated raw bam ====
	if ( "bam-consolidate" %in% step & !(file.exists(paste0(bam_path, "/", sample.id, ".consolidated.bam")) & file.size(paste0(bam_path, "/", sample.id, ".consolidated.bam")) !=0)) {
		system(paste0("bash ",SplitFusionPath,"/exec/SplitFusion.consolidated-bam.sh ",runInfo, " SA"))
	}

	##==== 2: consolidated.bam to breakpoint candidates ====
	if ( "breakpoint-consolidate" %in% step & !(file.exists("_breakpoint.noFilter3") & file.size("_breakpoint.noFilter3") !=0)) {

		system(paste0("bash ",SplitFusionPath,"/exec/SplitFusion.consolidated-breakpoint.sh ",runInfo))  ###Modified by Baifeng###

	}

	##==== 3: filters ====
	if ( "breakpoint-filter" %in% step & !(file.exists("breakpoint.candidates") & file.size("breakpoint.candidates") !=0)) {

		system(paste0("bash ",SplitFusionPath,"/exec/SplitFusion.breakpoint-filter.sh ",runInfo))
	}

	##==== 4: Annotate breakpoint gene, exon, cDNA position
	#if [ ! -f anno.left.right ]; then  ###Baifeng###
	#if [ ! -s mid.anno2 ]; then  ###Baifeng###
	#	bash $SplitFusionPath/scripts/SplitFusion.breakpoint-anno.sh $runInfo  ###Modified by Baifeng###
	#fi
	if ( "breakpoint-anno" %in% step & !(file.exists("mid.anno2") & file.size("mid.anno2") !=0)) {
        
        	system(paste0("bash ",SplitFusionPath,"/exec/SplitFusion.breakpoint-anno.sh ",runInfo))
	}



	##==== 5: further processing (in-frame status determination, etc)
	if ( "breakpoint-anno-post" %in% step) {

#		system(paste0("Rscript ",SplitFusionPath,"/R/SplitFusion.breakpoint-anno.postscript.R ",runInfo))  ###Modified by Baifeng###
		SplitFusion.breakpoint.anno.postscript(runInfo = runInfo)
	}

        system("rm _*")
}
