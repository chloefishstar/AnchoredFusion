# 

#system("sed 's/: /=/g' config.txt | sed 's:{::' | sed 's:}::' | sed 's/, /; /g' > config.R"
#source("config.R")
runAnchoredFusion <- function(runInfo){


#library(AnchoredFusion)
options(width=204)

#system("sed 's/: /=/g' config.txt | sed 's:{::' | sed 's:}::' | sed 's/, /\n/g' | sed \"s/'\(.*\)'=/\1=/g\" > config.r")
source(runInfo)

dir.create(paste0(output,"/",sample_id), recursive = TRUE, showWarnings = FALSE)
setwd(paste0(output,"/",sample_id))


step <- unlist(strsplit(step,","))


####==== 1: consolidated raw bam ====
	if ( "bam-consolidate" %in% step & !(file.exists(paste0(bam_path, "/", sample_id, ".consolidated.bam")) & file.size(paste0(bam_path, "/", sample_id, ".consolidated.bam")) !=0)) {
		suppressMessages(system(paste0("bash ",AnchoredFusionPath,"/exec/AnchoredFusion.consolidated-bam.sh ",runInfo, " SA")))

		cat("==========\n", date(), ": bam-consolidate completed successfully.\n==========\n")

	}

	##==== 2: consolidated.bam to breakpoint candidates ====
	if ( "breakpoint-consolidate" %in% step & !(file.exists("_breakpoint.noFilter3") & file.size("_breakpoint.noFilter3") !=0)) {

		suppressMessages(system(paste0("bash ",AnchoredFusionPath,"/exec/AnchoredFusion.consolidated-breakpoint.sh ",runInfo)))  ###Modified by Baifeng###

		cat("==========\n", date(), ": breakpoint-consolidate completed successfully.\n==========\n")

	}

	##==== 3: filters ====
	if ( "breakpoint-filter" %in% step & !(file.exists("breakpoint.candidates") & file.size("breakpoint.candidates") !=0)) {

		suppressMessages(system(paste0("bash ",AnchoredFusionPath,"/exec/AnchoredFusion.breakpoint-filter.sh ",runInfo)))

		cat("==========\n", date(), ": breakpoint-filter completed successfully.\n==========\n")

	}

	##==== 4: Annotate breakpoint gene, exon, cDNA position
	#if [ ! -f anno.left.right ]; then  ###Baifeng###
	#if [ ! -s mid.anno2 ]; then  ###Baifeng###
	#	bash $AnchoredFusionPath/scripts/AnchoredFusion.breakpoint-anno.sh $runInfo  ###Modified by Baifeng###
	#fi
	if ( "breakpoint-anno" %in% step & !(file.exists("mid.anno2") & file.size("mid.anno2") !=0)) {

        	suppressMessages(system(paste0("bash ",AnchoredFusionPath,"/exec/AnchoredFusion.breakpoint-anno.sh ",runInfo)))

		cat("==========\n", date(), ": breakpoint-anno completed successfully.\n==========\n")

	}



	##==== 5: further processing (in-frame status determination, etc)
	if ( "breakpoint-anno-post" %in% step) {

#		system(paste0("Rscript ",AnchoredFusionPath,"/R/AnchoredFusion.breakpoint-anno.postscript.R ",runInfo))  ###Modified by Baifeng###
		suppressMessages(AnchoredFusion.breakpoint.anno.postscript(runInfo = runInfo))

		cat("==========\n", date(), ": breakpoint-anno-post completed successfully.\n==========\n")

	}

        system("rm _*")


# Step 1...
# Step 2...

cat("==========\n", date(), ": Anchored-Fusion completed successfully.\n==========\n")
## END
}
