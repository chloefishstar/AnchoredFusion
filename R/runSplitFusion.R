# 
runSplitFusion <- function(runInfo){

options(width=204)
source(runInfo)

dir.create(paste0(output,"/",sample_id), recursive = TRUE, showWarnings = FALSE)
setwd(paste0(output,"/",sample_id))

step <- unlist(strsplit(step,","))
cat("==========\n", date(), ": SplitFusion started successfully.\n==========\n")

####==== 1: consolidated raw bam ====
		cat("==========\n", date(), ": Starting...bam-consolidate...\n==========\n")
	if ( "bam-consolidate" %in% step & !(file.exists(paste0(bam_path, "/", sample_id, ".consolidated.bam")) & file.size(paste0(bam_path, "/", sample_id, ".consolidated.bam")) !=0)) {
		suppressMessages(system(paste0("bash ",SplitFusionPath,"/exec/SplitFusion.consolidated-bam.sh ",runInfo, " SA")))
		cat("==========\n", date(), ": bam-consolidate completed successfully.\n==========\n")
	}

##==== 2: consolidated.bam to breakpoint candidates (no filter) ====
		cat("==========\n", date(), ": Starting...BamToBreakpoint...\n==========\n")
	if ( "BamToBreakpoint" %in% step & !(file.exists("breakpoint.candidates.preFilter") & file.size("breakpoint.candidates.preFilter") !=0)) {
		suppressMessages(system(paste0("bash ",SplitFusionPath,"/exec/SplitFusion.BamToBreakpoint.sh ",runInfo)))  ###Modified by Baifeng###
		cat("==========\n", date(), ": BamToBreakpoint completed successfully.\n==========\n")
	}

##==== 3: filters ====
		cat("==========\n", date(), ": Starting...breakpoint-filter...\n==========\n")
	if ( "breakpoint-filter" %in% step & !(file.exists("breakpoint.candidates") & file.size("breakpoint.candidates") !=0)) {
		suppressMessages(system(paste0("bash ",SplitFusionPath,"/exec/SplitFusion.breakpoint-filter.sh ",runInfo)))
		cat("==========\n", date(), ": breakpoint-filter completed successfully.\n==========\n")
	}

##==== 4: Annotate breakpoint gene, exon, cDNA position
		cat("==========\n", date(), ": Starting...breakpoint-anno...\n==========\n")
	if ( "breakpoint-anno" %in% step & !(file.exists("anno.left.right") & file.size("anno.left.right") !=0)) {
        	suppressMessages(system(paste0("bash ",SplitFusionPath,"/exec/SplitFusion.breakpoint-anno.sh ",runInfo)))
		cat("==========\n", date(), ": breakpoint-anno completed successfully.\n==========\n")
	}

##==== 5: further processing (in-frame status determination, etc)
		cat("==========\n", date(), ": Starting...breakpoint-anno-post...\n==========\n")
	if ( "breakpoint-anno-post" %in% step) {
#		system(paste0("Rscript ",SplitFusionPath,"/R/SplitFusion.breakpoint-anno.postscript.R ",runInfo))  ###Modified by Baifeng###
		suppressMessages(SplitFusion.breakpoint.anno.postscript(runInfo = runInfo))
		cat("==========\n", date(), ": breakpoint-anno-post completed successfully.\n==========\n")
	}

        system("rm _*")


# Step 1...
# Step 2...

cat("==========\n", date(), ": Split-Fusion completed successfully.\n==========\n")
## END
}
