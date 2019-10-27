# 
runSplitFusion <- function(runInfo){

options(width=204)
source(runInfo)

dir.create(paste0(output,"/",sample_id), recursive = TRUE, showWarnings = FALSE)
setwd(paste0(output,"/",sample_id))

step <- unlist(strsplit(step,","))
cat("==========\n", date(), ": SplitFusion started successfully.\n==========\n")

##==== Step 1: fastq to consolidated bam (de-duplication based on unique UMI and adatptor ligation site) ====
		cat("==========\n", date(), ": Starting...1_fastq-bam...\n==========\n")
	#if ("bam-consolidate" %in% step & !(file.exists(paste0(bam_path, "/", sample_id, ".consolidated.bam")) & file.size(paste0(bam_path, "/", sample_id, ".consolidated.bam")) !=0)) {
	if ("1_fastq-bam" %in% step){
		suppressMessages(system(paste0("bash ", SplitFusionPath, "/exec/SplitFusion.1_fastq-bam.sh ",runInfo, " SA")))
		cat("==========\n", date(), ": 1_fastq-bam completed successfully.\n==========\n")
	}

##==== 2: consolidated.bam to breakpoint candidates (no filter) ====
		cat("==========\n", date(), ": Starting...2_bam-breakpoint...\n==========\n")
	#if ( "BamToBreakpoint" %in% step & !(file.exists("breakpoint.candidates.preFilter") & file.size("breakpoint.candidates.preFilter") !=0)) {
	if ("2_bam-breakpoint" %in% step) {
		suppressMessages(system(paste0("bash ", SplitFusionPath, "/exec/SplitFusion.2_bam-breakpoint.sh ",runInfo)))
		cat("==========\n", date(), ": 2_bam-breakpoint completed successfully.\n==========\n")
	}

##==== 3: filters ====
		cat("==========\n", date(), ": Starting...3_breakpoint-filter...\n==========\n")
	if ("3_breakpoint-filter" %in% step) {
		suppressMessages(system(paste0("bash ", SplitFusionPath, "/exec/SplitFusion.3_breakpoint-filter.sh ",runInfo)))
		cat("==========\n", date(), ": 3_breakpoint-filter completed successfully.\n==========\n")
	}

##==== 4: Annotate breakpoint gene, exon, cDNA position
		cat("==========\n", date(), ": Starting...4_breakpoint-anno...\n==========\n")
	if ("4_breakpoint-anno" %in% step) {
        	suppressMessages(system(paste0("bash ", SplitFusionPath, "/exec/SplitFusion.4_breakpoint-anno.sh ",runInfo)))
		cat("==========\n", date(), ": 4_breakpoint-anno completed successfully.\n==========\n")
	}

##==== 5: further processing (in-frame status determination, etc)
		cat("==========\n", date(), ": Starting...5_breakpoint-anno-post...\n==========\n")
	if ( "5_breakpoint-anno-post" %in% step) {
		suppressMessages(SplitFusion.breakpoint.anno.postscript(runInfo = runInfo))
		cat("==========\n", date(), ": 5_breakpoint-anno-post completed successfully.\n==========\n")
	}

        system("rm _*")

cat("==========\n", date(), ": Split-Fusion completed successfully.\n==========\n")
## END
}
