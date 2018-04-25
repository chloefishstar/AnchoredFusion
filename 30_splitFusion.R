
###################
## Structure variation detection, including fusion, large indel, structure variation...
##	-- based on BWA MEM SA reads
###########################
## To run as a plugin independent of the whole pipeline:
## 1). bsub -Is -q interact /bin/bash
## 2). cd to analysis run folder
## 3). if background SA needs to be filtered, copy the filter file to 'bam' folder or
##	, if not applied to all or different filter files for different samples
##	, bam/A##-P7## folder and renamed as 'sa.filter', e.g.:
##	cp ~/projects/criStr/rd2.filter bam/sa.filter  ## or
##	cp ~/projects/criStr/A01.filter bam/A01-P701/sa.filter
## then:
##      . ./run.info.sh; Rscript  $SplitFusionPath/30_splitFusion.R &

source('run.info')
options(width=204)

## the iiFreq.txt is the sequencing samplesheet file that contains list of samples with the field AP7 (adaptor+P7) barcode.
##
ii = read.table("iiFreq.txt",header=T, sep='\t', stringsAsFactors=F, fill=T)
head(ii)
dir.create('fusion', showWarnings = FALSE, recursive = TRUE)
setwd('fusion')

## prep jobs
# (i=ii$AP7[1])
for (i in ii$AP7){
	system(paste("cat ", SplitFusionPath, "/30_splitFusion.sh | sed 's/SUBID_/", i,"/g' > _job_.", i, sep=''))
        source(paste(SplitFusionPath, 'lite.control.job.1.R', sep=''))
}
        source(paste(SplitFusionPath, 'lite.control.job.2.R', sep=''))


## postscript summary
## example reads
system('cat */*.fusion.table.xls | head -n1 | cut -f1-5 > fusion.brief.summary.txt')
system('cat */*.brief.summary | grep -v num_unique_reads | cut -f1-5 >> fusion.brief.summary.txt')

## filter
(fu.filter = readLines(paste(DEPATH, '/fusion.filters.txt', sep='')))
fu.filter

fu0 = read.table('fusion.brief.summary.txt', sep='\t', header=T, stringsAsFactors=F)
head(fu0)

fu0$ge1 = sub('---.*', '', fu0$GeneExon5...GeneExon3)
fu0$ge2 = sub('.*---', '', fu0$GeneExon5...GeneExon3)
fu0$g1 = sub('_.*', '', fu0$ge1)
fu0$g2 = sub('_.*', '', fu0$ge2)
fu0$g1g2 = paste(fu0$g1, fu0$g2, sep = '---')

fu1 = subset(fu0, !(g1g2 %in% fu.filter))
table(fu1$g1g2)
head(fu1)
fu2 = subset(fu1, select=c(AP7, GeneExon5...GeneExon3, num_unique_reads, frame))

iiFreq = read.table('../iiFreq.txt', header=T, stringsAsFactors=F)
head(iiFreq)

fu3 = merge(iiFreq[, c('AP7', 'Sample_ID', 'RNA2DNA')], fu2, by='AP7', all=TRUE)
	fu3[is.na(fu3)] = '-'
	head(fu3)
write.table(fu3, 'fusion.summary.txt', sep='\t', col.names=T, row.names=F)

system(paste("ssconvert fusion.summary.txt ../", runID, "_Fusion_Summary.xls 2>/dev/null", sep=''))

system('mv */*---* .')

### DONE

