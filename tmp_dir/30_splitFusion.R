
###################
## Structure variation detection, including fusion/cDNA/large indel
##	-- based on BWA MEM SA reads
###########################
## To run as a plugin independent of the whole pipeline:
## 1). bsub -Is -q interact /bin/bash
## 2). cd to analysis run folder
## 3). if control SA needs to be filtered, copy the filter file to 'bam' folder or, if not applied to all or different filter files for different samples, bam/A##-P7## folder and renamed as 'sa.filter'
##	cp ~/projects/criStr/rd2.filter bam/sa.filter  ## or
##	cp ~/projects/criStr/A01.filter bam/A01-P701/sa.filter
## then:
## . ./run.info.sh; Rscript --vanilla $pipelinePATH/30_splitFusion.R > log.30_splitFusion.Rout 2>log.30_splitFusion.ERROR.Rout &

source('run.info')
options(width=204)

ii = read.table("iiFreq.txt",header=T, sep='\t', stringsAsFactors=F, fill=T)
head(ii)
dir.create('fusion', showWarnings = FALSE, recursive = TRUE)
setwd('fusion')

## prep jobs
(i=ii$AP7[1])
for (i in ii$AP7){
	system(paste("cat ", pipelinePATH, "/30_splitFusion.sh | sed 's/SUBID_/", i,"/g' > _job_.", i, sep=''))
        source(paste(pipelinePATH, '/lite.control.job.1.R', sep=''))
}
        source(paste(pipelinePATH, '/lite.control.job.2.R', sep=''))

## fusion calls
system('cat */*.brief.summary | head -n1 | cut -f1-12 > fusion.brief.summary.txt')
system('cat */*.brief.summary | grep -v AP7 | cut -f1-12 >> fusion.brief.summary.txt')
fu0 = read.table('fusion.brief.summary.txt', sep='\t', header=T, stringsAsFactors=F)
    head(fu0)

## sample info
iiFreq = read.table('../iiFreq.txt', header=T, stringsAsFactors=F)
	head(iiFreq)
	ii.names = names(iiFreq)
	if (!('RNA2DNA' %in% ii.names)){
		stop("Please calculate RNA2DNA first.")
	}

iiFreq2=iiFreq[, colnames(iiFreq) %in% c('AP7', 'Sample_ID', 'RNA2DNA')]
	iiFreq2$RNA_QC='-'
	iiFreq2$RNA_QC[iiFreq2$RNA2DNA>0.1]='PASS'
	iiFreq2$RNA_QC[iiFreq2$RNA2DNA<0.1]='FAILED'
	iiFreq2$RNA_QC[iiFreq2$RNA2DNA=='-']='N.A.'

## merge
fu.n = as.numeric(system('wc -l fusion.brief.summary.txt | sed "s/ .*//"', intern=T))
if (fu.n > 0) {
    fu2 = merge(iiFreq2, fu0, by='AP7', all=TRUE)
	fu2[is.na(fu2)] = '-'
	head(fu2)
} else {
    fu2 = iiFreq2
	fu2$GeneExon5...GeneExon3 = '-'
	fu2$num_unique_reads = '-'
	fu2$frame = '-'
}

## remove redundant name AP7
if (unique(fu2$AP7 == fu2$Sample_ID) == TRUE){
	fu2 = fu2[, !(names(fu2) %in% 'AP7')]
}

write.table(fu2, 'fusion.summary.txt', sep='\t', col.names=T, row.names=F)
system(paste("ssconvert fusion.summary.txt ../", runID, "_Fusion_Summary.xls 2>/dev/null", sep=''))
system(paste("cp ../*_Fusion_Summary.xls ", outputRunFull, "/", sep=''))
system('mv */*---* .')

### DONE
