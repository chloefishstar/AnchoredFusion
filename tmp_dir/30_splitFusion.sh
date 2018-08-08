#!/bin/bash

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
. ../run.info.sh
subii='SUBID_'
mkdir -p $subii
cd  $subii

##==== 1: consolidated.bam to breakpoint candidates ====
if [ ! -f breakpoint ]; then
	bash $pipelinePATH/scripts/SplitFusion.consolidated-breakpoint.sh
fi

##==== 2: filters ====
if [ ! -f fusion.candidates ]; then
	bash $pipelinePATH/scripts/SplitFusion.breakpoint-filter.sh
fi

##==== 3: Annotate breakpoint gene, exon, cDNA position
if [ ! -f anno.left.right ]; then
	bash $pipelinePATH/scripts/SplitFusion.breakpoint-anno.sh
fi

##==== 4: further processing (in-frame status determination, etc)
	Rscript $pipelinePATH/scripts/SplitFusion.breakpoint-anno.postscript.R

rm _*

rm ../_job_.$subii

