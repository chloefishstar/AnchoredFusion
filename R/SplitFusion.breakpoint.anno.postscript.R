#' SplitFusion.breakpoint.anno.postscript
#'
#' Post annotation of breakpoint detected by SplitFusion.
#'
#' @param runInfo config file which sets the path and parameters of depended tools.
#'
#' @return NULL
#' @export
#'
#' @examples
#'
#' #"anno.left.right" file outputted by last step "breakpoint-anno" of SplitFusion should be exist.
#'
#' SplitFusion.breakpoint.anno.postscript(runInfo = "/path/example.runInfo")
#'
SplitFusion.breakpoint.anno.postscript <- function(runInfo){
# breakpoint post annotation processing
#
# required:
#	- fu.anno.right
#	- fu.anno.left
#	- fu.breakpointID.stgLeft
#	- fu.breakpointID.stgRight
#
options(width=204)
options(scipen=999)
#source('../../run.info')
#.libPaths(c(paste(DEPATH, '/Rlib', sep=''), .libPaths()))

#runInfo = commandArgs(TRUE)[1]  ###Modified by Baifeng###
source(runInfo)  ###Modified by Baifeng###

library(plyr)
subii = sub('.*/', '', getwd())

## target genes
ii = read.table(sampleInfo, sep='\t', stringsAsFactors=F, header=T)
head(ii)
(panel = sub("[vV][0-9]", "", ii$Panel[ii$AP7==subii][1]))
#(genes = readLines(paste(DEPATH, '/panel/', panel, '.target.genes', sep='')))
#(genes = readLines(paste(SplitFusionPath, '/scripts/', panel, '.target.genes', sep=''))) ###Modified by Baifeng###
#if (file.exists(paste(Panel_path, '/', panel,'.target.genes.txt', sep=''))){
#	(genes = readLines(paste(Panel_path, '/', panel,'.target.genes.txt', sep=''))) ###R version###
#}

##====1: left-right annotations
    colnames = c('readID', 'chrorp_L', 'chr_L', 'orp_L', 'pos_L', 'strand_L'
		, 'num_unique_molecules', 'num_start_site', 'num_start_site2', 'breakpoint', 'overlap'
		, 'gene_L', 'geneStrand_L', 'inEx_L', 'functiontype_L', 'nm_L', 'exon_L', 'cdna_L'
		, 'chrorp_R', 'chr_R', 'orp_R', 'pos_R', 'strand_R', 'gene_R', 'geneStrand_R', 'inEx_R', 'functiontype_R', 'nm_R', 'exon_R', 'cdna_R')
    lr1 = read.table('anno.left.right', sep=' ', header=F, fill=T, stringsAsFactors=F, col.names=colnames)
    lr1$cdna_L = suppressWarnings(as.numeric(lr1$cdna_L))
    lr1$cdna_R = suppressWarnings(as.numeric(lr1$cdna_R))
    head(lr1)
    (n.lr1 = nrow(lr1))
    lr1$exonn_L = suppressWarnings(as.numeric(sub('exon|intron','',lr1$exon_L)))
    lr1$exonn_R = suppressWarnings(as.numeric(sub('exon|intron','',lr1$exon_R)))

    ###add transcript orientation
    orien <- data.frame(fread(paste0(Panel_path,'/ENSEMBL.orientation.txt')))
    rownames(orien) <- orien[,"ID"]
    nm_L <- lr1[,"nm_L"]
    lr1[,"geneStrand_L"] <- orien[nm_L,"Orientation"]
    nm_R <- lr1[,"nm_R"]
    lr1[,"geneStrand_R"] <- orien[nm_R,"Orientation"]

#    lr1.rev = subset(lr1, (strand_L != geneStrand_L & strand_R != geneStrand_R)
#                       | (is.na(geneStrand_L) & strand_R != geneStrand_R)
#                       | (is.na(geneStrand_R) & strand_L != geneStrand_L))
#    lr1.nor = subset(lr1, !((strand_L != geneStrand_L & strand_R != geneStrand_R)
#                       | (is.na(geneStrand_L) & strand_R != geneStrand_R)
#                       | (is.na(geneStrand_R) & strand_L != geneStrand_L)))

#    if (nrow(lr1.rev)>0){
#        lr1.rev.r = lr1.rev
#        tmp.n = gsub('_L', '_Tmp000R', names(lr1.rev))
#        tmp.n2 = gsub('_R', '_Tmp000L', tmp.n)
#        tmp.n3 = gsub('_Tmp000', '_', tmp.n2)
#        names(lr1.rev.r) = tmp.n3

#        lr1 = rbind(lr1.nor, lr1.rev.r)
#   } else {lr1 = lr1}



##====2: filter targeted site (non-ligation end) being non-targets
if (n.lr1 >0){
	# extract ligation end
	lr1$ligEnd0 = sub('.*umi:C', '', lr1$readID)
	lr1$ligEndChr = sub('P.*', '', lr1$ligEnd0)
	lr1$ligEndPos = as.numeric(sub('-.*', '', sub('.*P', '', lr1$ligEnd0))) - 100000000

	lr1$gene_T = ifelse((lr1$chr_L == lr1$ligEndChr) & (abs(lr1$ligEndPos - as.numeric(lr1$pos_L)) < 300)
				, lr1$gene_R
				, ifelse(lr1$chr_R == lr1$ligEndChr & (abs(lr1$ligEndPos - as.numeric(lr1$pos_R)) < 300)
                        	, lr1$gene_L
				, "-"
			))
    lr1$exclude = 0

    if (file.exists(paste(Panel_path, '/', panel,'.target.genes.txt', sep=''))){
        (genes = readLines(paste(Panel_path, '/', panel,'.target.genes.txt', sep=''))) ###R version###
	lr1$exclude[!(lr1$gene_T %in% genes)] = 1
    }

	table(lr1$exclude)
    lr2 = subset(lr1, exclude==0)
	}else{
	    lr2 = lr1
}

##== remove recurrent breakpoints that are not significant (mannually curated)
#filter.breakpointID = readLines(paste(DEPATH, '/fusion.filter.breakpointID', sep=''))
#lr = subset(lr, !(fusionID %in% filter.breakpointID))

## kepp original cdna pos for later frameness calculation
lr2$cdna_L0 = lr2$cdna_L
lr2$cdna_R0 = lr2$cdna_R

## lmr
## 	correct breakpoint exon number, cdna, gdna by add/minus middle split size
if (file.exists('mid.anno2')){
    colnamesM = c('annoPos', 'readID', 'chr_M', 'start_M', 'end_M', 'gene_M', 'geneStrand_M', 'inEx_M', 'functiontype_M', 'nm_M', 'exon_M', 'cdna_M')
    mid = read.table('mid.anno2', sep=' ', header=F, fill=T, stringsAsFactors=F, col.names=colnamesM)
	mid$exonn_M = suppressWarnings(as.numeric(sub('exon|intron', '', mid$exon_M)))
	head(mid)

	lmr = merge(lr2, mid, by="readID")
	if (nrow(lmr) >0){
		lmr_L = subset(lmr, nm_L==nm_M)
		lmr_R = subset(lmr, nm_R==nm_M & nm_L != nm_M)

		# mid belong to Left
		if (nrow(lmr_L)>0){
		    	lmr_L$exon_L = lmr_L$exon_M
		    	lmr_L$exonn_L = lmr_L$exonn_M
			lmr_L$cdna_L = lmr_L$cdna_L + (1 - (lmr_L$exonn_L > lmr_L$exonn_M)*2) * lmr_L$overlap

			lmr_L$absMstart = abs(lmr_L$start_M - lmr_L$orp_L)
			lmr_L$absMend = abs(lmr_L$end_M - lmr_L$orp_L)
			lmr_L$pos_M = ifelse(
				lmr_L$absMstart > lmr_L$absMend
				, lmr_L$start_M
				, lmr_L$end_M
				)
			lmr_L$chrpos_M = paste(lmr_L$chr_L, lmr_L$pos_M, sep='_')
			lmr_L$chrpos_Target = paste(lmr_L$chr_R, lmr_L$pos_R, sep='_')
			lmr_L$breakpoint = ifelse(
				lmr_L$chrpos_M < lmr_L$chrpos_Target
				, paste(lmr_L$chrpos_M, lmr_L$chrpos_Target, sep='__')
				, paste(lmr_L$chrpos_Target, lmr_L$chrpos_M, sep='__')
				)
			}

		# mid belong to Right
		if (nrow(lmr_R)>0){
		    	lmr_R$exon_R = lmr_R$exon_M
		    	lmr_R$exonn_R = lmr_R$exonn_M
			lmr_R$cdna_R = lmr_R$cdna_R + (1 - (lmr_R$exonn_R > lmr_R$exonn_M)*2) * lmr_R$overlap

			lmr_R$absMstart = abs(lmr_R$start_M - as.numeric(lmr_R$orp_R))
			lmr_R$absMend = abs(lmr_R$end_M - as.numeric(lmr_R$orp_R))
			lmr_R$pos_M = ifelse(
				lmr_R$absMstart > lmr_R$absMend
				, lmr_R$start_M
				, lmr_R$end_M
				)
			lmr_R$chrpos_M = paste(lmr_R$chr_R, lmr_R$pos_M, sep='_')
			lmr_R$chrpos_Target = paste(lmr_R$chr_L, lmr_R$pos_L, sep='_')
			lmr_R$breakpoint = ifelse(
				lmr_R$chrpos_M < lmr_R$chrpos_Target
				, paste(lmr_R$chrpos_M, lmr_R$chrpos_Target, sep='__')
				, paste(lmr_R$chrpos_Target, lmr_R$chrpos_M, sep='__')
				)
			}

		lmr2 = rbind(lmr_L, lmr_R)
		lmr2.keep = lmr2[, c(names(lr2))]
		    lr0 = lr2[!(lr2$readID %in% lmr2.keep$readID),]
		    nrow(lr0); nrow(lmr2.keep); nrow(lr2)
		lr2b = rbind(lmr2.keep, lr0)
		} else {
			lr2b = lr2
		}
	} else {
		lr2b = lr2
}

## make Left as 5' and Right as 3'
# tab.sen = ddply(lr2b, .(strand_L, geneStrand_L, strand_R, geneStrand_R), summarize, n=length(readID))
#sense = subset(lr2b, (strand_L == geneStrand_L & strand_R == geneStrand_R)
#			| (is.na(geneStrand_L) & strand_R == geneStrand_R)
#			| (is.na(geneStrand_R) & strand_L == geneStrand_L)
#)
#antisense = subset(lr2b, (strand_L != geneStrand_L & strand_R != geneStrand_R)
#			| (is.na(geneStrand_L) & strand_R != geneStrand_R)
#			| (is.na(geneStrand_R) & strand_L != geneStrand_L)
#)
#nosense = subset(lr2b, is.na(geneStrand_L) & is.na(geneStrand_R)
#			| (strand_L == geneStrand_L & strand_R != geneStrand_R)
#			| (strand_L != geneStrand_L & strand_R == geneStrand_R)
#			)
#nrow(sense); nrow(antisense); nrow(nosense); nrow(lr2b)
#lr2c = rbind(sense, anti2, nosense)
	#nrow(lr2b); nrow(lr2c)

# reverse strand for sense Read2
r1 = lr2b[!(grepl("/2", lr2b$readID)),]
r2 = lr2b[grepl("/2", lr2b$readID),]


## switch L/R for Read2
if (nrow(r2)>0){
    r2r = r2
	tmp.n = gsub('_L', '_Tmp000R', names(r2))
	tmp.n2 = gsub('_R', '_Tmp000L', tmp.n)
	tmp.n3 = gsub('_Tmp000', '_', tmp.n2)
	names(r2r) = tmp.n3

    lr3 = rbind(r1, r2r)
} else {lr3 = lr2b}



n.lr3 = nrow(lr3); n.lr3


#=========================================================
#=========================================================
if (n.lr3 >0){
	## add the 6 base back to cDNA position
	lr3$cdna_L = lr3$cdna_L0 + 6
	lr3$cdna_R = lr3$cdna_R0 - 6

	lr3.rev = subset(lr3, ((((strand_L != geneStrand_L & strand_R != geneStrand_R)
                       | (is.na(geneStrand_L) & strand_R != geneStrand_R)
                       | (is.na(geneStrand_R) & strand_L != geneStrand_L)) & grepl("/1",readID))
		       | (((strand_L == geneStrand_L & strand_R == geneStrand_R)
                       | (is.na(geneStrand_L) & strand_R == geneStrand_R)
                       | (is.na(geneStrand_R) & strand_L == geneStrand_L)) & grepl("/2",readID))))
        lr3.nor = subset(lr3, !(((((strand_L != geneStrand_L & strand_R != geneStrand_R)
                       | (is.na(geneStrand_L) & strand_R != geneStrand_R)
                       | (is.na(geneStrand_R) & strand_L != geneStrand_L)) & grepl("/1",readID))
                       | (((strand_L == geneStrand_L & strand_R == geneStrand_R)
                       | (is.na(geneStrand_L) & strand_R == geneStrand_R)
                       | (is.na(geneStrand_R) & strand_L == geneStrand_L)) & grepl("/2",readID)))))

    if (nrow(lr3.rev)>0){
        lr3.rev.r = lr3.rev
        tmp.n = gsub('_L', '_Tmp000R', names(lr3.rev))
        tmp.n2 = gsub('_R', '_Tmp000L', tmp.n)
        tmp.n3 = gsub('_Tmp000', '_', tmp.n2)
        names(lr3.rev.r) = tmp.n3

        lr3 = rbind(lr3.nor, lr3.rev.r)
   } else {lr3 = lr3}

	## intra-gene
	lr3$intragene = ifelse(lr3$gene_L == lr3$gene_R, 1, 0)
	lr3$exon = substr(lr3$exon_L,1,4) ### NULL?
	lr3$exonD = 0
	lr3$exonD = ifelse(lr3$gene_L == lr3$gene_R, lr3$exonn_R - lr3$exonn_L, lr3$exonD)
	lr3$gec_L = paste(lr3$gene_L, ' ', lr3$exon_L, ' c.', lr3$cdna_L, ' (', lr3$nm_L, ')', sep='')
	lr3$gec_R = paste(lr3$gene_R, ' ', lr3$exon_R, ' c.', lr3$cdna_R, ' (', lr3$nm_R, ')', sep='')
	lr3$gec_LR = paste(lr3$gec_L, lr3$gec_R, sep='---')
	lr3$ge1 = paste(lr3$gene_L, '_', lr3$exon_L, sep='')
	lr3$ge2 = paste(lr3$gene_R, '_', lr3$exon_R, sep='')
	lr3$ge1ge2 = paste(lr3$ge1,'---', lr3$ge2, sep='')


	##======================
	## frame status
	##======================
	    lr3$frameD = (lr3$cdna_R0 %% 3 - (lr3$cdna_L0 + abs(lr3$overlap)) %% 3) %% 3
	    lr3$frame = ifelse(lr3$frameD ==1, 'in-frame', 'out-frame')
	#     lr3$frameD = (lr3$cdna_R0 + lr3$cdna_L0 - lr3$overlap) %% 3
         #    lr3$frame = ifelse(lr3$frameD == 0, 'in-frame', 'out-frame')
	    lr3$frame[is.na(lr3$frame)] = 'N.A.'
		#head(lr3)
		## expect mostly in-frame or _NA_. Otherwise, might need inspection
		#table(lr3$frame,useNA='always')

	    lr3$neighb = ifelse(lr3$exonD==1, 1, 0)
	    lr3$neighb[is.na(lr3$neighb)] = 0

		## order
		lr3 = lr3[order(lr3$intragene, lr3$frame, lr3$gec_LR),]

	##======================
	## fitler num_start_site: >=3 or known recurrent

#	(known.partners = readLines(paste(DEPATH, '/fusion.partners.txt', sep='')))
#	(known.ge = readLines(paste(DEPATH, '/fusion.gene-exon.txt', sep='')))
#	(known.filter = readLines(paste(DEPATH, '/fusion.gene-exon.filter.txt', sep='')))
##(known.partners = readLines(paste(SplitFusionPath, '/scripts/','/fusion.partners.txt', sep='')))   ###Modified by Baifeng###
##(known.ge = readLines(paste(SplitFusionPath, '/scripts/','/fusion.gene-exon.txt', sep='')))   ###Modified by Baifeng###
##(known.filter = readLines(paste(SplitFusionPath,'/scripts/','/fusion.gene-exon.filter.txt', sep='')))
(known.partners = readLines(paste(Panel_path,'/fusion.partners.txt', sep='')))   ###R version###
(known.ge = readLines(paste(Panel_path,'/fusion.gene-exon.txt', sep='')))   ###R version###
(known.filter = readLines(paste(Panel_path,'/fusion.gene-exon.filter.txt', sep=''))) ###R version###


	lr3$known=0
	lr3$known[ (lr3$gene_T == lr3$gene_L & lr3$gene_R %in% known.partners & lr3$intragene==0)
		  |(lr3$gene_T == lr3$gene_R & lr3$gene_L %in% known.partners & lr3$intragene==0)
		  | lr3$ge1ge2 %in% known.ge
		  | (lr3$ge1 %in% c('FGFR1_exon17', 'FGFR2_exon17', 'FGFR3_exon17') & lr3$intragene==0)
			] = 1
	lr32 = subset(lr3, (known==1 | (num_start_site >=3 & intragene==0)) & !(ge1ge2 %in% known.filter))
		nrow(lr3); nrow(lr32)

	if (nrow(lr32)>0){
		# output for furture research
		write.table(lr32, paste(subii, '.fusion.list.pre-processing.txt', sep=''), row.names=F, quote=F, sep='\t')

	##======================
	## Further processing
	##======================
	## To remove:
	#	neighbour exons
	#	out-frame (keep out-frame of different genes. FGFRs fusion can be functional even out-frame)
	#	same exon
		lr40 = subset(lr32, !(neighb==1 | (intragene==1 & exon_L==exon_R)
					| (intragene==1 & frame=='out-frame')
					))
		lr4 = lr40[!duplicated(lr40$readID),]
		(n.lr40 = nrow(lr40))
		(n.lr4 = nrow(lr4))

		if (n.lr4>0){
		lr4$AP7 = subii
		# remove illegal symbol
		lr4$gec_LR = gsub('\\(', '.', lr4$gec_LR)
		lr4$gec_LR = gsub('\\)', '.', lr4$gec_LR)
		lr4$ge1ge2 = gsub('\\(', '.', lr4$ge1ge2)
		lr4$ge1ge2 = gsub('\\)', '.', lr4$ge1ge2)

	## output
		lr4 = lr4[order(-lr4$known, -lr4$num_start_site, lr4$intragene, lr4$frame),]
		fusion.list = subset(lr4, select=c('AP7', 'ge1ge2', 'frame','gec_LR','chr_L','pos_L','inEx_L','gene_L','nm_L','exon_L','cdna_L'
					, 'overlap','chr_R','pos_R','inEx_R','gene_R','nm_R','exon_R','cdna_R','intragene','readID'
					, 'breakpoint', 'num_start_site', 'known'))
			write.table(fusion.list, paste(subii, '.fusion.list.post-processing.txt', sep=''), row.names=F, quote=F, sep='\t')

		## consolidate breakpoint
			fusion.list[is.na(fusion.list)]='-'
			fusion.list$readID2 = sub(':umi:.*', '', fusion.list$readID)
		fusion.tab1 = ddply(fusion.list, .(breakpoint), here(summarize)
						, 'AP7'=subii
						, 'ge1ge2'=names(sort(table(ge1ge2),decreasing=TRUE)[1])
						, 'num_unique_reads'=length(unique(readID2))
						, 'frame'=frame[1]
						, 'transcript_L'=nm_L[1], 'transcript_R'=nm_R[1]
						, 'function_L'=inEx_L[1], 'function_R'=inEx_R[1]
						, 'cdna_L'=cdna_L[1], 'cdna_R'=cdna_R[1]
						, 'intragene'=intragene[1]
						, 'overlap'=overlap[1]
						, 'num_start_site'=num_start_site[1]
						, 'known'=known[1]
						)
		fusion.tab1 = fusion.tab1[order(-fusion.tab1$known, -fusion.tab1$num_start_site, fusion.tab1$intragene, fusion.tab1$frame, -fusion.tab1$num_unique_reads),]

		##==============================
		## filter: overlap length and num of reads
		##==============================
		# Deprecated: fitler those with too few unique reads given long L/R overlap, when overlap >6
			# fusion.tab2 = subset(fusion.tab1, ((num_unique_reads > overlap) | overlap <=6))
			fusion.tab2=fusion.tab1
				head(fusion.tab2)
			out.names = c("AP7", "ge1ge2", "frame", "num_start_site", "num_unique_reads", "breakpoint"
					, "transcript_L", "transcript_R", "function_L", "function_R", "cdna_L", "cdna_R", "intragene", "known")

		fusion.table = fusion.tab2[, out.names]
			fusion.table = fusion.table[order(-fusion.table$known, fusion.table$frame, fusion.table$intragene, -fusion.table$num_start_site),]
			name1 = names(fusion.table)
			name2 = sub('_L', "_5", name1)
			name3 = sub('_R', "_3", name2)
			name4 = sub('ge1ge2', "GeneExon5_GeneExon3", name3)
			names(fusion.table) = name4
		write.table(fusion.table, paste(subii, '.fusion.table.txt', sep=''), row.names=F, quote=F, sep='\t')

		##=======================================
		## export max 10 example fusion reads
		##=======================================
		    # inter-gene, in-frame
		    # or known
		ex = subset(fusion.table, known ==1
					| (intragene==0 & frame =='in-frame' & num_start_site >=3)
				)[, 1:14]
		# keep 1st
		ex1 = ex[order(ex$frame, -ex$num_start_site, -ex$num_unique_reads),]
		ex2 = ex1[!duplicated(ex1$"GeneExon5_GeneExon3"),]

			write.table(ex2, paste(subii, '.brief.summary', sep=''), row.names=F, quote=F, sep='\t')


	ex2 = read.table(paste(subii, '.brief.summary', sep=''), sep='\t', stringsAsFactors=F, header=T)
		(nex = min(10, nrow(ex2)))
			if (nex >0){
			## show max 25 example reads
				for (i in 1:nex){
					## get read ID
					(gei = ex2$"GeneExon5_GeneExon3"[i])
					(ngeci = ex2$num_unique_reads[i])
					readids = unique(lr4$readID[lr4$ge1ge2 == gei])
					(sample.n = min(10, ngeci, length(readids)))
					(readid = sample(unique(lr4$readID[lr4$ge1ge2 == gei]), sample.n))
					writeLines(readid, 'tmp.readid')
					system('sed -e "s/:umi.*//" -e "s:/1.*::" -e "s:/2.*::" tmp.readid | sort -u > tmp.readid2')

				## get fa
					if (file.exists(paste("../../fq.demux/", subii, ".R2.fq", sep=''))){
#						system(paste("grep -f tmp.readid2 -A1 ../../fq.demux/", subii, ".R2.fq | sed 's: :/1 :' | sed 's:^@:>:' | grep -v '\\-\\-' > ", subii, '.', gei, ".txt", sep=''))
					} else if (file.exists(paste("../../fq.demux/", subii, ".R2.fq.gz", sep=''))){
#						system(paste("zcat ../../fq.demux/", subii, ".R2.fq.gz | grep -f tmp.readid2 -A1 | sed 's: :/1 :' | sed 's:^@:>:' | grep -v '\\-\\-' > ", subii, '.', gei, ".txt", sep=''))
					} else if (file.exists(paste(bam_path, subii, ".consolidated.bam", sep=''))){
#						system(paste(REPPATH, "/samtools/bin/samtools view ../../bam/", subii
						system(paste(samtools, " view ", bam_path, subii   ###Modified by Baifeng###
#							, ".consolidated.bam | grep -f tmp.readid2 | cut -f1,10 | sed 's/^/>/' | tr '\t' '\n' > ", subii, '.', gei, ".txt", sep=''))
							, ".consolidated.bam | grep -f tmp.readid2 | cut -f1,2,10 > ", subii, '.', gei, ".txt", sep=''))
						system(paste(samtools, " view ", bam_path, subii   ###Modified by Baifeng###
#                                                       , ".consolidated.bam | grep -f tmp.readid2 | cut -f1,10 | sed 's/^/>/' | tr '\t' '\n' > ", subii, '.', gei, ".txt", sep=''))
                                                        , ".consolidated.bam | grep -f tmp.readid2 | ",samtools, " view -T ", hgRef, " -bS - > ", subii, '.', gei, ".bam", sep=''))
						fq2png(seq=paste0(subii,".",gei,".txt"),aln="breakpoint.reads",nam=paste0(subii,".",gei))
					}
					system('rm tmp.readid*')
				}
	       		}
		}
	}
}

	if (!file.exists(paste(subii, '.brief.summary', sep=''))){
			file.create(paste(subii, '.brief.summary', sep=''))

}
}
