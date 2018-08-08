# breakpoint post annotation processing
# 
# required: 
#	- fu.anno.right
#	- fu.anno.left
#	- fu.breakpointID.stgLeft
#	- fu.breakpointID.stgRight
#
options(width=204)
source('../../run.info')
.libPaths(c(paste(DEPATH, '/Rlib', sep=''), .libPaths()))
library(plyr)
subii = sub('.*/', '', getwd())

## target genes
ii = read.table('../../iiFreq.txt', sep='\t', stringsAsFactors=F, header=T)
head(ii)
(panel = sub("[vV][0-9]", "", ii$Panel[ii$AP7==subii][1]))
(genes = readLines(paste(DEPATH, '/panel/', panel, '.target.genes', sep='')))

##====1: left-right annotations
    colnames = c('readID', 'chrorp_L', 'chr_L', 'orp_L', 'pos_L', 'strand_L', 'overlap', 'fusion_point', 'num_start_site', 'stgSide'
		, 'gene_L', 'geneStrand_L', 'inEx_L', 'functiontype_L', 'nm_L', 'exon_L', 'cdna_L'
		, 'chrorp_R', 'chr_R', 'orp_R', 'pos_R', 'strand_R', 'gene_R', 'geneStrand_R', 'inEx_R', 'functiontype_R', 'nm_R', 'exon_R', 'cdna_R')
    lr1 = read.table('anno.left.right', sep=' ', header=F, fill=T, stringsAsFactors=F, col.names=colnames)
    lr1$cdna_L = as.numeric(lr1$cdna_L)
    lr1$cdna_R = as.numeric(lr1$cdna_R)
    head(lr1)
    (n.lr1 = nrow(lr1))
    lr1$exonn_L = suppressWarnings(as.numeric(sub('exon','',lr1$exon_L)))
    lr1$exonn_R = suppressWarnings(as.numeric(sub('exon','',lr1$exon_R)))

##====2: filter
## Among one-sided stagger fusion, remove fusion candidates with non-staggering side not on target genes
if (n.lr1 >0){
	    lr1$exclude = 0
		#==== filter non-stg side being non-targets
	    lr1$exclude[     (lr1$stgSide =='stgLeft' & !(lr1$gene_R %in% genes))
			    | (lr1$stgSide =='stgRight' & !(lr1$gene_L %in% genes))
			    ] = 1
		table(lr1$exclude)
	    lr2 = subset(lr1, exclude==0)
	}else{
	    lr2 = lr1
}

## exlude off targets (neither of the partners is on target)
lr = subset(lr2, (lr2$gene_L %in% genes) | (lr2$gene_R %in% genes))


##== remove recurrent breakpoints that are not significant (mannually curated)
#filter.breakpointID = readLines(paste(DEPATH, '/fusion.filter.breakpointID', sep=''))
#lr = subset(lr, !(fusionID %in% filter.breakpointID))

## kepp original cdna pos for later frameness calculation
lr$cdna_L0 = lr$cdna_L
lr$cdna_R0 = lr$cdna_R
head(lr,3)

## lmr
if (file.exists('mid.anno.ext')){
	    colnamesM = c('readID', 'start_M', 'end_M', 'mq_M', 'gene_M', 'geneStrand_M', 'inEx_M', 'functiontype_M', 'nm_M', 'exon_M', 'cdna_M')
	    mid = read.table('mid.anno2', sep=' ', header=F, fill=T, stringsAsFactors=F, col.names=colnamesM)
		mid2 = subset(mid, inEx_M == 'exonic')
		head(mid2)
		#head(subset(mid2, gene_M ==chk), 6)

		lmr = merge(lr, mid2, by="readID")
		if (nrow(lmr) >0){
			lmr$exonn_M = as.numeric(sub('exon', '', lmr$exon_M))
			head(lmr,4)
			#head(subset(lmr, gene_M ==chk), 6)

		    ## correct breakpoint exon number by add/minus middle split size
		    lmr$nm_L[is.na(lmr$nm_L)] = '-l'
		    lmr$nm_R[is.na(lmr$nm_R)] = '-r'
		    lmr$exon_L[lmr$nm_L == lmr$nm_M] = lmr$exon_M[lmr$nm_L == lmr$nm_M]
		    lmr$exon_R[lmr$nm_R == lmr$nm_M] = lmr$exon_M[lmr$nm_R == lmr$nm_M]

		    ## correct breakpoint cdna position by add/minus middle split size
		    lmr$cdna_L[lmr$nm_L==lmr$nm_M] = lmr$cdna_L[lmr$nm_L == lmr$nm_M] + (1 - (lmr$exonn_M[lmr$nm_L == lmr$nm_M] > lmr$exonn_L[lmr$nm_L == lmr$nm_M])*2) * lmr$overlap_L[lmr$nm_L == lmr$nm_M]
		    lmr$cdna_R[lmr$nm_R==lmr$nm_M] = lmr$cdna_R[lmr$nm_R == lmr$nm_M] + (1 - (lmr$exonn_M[lmr$nm_R == lmr$nm_M] > lmr$exonn_R[lmr$nm_R == lmr$nm_M])*2) * lmr$overlap_R[lmr$nm_R == lmr$nm_M]

		    lmr$exonn_L = as.numeric(sub('exon', '', lmr$exon_L))
		    lmr$exonn_R = as.numeric(sub('exon', '', lmr$exon_R))
		    
		    lmr.keep = lmr[, c(names(lr))]
		    lr0 = lr[!(lr$readID %in% lmr.keep$readID),]

		    nrow(lr0); nrow(lmr.keep); nrow(lr)
		    lr2 = rbind(lmr.keep, lr0)
		} else {
			lr2 = lr
		}
	} else {
		lr2 = lr
}

## make Left as 5' and Right as 3'
sense = subset(lr2, (strand_L == geneStrand_L & strand_R == geneStrand_R)
			| (is.na(geneStrand_L) & strand_R == geneStrand_R)
			| (is.na(geneStrand_R) & strand_L == geneStrand_L)
)
antisense = subset(lr2, (strand_L != geneStrand_L & strand_R != geneStrand_R)
			| (is.na(geneStrand_L) & strand_R != geneStrand_R)
			| (is.na(geneStrand_R) & strand_L != geneStrand_L)
)
nosense = subset(lr2, is.na( geneStrand_L) & is.na(geneStrand_R))
#head(nosense[nosense$gene_L =='PIK3CA',])
#head(sense)
#head(nosense)

nrow(sense); nrow(antisense); nrow(nosense); nrow(lr2)

## switch L/R
    tmp.n = gsub('_L', '_000R', names(antisense))
    tmp.n2 = gsub('_R', '_000L', tmp.n)
    tmp.n3 = gsub('_000', '_', tmp.n2)
    tmp.n3
    anti2 = antisense
    names(anti2) = tmp.n3
    #head(antisense,3)
    #head(anti2,3)

lr3 = rbind(sense, anti2, nosense)
lr3 = lr3[!duplicated(lr3$readID),]
nrow(lr2)
n.lr3 = nrow(lr3)

#=========================================================
#=========================================================
if (n.lr3 >0){
	## add the 6 base back to cDNA position
	lr3$cdna_L = lr3$cdna_L0 + 6
	lr3$cdna_R = lr3$cdna_R0 - 6

	## intra-gene
	lr3$intragene = ifelse(lr3$gene_L == lr3$gene_R, 1, 0)
	lr3$exon = substr(lr3$exon_L,1,4)
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
	
	(known.partners = readLines(paste(DEPATH, '/fusion.partners.txt', sep='')))
	(known.ge = readLines(paste(DEPATH, '/fusion.gene-exon.txt', sep='')))
	#(known.skipping = readLines(paste(DEPATH, '/fusion.skipping.txt', sep='')))
	lr3$known=0
	lr3$known[ (lr3$stgSide =='stgLeft'  & lr3$gene_L %in% known.partners & lr3$intragene==0)
		  |(lr3$stgSide =='stgRight' & lr3$gene_R %in% known.partners & lr3$intragene==0)
		  | lr3$ge1ge2 %in% known.ge
		#  | (lr3$ge1 %in% known.skipping & lr3$neighb==0 & lr3$intragene==1)
		#  | (lr3$ge2 %in% known.skipping & lr3$neighb==0 & lr3$intragene==1)
		  | (lr3$ge1 %in% c('FGFR1_exon17', 'FGFR2_exon17', 'FGFR3_exon17') & lr3$intragene==0)
			] = 1
	lr32 = subset(lr3, known==1 | (num_start_site >=3 & intragene==0))
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
		lr4$gec_LR = gsub('\\(', '.', lr4$gec_LR)
		lr4$gec_LR = gsub('\\)', '.', lr4$gec_LR)
		lr4$ge1ge2 = gsub('\\(', '.', lr4$ge1ge2)
		lr4$ge1ge2 = gsub('\\)', '.', lr4$ge1ge2)

	## output
		lr4 = lr4[order(-lr4$known, -lr4$num_start_site, lr4$intragene, lr4$frame),]
		fusion.list = subset(lr4, select=c('AP7', 'ge1ge2', 'frame','gec_LR','chr_L','pos_L','inEx_L','gene_L','nm_L','exon_L','cdna_L'
					, 'overlap','chr_R','pos_R','inEx_R','gene_R','nm_R','exon_R','cdna_R','intragene','readID'
					, 'fusion_point', 'num_start_site', 'known'))
			write.table(fusion.list, paste(subii, '.fusion.list.post-processing.txt', sep=''), row.names=F, quote=F, sep='\t')

		## consolidate fusion_point
		fusion.tab1 = ddply(fusion.list, .(fusion_point), summarize
						, fusion_point=fusion_point[1], 'AP7'=subii, 'ge1ge2'=ge1ge2[1]
						, 'num_unique_reads'=length(readID), 'frame'=frame[1]
						, 'transcript_L'=nm_L[1], 'transcript_R'=nm_R[1]
						, 'function_L'=inEx_L[1], 'function_R'=inEx_R[1]
						, 'cdna_L'=cdna_L[1], 'cdna_R'=cdna_R[1]
						, 'intragene'=unique(intragene)
						, 'overlap'=overlap[1]
						, 'num_start_site'=num_start_site[1]
						, 'known'=known[1]
						)
		fusion.tab1 = fusion.tab1[order(-fusion.tab1$known, -fusion.tab1$num_start_site, fusion.tab1$intragene, fusion.tab1$frame, -fusion.tab1$num_unique_reads),]

		##==============================
		## filter: overlap length and num of reads
		##==============================
		# fitler those with too few unique reads given long L/R overlap, when overlap >6
		fusion.tab2 = subset(fusion.tab1, ((num_unique_reads > overlap) | overlap <=6))
			head(fusion.tab2)
			out.names = c("AP7", "ge1ge2", "frame", "num_start_site", "num_unique_reads", "fusion_point"
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
		    # inter-gene, NA frame and num_start_site >10 reads  
		    # or known
		ex = subset(fusion.table, known ==1
					| (intragene==0 & frame =='in-frame' & num_start_site >=3)
					| (intragene==0 & num_start_site >=10)
				)[, 1:12]
			head(ex)

		# keep 1st
		ex1 = ex[order(ex$"GeneExon5_GeneExon3", ex$frame, -ex$num_start_site, -ex$num_unique_reads),]
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
				(sample.n = min(10, ngeci))
				(readid = sample(unique(lr4$readID[lr4$ge1ge2 == gei]), sample.n))
				writeLines(readid, 'tmp.readid')
				system('sed -e "s/:umi.*//" -e "s:/1.*::" -e "s:/2.*::" tmp.readid | sort -u > tmp.readid2')

				## get fa
				if (file.exists(paste("../../fq.demux/", subii, ".R2.fq", sep=''))){
					system(paste("grep -f tmp.readid2 -A1 ../../fq.demux/", subii, ".R2.fq | sed 's: :/1 :' | sed 's:^@:>:' | grep -v '\\-\\-' > ", subii, '.', gei, ".txt", sep=''))
				} else if (file.exists(paste("../../fq.demux/", subii, ".R2.fq.gz", sep=''))){
				system(paste("zcat ../../fq.demux/", subii, ".R2.fq.gz | grep -f tmp.readid2 -A1 | sed 's: :/1 :' | sed 's:^@:>:' | grep -v '\\-\\-' > ", subii, '.', gei, ".txt", sep=''))
				} else if (file.exists(paste("../../bam/", subii, ".consolidated.bam", sep=''))){
					system(paste(REPPATH, "/samtools/bin/samtools view ../../bam/", subii
							, ".consolidated.bam | grep -f tmp.readid2 | cut -f1,10 | sed 's/^/>/' | tr '\t' '\n' > ", subii, '.', gei, ".txt", sep=''))
				}
			system('rm tmp.readid*')
			}
	       }
	}
	}
}
