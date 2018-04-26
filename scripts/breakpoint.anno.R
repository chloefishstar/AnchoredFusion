options(width=204)
source('../../run.info')
library(plyr)
subii = sub('.*/', '', getwd())

## target genes
ii = read.table('../../iiFreq.txt', sep='\t', stringsAsFactors=F, header=T)
head(ii)
(panel = sub("\\..*", '', ii$Panel[ii$AP7 == subii][1]))
(genes = readLines(paste(SplitFusionPath, '/scripts/', panel, '.target.genes', sep='')))
## left
    colnames = c('chrorp', 'chr', 'pos', 'strand', 'overlap', 'readID', 'gene', 'geneStrand', 'inEx', 'functiontype', 'nm', 'exon', 'cdna')
    left = read.table('fu.anno.left', sep=' ', header=F, fill=T, stringsAsFactors=F, col.names=colnames)
    left$p = paste(left$chr, left$pos, sep='_')
    left$cdna = as.numeric(left$cdna)
    head(left)
    nrow(left)
    #(funtab = table(left$functiontype))
    left$exonn = suppressWarnings(as.numeric(sub('exon','',left$exon)))
    #head(left[left$gene=='EML4',])
    #nrow(left[left$gene=='CAPN12',])

## right
    right = read.table('fu.anno.right', sep=' ', header=F, fill=T, stringsAsFactors=F, col.names=colnames)
    right$p = paste(right$chr, right$pos, sep='_')
    right$cdna = as.numeric(right$cdna)
    head(right)
    nrow(right)
    #(funtab = table(right$functiontype))
    right$exonn = suppressWarnings(as.numeric(sub('exon', '', right$exon)))
    #head(right[right$gene=='CAPN12',])

## left-right
lr0 = merge(left, right, by='readID', suffixes = c("_L","_R"))
head(lr0, 3)

## Among one-sided stagger fusion, remove fusion candidates with cliff side not on target genes
system("join -a1 -a2 fu.breakpointID.stgLeft fu.breakpointID.stgRight | grep -ve 'stgLeft stgRight' > stgSingle")
(n.single = as.numeric(system("wc -l stgSingle | sed 's: .*::' ", intern=T)))
if (n.single >0){
	stgSingle = read.table('stgSingle', header=F, stringsAsFactors=F, sep=' ')
    head(stgSingle)
    names(stgSingle) = c('fusionID', 'stgSide')
    lr0$fusionID = ifelse(lr0$p_L < lr0$p_R
			    , paste(lr0$p_L, lr0$p_R, sep='-')
			    , paste(lr0$p_R, lr0$p_L, sep='-')
			    )
    	#head(lr0, 3)
    lr00 = merge(lr0, stgSingle, by ='fusionID', all.x=T)
    	#head(lr00, 3)
	#tail(lr00, 3)
    lr00$read12 = sub('.*/', '', lr00$readID)
	#table(lr00$read12)
    lr00$exclude = 0
    lr00$exclude[     (lr00$stgSide =='stgLeft' & !(lr00$gene_R %in% genes))
		    | (lr00$stgSide =='stgRight' & !(lr00$gene_L %in% genes))
		    ] = 1
	table(lr00$exclude)
    lr000 = subset(lr00, exclude==0)
}else{
    lr000 = lr0
}

## exlude off targets (neither of the partners is on target)
lr = subset(lr000, (lr000$gene_L %in% genes) | (lr000$gene_R %in% genes))


##== remove recurrent breakpoints that are not significant (mannually curated)
filter.breakpointID = readLines(paste(SplitFusionPath, '/scripts/fusion.filter.breakpointID', sep=''))
lr = subset(lr, !(fusionID %in% filter.breakpointID))

#head(subset(lr00, exclude==1))
#chk='CD74'
#head(subset(lr, gene_L ==chk), 6)
#head(subset(lr, gene_R ==chk), 6)

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
#head(lr2)
#table(lr2$exonn_R[lr2$gene_R == 'EML4'])
#head(subset(lr2, exonn_R ==5  & gene_R == 'EML4'), 3)
#head(subset(lr2, exonn_R ==6  & gene_R == 'EML4'), 3)

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
nrow(lr2)
n.lr3 = nrow(lr3)

if (n.lr3 >0){
#head(lr3)
#nrow(subset(lr3, gene_R=='CAPN12'))

## add the 6 base back to cDNA position
lr3$cdna_L = lr3$cdna_L + 6
lr3$cdna_R = lr3$cdna_R - 6

#head(lr3[lr3$gene_L=='EML4',], 6)
#tail(lr3[lr3$gene_R=='CAPN12',], 6)

## intra-gene
lr3$intragene = ifelse(lr3$gene_L == lr3$gene_R, 1, 0)
#head(lr3,3)
#str(lr3)
lr3$exon = substr(lr3$exon_L,1,4)
#tail(lr3[lr3$gene_L=='ERBB2',])
#head(lr3[lr3$exon=='ERBB',])
#table(lr3$inEx_L, lr3$exon_L, useNA='always')

lr3$exonD = 0
lr3$exonD = ifelse(lr3$gene_L == lr3$gene_R, lr3$exonn_R - lr3$exonn_L, lr3$exonD)
lr3$gec_L = paste(lr3$gene_L, ' ', lr3$exon_L, ' c.', lr3$cdna_L, ' (', lr3$nm_L, ')', sep='')
lr3$gec_R = paste(lr3$gene_R, ' ', lr3$exon_R, ' c.', lr3$cdna_R, ' (', lr3$nm_R, ')', sep='')
lr3$gec_LR = paste(lr3$gec_L, lr3$gec_R, sep='---')
lr3$ge1ge2 = paste(lr3$gene_L, '_', lr3$exon_L,'---', lr3$gene_R, '_', lr3$exon_R, sep='')
#head(lr3,3)


##======================
## frame status
##======================
    lr3$frameD = (lr3$cdna_R0 %% 3 - (lr3$cdna_L0 + abs(lr3$overlap_L)) %% 3) %% 3
    lr3$frame = ifelse(lr3$frameD ==1, 'in-frame', 'out-frame')
	#head(lr3)

	## expect mostly in-frame or _NA_. Otherwise, might need inspection
	#table(lr3$frame,useNA='always')
    lr3$frame[is.na(lr3$frame)] = '_NA_'
	#table(lr3$frame,useNA='always')
			     
    lr3$neighb = ifelse(lr3$exonD==1, 1, 0)
    lr3$neighb[is.na(lr3$neighb)] = 0
## order
	lr3 = lr3[order(lr3$intragene, lr3$frame, lr3$gec_LR),]
	#head(lr3)

	## chk
	# head(subset(lr3, gene_L=='EML4' & gene_R=='ALK'),6)
	# tail(subset(lr3, gene_L=='EML4' & gene_R=='ALK'),6)
	# table(subset(lr3, gene_L=='EML4' & gene_R=='ALK')$frame)

## remove:
#	neighbour exons
#	out-frame
#	same exon
#chk=lr3
#chk.g = 'ROS1'
#nrow(chk[grep(chk.g, chk$gene_R),])
	#lr40 = subset(lr3, !(frame == 'out-frame' | neighb==1 | (intragene==1 & exon_L==exon_R)))
	## keep out-frame of different genes. FGFRs fusion can be functional even out-frame
	lr40 = subset(lr3, !(neighb==1 | (intragene==1 & exon_L==exon_R)
				| (intragene==1 & frame=='out-frame')
				))
	# tail(subset(lr3, (intragene==1 & exon_L==exon_R)), 6)
	lr4 = lr40[!duplicated(lr40$readID),]	
	#head(lr4,4)
	nrow(lr3); nrow(lr4)
	lr4$AP7 = subii
	lr4$gec_LR = gsub('\\(', '.', lr4$gec_LR)
	lr4$gec_LR = gsub('\\)', '.', lr4$gec_LR)
	lr4$ge1ge2 = gsub('\\(', '.', lr4$ge1ge2)
	lr4$ge1ge2 = gsub('\\)', '.', lr4$ge1ge2)

## output
	# fusion list
	fusion.list = subset(lr4, select=c('AP7', 'ge1ge2', 'frame','gec_LR','chr_L','pos_L','inEx_L','gene_L','nm_L','exon_L','cdna_L'
				,'chr_R','pos_R','inEx_R','gene_R','nm_R','exon_R','cdna_R','intragene','readID'))
	head(fusion.list)
	tail(fusion.list)
	    #	head(subset(fusion.list, gene_L=='KIF5B' & gene_R=='RET'))
	    #	head(subset(fusion.list, gene_L=='FGFR2' & gene_R=='CAPN12'))
	    #	nrow(fusion.list[grep('CAPN12', fusion.list$gene_R),])

	name1 = sub('gec_LR', "Gene_Exon_cDNA_5'_3'", names(fusion.list))
	name2 = sub('_L', "_5'", name1)
	name3 = sub('_R', "_3'", name2)
	name4 = sub('ge1ge2', "GeneExon5'---GeneExon3'", name3)
	fusion.list.xls = fusion.list
	names(fusion.list.xls) = name4
    	head(fusion.list.xls)
    	write.table(fusion.list.xls, paste(subii, '.fusion.list.xls', sep=''), row.names=F, quote=F, sep='\t')

	# fusion tab
	fusion.tab1 = ddply(fusion.list, .(gec_LR, frame), summarize, 'AP7'=subii, 'ge1ge2'=unique(ge1ge2)
					   , 'num_unique_reads'=length(gec_LR), 'frame'=frame[1]
					    , 'transcript_L'=nm_L[1], 'transcript_R'=nm_R[1]
					    , 'function_L'=inEx_L[1], 'function_R'=inEx_R[1]
					    , 'cdna_L'=cdna_L[1], 'cdna_R'=cdna_R[1]
					    , 'intragene'=unique(intragene)
	)
	fusion.tab1 = fusion.tab1[order(fusion.tab1$intragene, fusion.tab1$frame, -fusion.tab1$num_unique_reads),]
	fusion.tab1[1:6, 3:12]
	    ## keep one per fusionExon
	    fusion.tab2 = subset(fusion.tab1, num_unique_reads >1) ## more than 1 unique read
	    head(fusion.tab2)
	name1 = sub('gec_LR', "Gene_Exon_cDNA_5'_3'", names(fusion.tab2))
	name2 = sub('_L', "_5'", name1)
	name3 = sub('_R', "_3'", name2)
	name4 = sub('ge1ge2', "GeneExon5'---GeneExon3'", name3)
	fusion.tab.xls = fusion.tab2
	names(fusion.tab.xls) = name4
    	head(fusion.tab.xls,3)
	fusion.tab.xls2 = fusion.tab.xls[,c(2:5,1,6:ncol(fusion.tab.xls))]
    	head(fusion.tab.xls2)
	
	# order
	fusion.tab.xls2$order[fusion.tab.xls2$frame=="in-frame"] = 1
	fusion.tab.xls2 = fusion.tab.xls2[order(fusion.tab.xls2$order, fusion.tab.xls2$intragene, -fusion.tab.xls2$num_unique_reads),]
    	head(fusion.tab.xls2)
	
	write.table(fusion.tab.xls2, paste(subii, '.fusion.table.xls', sep=''), row.names=F, quote=F, sep='\t')


## extract fusion reads fasta, only for
	# inter-gene, in-frame
	# inter-gene, NA frame and >100 reads  
	# or forced 
	(forced = readLines(paste(SplitFusionPath, '/scripts/fusion.forced.txt', sep='')))
	# table(fusion.tab2$ge1ge2)	
	# head(fusion.tab2)
	ex = subset(fusion.tab2, (intragene==0 & frame=='in-frame')
				| (intragene==0 & frame =='_NA_' & num_unique_reads >=100)
				| ge1ge2 %in% forced
				| (intragene==0 & (sub('---.*', '---', ge1ge2) %in% forced))
				| (intragene==0 & (sub('.*---', '---', ge1ge2) %in% forced))
			)	
	# order
	ex$order[ex$ge1ge2 %in% forced] = 1
	ex$order[ex$frame=="in-frame"] = 2
	ex = ex[order(ex$order, ex$intragene, -ex$num_unique_reads),]
    		head(ex)
	ex = ex[order(ex$order, ex$num_unique_reads),]
		ex.xls = ex[,c(2:5,1,6:ncol(ex))]
		head(ex.xls)
		write.table(ex.xls, paste(subii, '.brief.summary', sep=''), row.names=F, quote=F, sep='\t')

	##=======================================
        ## output FASTA reads for max 10 fusion candidates
	##=======================================
	(nex = min(10, nrow(ex)))
            if (nex >0){
                ## show max 25 example reads per candidate
		i=1
                for (i in 1:nex){
                        ## get read ID
			(gei = ex$ge1ge2[i])
			(geci = ex$gec_LR[i])
			(ngeci = ex$num_unique_reads[i])
			(sample.n = min(25, ngeci))
			(readid = sample(lr4$readID[lr4$gec_LR == geci], sample.n))
			writeLines(readid, 'tmp.readid')
			system('sed "s/::umi.*//" tmp.readid | sort -u > tmp.readid2')

                        ## get fa
                        system(paste("grep -f tmp.readid2 -A1 _sa.sam | cut -f1,2,10 | sed -e 's: :/1 :' -e 's:^:>:' | sed 's/\t/_flag_/' | tr '\t' '\n' | grep -v '\\-\\-' >> ", subii, '.', gei, ".fa", sep=''))
                system('rm tmp.readid*')
                }
       }
}
