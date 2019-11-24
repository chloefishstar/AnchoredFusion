# breakpoint post annotation processing
# 
# required: 
#	- fu.anno.right
#	- fu.anno.left
#	- fu.breakpointID.stgLeft
#	- fu.breakpointID.stgRight
#
SplitFusion.breakpoint.anno.postscript = function(configFile){

library(data.table)
options(width=204)
options(scipen=999)
library(plyr)
(subii = sub('.*/', '', getwd()))
source(configFile)

##====1: left-right annotations
    colnames = c('readID', 'chrorp_L', 'chr_L', 'orp_L', 'pos_L', 'strand_L'
		, 'num_unique_molecules', 'num_partner_ends', 'num_partner_ends2', 'breakpoint', 'overlap'
		, 'gene_L', 'geneStrand_L', 'inEx_L', 'functiontype_L', 'nm_L', 'exon_L', 'cdna_L'
		, 'chrorp_R', 'chr_R', 'orp_R', 'pos_R', 'strand_R', 'gene_R', 'geneStrand_R', 'inEx_R', 'functiontype_R', 'nm_R', 'exon_R', 'cdna_R')
    lr1 = fread('anno.left.right', sep=' ', header=F, fill=T, stringsAsFactors=F, col.names=colnames)
    lr1$cdna_L = suppressWarnings(as.numeric(lr1$cdna_L))
    lr1$cdna_R = suppressWarnings(as.numeric(lr1$cdna_R))
    lr1$pos_L = suppressWarnings(as.numeric(lr1$pos_L))
    lr1$pos_R = suppressWarnings(as.numeric(lr1$pos_R))
    #head(lr1)
    n.lr1 = nrow(lr1)
    lr1$exonn_L = suppressWarnings(as.numeric(sub('exon','',lr1$exon_L)))
    lr1$exonn_R = suppressWarnings(as.numeric(sub('exon','',lr1$exon_R)))

##====2: For targeted RNA-seq: filter off-targets, i.e. non-ligation end (anchored end) being non-targets

if (panel == 'NA'){
	lr2 = lr1
}else{
	genes = readLines(paste0(path.package("SplitFusion"), '/data/', panel, '.target.genes'))
	if (n.lr1 >0){
	    # extract ligation end
	    lr1$ligEnd0 = sub('.*umi:C', '', lr1$readID)
	    lr1$ligEndChr = sub('P.*', '', lr1$ligEnd0)
	    lr1$ligEndPos = as.numeric(sub('-.*', '', sub('.*P', '', lr1$ligEnd0))) - 100000000

	    lr1$gene_T = ifelse((lr1$chr_L == lr1$ligEndChr) & (abs(lr1$ligEndPos - lr1$pos_L) < 300)
				    , lr1$gene_R
				    , ifelse(lr1$chr_R == lr1$ligEndChr & (abs(lr1$ligEndPos - lr1$pos_R) < 300)
				    , lr1$gene_L
				    , "-"
			    ))
	lr1$exclude = 0
	    if (panel != 'NA'){
		lr1$exclude[!(lr1$gene_T %in% genes)] = 1
	    }
	    #table(lr1$exclude)
	lr2 = subset(lr1, exclude==0)
    }
}

    ## kepp original cdna pos for later frameness calculation
    lr2$cdna_L0 = lr2$cdna_L
    lr2$cdna_R0 = lr2$cdna_R

    ##==== Optional: remove recurrent breakpoints that are not significant (mannually curated)
    #filter.breakpointID = readLines(paste0(path.package("SplitFusion"), '/data/fusion.filter.breakpointID'))
    #lr = subset(lr, !(fusionID %in% filter.breakpointID))

##==== 3: connect left-mid-right
##==== For breakpoints with middle split: to correct breakpoint exon number, cdna, gdna by add/minus middle split size
lr2b = SplitFusion.breakpoint.anno.postscript.mid.anno(configFile = configFile, lr2)

##==== 4: sorting direction
lr3 = SplitFusion.breakpoint.anno.postscript.direction(configFile = configFile, lr2b)
SplitFusion.breakpoint.anno.postscript.direction.sub(configFile = configFile, lr3, sampleID=subii)
}
