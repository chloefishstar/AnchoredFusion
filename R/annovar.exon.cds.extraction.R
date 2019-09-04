
## Keep one of ANNOVAR annotations (NM, exon#, cdna#) by the desried transcript -- here use the refGene transript with most exons

# setwd('"$analysisRunFull"/bam/A23-P706')
# in.anno = 'mid.anno'

annovar.exon.cds.extraction <- function(input){

options(width=204)
#source(runInfo)
#.libPaths(DEPATH)

in.anno = input
#in.anno = 'fu.anno'
refNM = read.table(paste(path.package("AnchoredFusion"), '/data/refGene.mostExon.NM', sep=''), sep='\t', header=T, stringsAsFactors=F)
names(refNM) = c('Gene.refGene', 'refNM')
#head(refNM)

   d1 = read.table(paste(in.anno, '.hg19_multianno.txt', sep=''), sep='\t', header=T, stringsAsFactors=F)
	d1$Gene.refGene = sub(';.*', '', d1$Gene.refGene)
	d1$Gene.refGene = sub(',.*', '', d1$Gene.refGene)
	d2 = merge(d1, refNM, by='Gene.refGene', all.x=T)	
	d2$tmp1 = mapply(sub, pattern=paste('.*', d2$refNM, ':', sep=''), replacement='', x=d2$AAChange.refGene)
	
	    # when ANNOVAR.refGenes do not match refNM (most exons refGene)
	    d2$s1.4 = substr(d2$tmp1, 1, 4)
	    d2$notMatch = ifelse(d2$Func.refGene == 'exonic' & d2$s1.4 != 'exon', 1, 0) 
		    table(d2$notMatch)
	    d2$tmp1[d2$notMatch==1] = sub('.*NM_', 'NM_', d2$AAChange.refGene[d2$notMatch==1])
		    head(d2[d2$notMatch==1,])
	    d2$refNM[d2$notMatch==1] = sub(':.*', '', d2$tmp1[d2$notMatch==1])
	    d2$tmp1[d2$notMatch==1] = sub('.*:exon', 'exon', d2$tmp1[d2$notMatch==1])
		    head(d2[d2$notMatch==1,])

	d2$tmp2 = sub(',.*', '', d2$tmp1)
	d2$exon = sub(':.*', '', d2$tmp2)
	d2$tmp3 = sub('.*c.[ACTG]', '', d2$tmp2)
	d2$cdna = sub('[ACTG].*', '', d2$tmp3)
	d2$tmp4 = sub(':p.*', '', d2$tmp3) 
	d2$tmp5 = sub("\\d+", '', d2$tmp4) 
	d2$geneStrand = ifelse(d2$tmp5 == 'A', '+', '-')
	d2$chrpos = paste(d2$Chr, d2$Start, sep='_')
	d2 = d2[order(d2$chrpos),]
	d2$ExonicFunc.refGene = gsub(' ', '_', d2$ExonicFunc.refGene)

		table(d2$Func, d2$exon, useNA='always')
	d2$exon[is.na(d2$exon)] = d2$Func.refGene[is.na(d2$exon)]

	d3 = subset(d2, select=c('chrpos', 'Gene.refGene', 'geneStrand', 'Func.refGene','ExonicFunc.refGene', 'refNM', 'exon', 'cdna'))
		head(d3)

write.table(d3, paste(in.anno, '.ext0', sep=''), col.names=F, row.names=F, quote=F, sep='\t')

}
