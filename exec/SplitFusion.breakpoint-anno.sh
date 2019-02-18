#!/bin/bash

. $1

##==== 3: Annotate breakpoint gene, exon, cDNA position

## Adjust breakpoint for annotation
#	- introduce overlap-removed point (orp)
#	- to correctly annotate exon and cDNA position (some alignment ends in intron could be overlapped with,
#		 and should belong to, next (on the sense strand if on cDNA) split alignment)
#	- Add 6 bases to breakpoint to make orp, which needs to be accounted for in later frame determination
    awk '{OFS="\t"; if ($10=="+") {
		    orpL = $8-$26-6
	    } else if ($10=="-") {
		    orpL = $8+$26+6
	    };
	    if ($20=="+") {
		    orpR = $17+$26+6
	    } else if ($20=="-") {
		    orpR = $17-$26-6
	    };
	    print $6"_"orpL,$6,orpL,$8,$10,$16"_"orpR,$16,orpR,$17,$20,$4,$1,$2,$3,$25,$26
    }' breakpoint.candidates > _orp

## Make mock variants for annotation to get gene, exon and fucntion info
## Left fileds 1-5, Right fileds 6-10, readID 11
    cut -f 1-5,11- _orp > __orpLeft 
    cut -f 6-11 _orp > __orpRight
    
    cat __orpLeft __orpRight | cut -f 2,3 | sort -u > __breakpoint.for.anno0
#    awk '{print $1,$2,$2,"A","A"}' __breakpoint.for.anno0 > __breakpoint.for.anno
    awk '{OFS="\t"; print $1,$2,".",".","A",".",".","."}' __breakpoint.for.anno0 > __breakpoint.for.anno

## Annotate
    #perl $REPPATH/annovar/table_annovar.pl __breakpoint.for.anno $REPPATH/annovar/humandb/ -buildver hg19 -out __breakpoint.annotated -remove -protocol refGene -operation g -nastring NA
    #Rscript $pipelinePATH/scripts/chr.pos.anno.extraction.R __breakpoint.annotated ## generate .ext0
    $java -jar $snpEff -canon $snpEff_ref __breakpoint.for.anno > __breakpoint.annotated  ### Note: annotation result varys with different version of $snpEff_ref
#    $R -f $SplitFusionPath/R/exon.cds.extraction.R --args __breakpoint.annotated  ###R version###
    $R -e 'library(SplitFusion);exon.cds.extraction(input = "__breakpoint.annotated")'

    sort -k1,1b __breakpoint.annotated.ext0 > __breakpoint.annotated.extr

## Merge annotation back to read
	sort -k1,1b __orpLeft > __orpLeft.s
	sort -k1,1b __orpRight > __orpRight.s
	join __orpLeft.s  __breakpoint.annotated.extr > __anno.left
	join __orpRight.s __breakpoint.annotated.extr > __anno.right

	sort -k6,6b __anno.left > _anno.left
	sort -k6,6b __anno.right > _anno.right
	join -1 6 -2 6 _anno.left _anno.right > anno.left.right

    ## anno middle split
    awk '{mid = $5 + 10; print $4,mid,mid,"A","A",$1,$4,$5,$6}' split.mid > _mid.for.anno0
    tr ' ' '\t' < _mid.for.anno0 | cut -f1-5 | sort -u |awk '{OFS="\t"; print $1,$2,".",".","A",".",".","."}' > mid.for.anno

#    perl $REPPATH/annovar/table_annovar.pl mid.for.anno $REPPATH/annovar/humandb/ -buildver hg19 -out mid.anno -remove -protocol refGene -operation g -nastring NA
#    Rscript $pipelinePATH/scripts/chr.pos.anno.extraction.R mid.anno ## generate .ext0
    $java -jar $snpEff -canon $snpEff_ref mid.for.anno > mid.anno  ### Note: annotation result varys with different version of $snpEff_ref
#    $R -f $SplitFusionPath/R/exon.cds.extraction.R --args mid.anno  ###R version###
    $R -e 'library(SplitFusion);exon.cds.extraction(input = "mid.anno")'

    sed 's: :_:' _mid.for.anno0 | sort -k1,1b > _mid.for.anno1
    sort -k1,1b mid.anno.ext0 > _mid.anno.ext
    join _mid.for.anno1 _mid.anno.ext | cut -d ' ' -f2,5- > mid.anno2

#rm __*
