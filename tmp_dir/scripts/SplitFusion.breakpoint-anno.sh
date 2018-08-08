#!/bin/bash

. ../../run.info.sh

##==== 3: Annotate breakpoint gene, exon, cDNA position

## Adjust breakpoint for annotation
#	- introduce overlap-removed point (orp)
#	- to correctly annotate exon and cDNA position (some alignment ends in intron could be overlap with,
#		 and should belong to, next (sense if on cDNA) split alignment)
#	- Add 6 bases to breakpoint to make orp, which needs to be accounted for in later frame determination
    awk '{	if ($10=="+") {
		    orpL = $8-$35-6
	    } else if ($10=="-") {
		    orpL = $8+$35+6
	    };

	    if ($24=="+") {
		    orpR = $21+$35+6
	    } else if ($24=="-") {
		    orpR = $21-$35-6
	    };

	    print $6"_"orpL,$6,orpL,$8,$10,$20"_"orpR,$20,orpR,$21,$24,$35,$1,$2,$3,$4
    }' fusion.candidates > _orp

## Make mock variants for annotation to get gene, exon and fucntion info
## Left fileds 1-5, Right fileds 6-10, readID 15
    cut -d ' ' -f 1-5,11- _orp > __orpLeft 
    cut -d ' ' -f 6-10,15 _orp > __orpRight
    
    cat __orpLeft __orpRight | cut -d ' ' -f 2,3 | sort -u > __breakpoint.for.anno0
    awk '{print $1,$2,$2,"A","A"}' __breakpoint.for.anno0 > __breakpoint.for.anno

## Annotate
    perl $REPPATH/annovar/table_annovar.pl __breakpoint.for.anno $REPPATH/annovar/humandb/ -buildver hg19 -out __breakpoint.annotated -remove -protocol refGene -operation g -nastring NA
    Rscript $pipelinePATH/scripts/chr.pos.anno.extraction.R __breakpoint.annotated ## generate .ext0
    sort -k1,1b __breakpoint.annotated.ext0 > __breakpoint.annotated.extr

## Merge annotation back to read
	sort -k1,1b __orpLeft > __orpLeft.s
	sort -k1,1b __orpRight > __orpRight.s
	join __orpLeft.s  __breakpoint.annotated.extr > __anno.left
	join __orpRight.s __breakpoint.annotated.extr > __anno.right

	sort -k10,10b __anno.left > _anno.left
	sort -k6,6b __anno.right > _anno.right
	join -1 10 -2 6 _anno.left _anno.right > anno.left.right

    ## anno middle split
    awk '{mid = $4 + 10; print $3,mid,mid,"A","A",$1,$4,$5,$6}' split.mid > _mid.for.anno0
    tr ' ' '\t' < _mid.for.anno0 | cut -f1-5 | sort -u > mid.for.anno

    perl $REPPATH/annovar/table_annovar.pl mid.for.anno $REPPATH/annovar/humandb/ -buildver hg19 -out mid.anno -remove -protocol refGene -operation g -nastring NA
    Rscript $pipelinePATH/scripts/chr.pos.anno.extraction.R mid.anno ## generate .ext0

    sed 's: :_:' _mid.for.anno0 | sort -k1,1b > _mid.for.anno1
    sort -k1,1b mid.anno.ext0 > _mid.anno.ext
    join _mid.for.anno1 _mid.anno.ext | cut -d ' ' -f5- > mid.anno2

rm __*

