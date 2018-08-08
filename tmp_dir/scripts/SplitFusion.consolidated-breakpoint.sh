#!/bin/bash
. ../../run.info.sh

subii=$( pwd | sed "s:.*/::")
# get threads
	head -n 1 ../../iiFreq.txt > _head.1
        gawk '{for (i=1; i<=NF; i++) {if ($i ~ /cpuBWA/) {print $i,i}}}' _head.1 | sed 's/.* /cpuField=/' > _t.sh
        . _t.sh
        cpuBWA=$(grep $subii ../../iiFreq.txt | cut -f $cpuField | sed 's: .*::')

## get reads with SA 
## from Read2 reads only (20180504)
	if [ -f ../../bam/$subii/consolidated.sam ]; then
		grep 'SA:' ../../bam/$subii/consolidated.sam | grep -e '_r2' -e '/2' > _sa.sam
	else
		$REPPATH/samtools/bin/samtools view --threads $cpuBWA ../../bam/$subii.consolidated.bam | grep 'SA:' | grep -e '_r2' -e '/2' > _sa.sam
	fi

$REPPATH/samtools/bin/samtools view --threads $cpuBWA -T /dev/shm/Homo_sapiens_assembly19.fasta -bS _sa.sam > _sa.bam
$REPPATH/bedtools/bin/bedtools bamtobed -cigar -i _sa.bam > _sa.bed0
	## bedtools uses 0-base, change to 1-base:
	# awk '{start = $2+1; $2=start; print}' _sa.bed0 | tr ' ' '\t' > _sa.bed
	# add compatbility to processing bed with '/1' '/2' appended to readID from Bedtools
	awk '{start = $2+1; $2=start; print}' _sa.bed0 | tr ' ' '\t' | sed -e 's:/1::' -e 's:/2::' > _sa.bed
	sort -k4,4b _sa.bed > _sa.bed.s

    ## read length
    #awk '{n=length($10); print $1,n}' _sa.sam > _sa.len
        awk '{n=length($10); print $1,n}' _sa.sam | sed -e 's:/1::' -e 's:/2::' > _sa.len
	sort -k1,1b -k2,2nr _sa.len > _sa.len.s0
	sort -k1,1b -u _sa.len.s0 > _sa.len.s

    join -1 1 -2 4 _sa.len.s _sa.bed.s > _sa.len2

    ## MQ 
    echo | awk -v minMQ=$minMQ '{if ($6 >= minMQ) print}' _sa.len2 | tr ' ' '\t' > _sa.bed.mq

## corrected read length for small size clipping which is either
		## 1)  likely residual adaptor
		## 2)  shorter splice fragments (< minMapLength) that could not be mapped by BWA MEM

		# transform CIGAR strings
	sed 's:\([0-9]\+\)\([SH]\)\([0-9A-Z]\+\)$:\1\t\2\t\3:' _sa.bed.mq > _sa.SMH01
	sed 's:\([0-9]\+\)\([SH]\)$:\t\1\t\2:' _sa.SMH01 > _sa.SMH02

		# two splits (10 fileds only) to _tmpOK
		# complicated ones to _tmp1 for further processing
	awk '{if ($12 ~ /[HS]/) {print > "_tmp1"} else {print $1,$2,$3,$4,$5,$6,$7,$8$9$10 > "_tmpOK"} }' _sa.SMH02 

		# mid split
	#minMapLength=20
	echo | awk -v minMapLength=$minMapLength '{if ($8 >= minMapLength && $11 >= minMapLength) {print > "split.mid"} else {print > "_tmp2"} }' _tmp1

		# remove headi/tail HS and correct effective query read length
    	awk '{if ($8 >= $11) {$2 = $2-$11; newcigar=$8$9$10} else {$2 = $2-$8; newcigar=$10$11$12} print $0,newcigar}' _tmp2 | cut -d ' ' -f 1-7,13 >  _tmpOK2

		# Corrected SMH, with left size and SMH type separated for later processing. 
	cat _tmpOK _tmpOK2 | tr ' ' '\t' | sed 's:\([0-9]\+\)\([SMH]\)\([0-9A-Z]\+\)$:\1\t\2\t\3:' > _sa.SMHc


	##==========================================
	## calculate query start, end
    	awk '{if ($9 == "M") print}' _sa.SMHc | tr ' ' '\t' | cut -f 1-9 > _M
    	awk '{if ($9 != "M") print}' _sa.SMHc | tr ' ' '\t' | cut -f 1-9 > _HS
	
	## print chr.pos 1,2,3,4 in read order, re-calculate read query start and end for '-'
		## add read part (i.e left/right at field $12, breadkpoint chr ($13) and pos ($14)

	    ## if left of query is mapped
		    ## Left: $7+ and M; $7- and HS
		    ## Right: $7- and M; $7+ and HS
	    awk '{if ($7 =="+") {qstart=1; qend=$5-$4+1; print $0,qstart,qend,"left",$3,$5}}' _M | tr ' ' '\t' > _left
	    awk '{if ($7 =="-") {qstart=1; qend=$5-$4+1; mstart=$5; mend=$4; $4=mstart; $5=mend; print $0,qstart,qend,"left",$3,$5}}' _HS | tr ' ' '\t' >> _left

	    ## if right of query is mapped
	    awk '{if ($7 =="+") {qstart=$8+1; qend=$8+$5-$4; print $0,qstart,qend,"right",$3,$4}}' _HS | tr ' ' '\t' > _right
	    awk '{if ($7 =="-") {qstart=$2-$5+$4; qend=$2; mstart=$5; mend=$4; $4=mstart; $5=mend; print $0,qstart,qend,"right",$3,$4}}' _M | tr ' ' '\t' >> _right

	## adjust breakpoint position to group breakpoints fuzzy window (5 bases)
	cat _left _right > _left_right0
	sort -k13,13b -k14,14n _left_right0 > _left_right

			    awk '{if ($13==pre13) {
				    d=$14-pre14; if (d<5){
							    posAdj=preAdj
						    } else {posAdj=$14}	
				    }else{posAdj=$14}
				    ; pre13=$13; pre14=$14; preAdj=posAdj; print $0,posAdj
				    }' _left_right > _left_right.adj
			    awk '{print $0,$14}' _left_right > _left_right.adj

	    ## join left and right to form a query read
		## for pair-end read 2, reverse left <-> right
	    grep 'left'  _left_right.adj > _lefts
	    grep 'right' _left_right.adj > _rights

		sort -k1,1b _lefts > _leftsrt
		sort -k1,1b _rights > _rightsrt

		## $1: readID
		## $2 - $15: left
		## $16 - $29: right
	    join _leftsrt _rightsrt > _left_right.j

	    ##====== breakpoint 
	    ## add at $30, $31, $32
	    awk '{pLeft=$13"_"$15; pRight=$27"_"$29; if (pLeft < pRight){pp=pLeft"-"pRight} else {pp=pRight"-"pLeft}; print $0,pLeft,pRight,pp}' _left_right.j > breakpoint0
		sort -k1,1b breakpoint0 > breakpoint

# DONE breakpoint candidate 
