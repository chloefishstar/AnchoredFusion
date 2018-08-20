#!/bin/bash
. ../../run.info.sh

subii=$( pwd | sed "s:.*/::")
# get threads
	head -n 1 ../../iiFreq.txt > _head.1
        gawk '{for (i=1; i<=NF; i++) {if ($i ~ /cpuBWA/) {print $i,i}}}' _head.1 | sed 's/.* /cpuField=/' > _t.sh
        . _t.sh
        cpuBWA=$(grep $subii ../../iiFreq.txt | cut -f $cpuField | sed 's: .*::')


	## determine Panel field
        head -n 1 ../../iiFreq.txt > _head.1
        gawk '{for (i=1; i<=NF; i++) {if ($i ~ /Panel/) {print $i,i}}}' _head.1 | sed 's/.* /PanelField=/' > _t.sh
        . _t.sh
    	panel=$(grep $subii ../../iiFreq.txt | cut -f $PanelField | sed 's: .*::' | sed 's:[vV][0-9]::')
	
	if [[ $panel == HBV ]]; then
        	refGenome=hbv_hg19.fasta
		minMQ=0
	else
		refGenome=Homo_sapiens_assembly19.fasta
	fi

##==== 1.1. get reads with SA from both Read 1 and 2
	$REPPATH/samtools/bin/samtools view --threads $cpuBWA ../../bam/$subii.consolidated.bam | grep 'SA:' > _sa.sam
	$REPPATH/samtools/bin/samtools view --threads $cpuBWA -T /dev/shm/$refGenome -bS _sa.sam > _sa.bam
	$REPPATH/bedtools/bin/bedtools bamtobed -cigar -i _sa.bam > _sa.bed0
	## bedtools uses 0-base, change to 1-base:
	# awk '{start = $2+1; $2=start; print}' _sa.bed0 | tr ' ' '\t' > _sa.bed
	# add compatbility to processing bed with '/1' '/2' appended to readID from Bedtools
	awk '{start = $2+1; $2=start; print}' _sa.bed0 | tr ' ' '\t' > _sa.bed 

##==== 1.2. get read length
	# By default, number and order of reads in files _sa.sam, _sa.bam and _sa.bed0 are identical. So, paste.
        awk '{n=length($10); print $1,n}' _sa.sam > _sa.len

##==== 2. Join bed and read length, and keep the longest
    	paste _sa.len _sa.bed | tr ' ' '\t' | cut -f2- > _sa.len.bed
		
	cut -f1,5 _sa.len.bed > _len.readID1
	sort -k2,2b -k1,1nr _len.readID1 > _len.readID2
	sort -k2,2b -u _len.readID2 > _len.readID

	sort -k2,2b _len.readID > _len.readID.s
	sort -k4,4b _sa.bed > _sa.bed.s
	join -1 2 -2 4 _len.readID.s _sa.bed.s > _sa.len.bed

##==== 3. remove MQ low
	echo | awk -v minMQ=$minMQ '{OFS="\t"; if ($6 >= minMQ) print $0}' _sa.len.bed > _sa.len.bed.mq

##==== 4. corrected read length for small size clipping which is either
		# 1) likely residual adaptor
		# 2) shorter splice fragments (< minMapLength) that could not be mapped by BWA MEM
	# transform CIGAR strings
	sed 's:\([0-9]\+\)\([SH]\)\([0-9A-Z]\+\)$:\1\t\2\t\3:' _sa.len.bed.mq > _sa.SMH01
	sed 's:\([0-9]\+\)\([SH]\)$:\t\1\t\2:' _sa.SMH01 > _sa.SMH02

##==== 5. output two-split alignment reads (10 fileds only) to _tmpOK
		# more than two-split reads sent to _tmp1 for later processing
	awk '{if ($12 ~ /[HS]/) {print > "_tmp1"} else {print $1,$2,$3,$4,$5,$6,$7,$8$9$10 > "_tmpOK"} }' _sa.SMH02 

		# mid split
		echo | awk -v minMapLength=$minMapLength '{if ($8 >= minMapLength && $11 >= minMapLength) {print > "split.mid"} else {print > "_tmp2"} }' _tmp1

		# remove head/tail HS and correct effective query read length
    		awk '{if ($8 >= $11) {$2 = $2-$11; newcigar=$8$9$10} else {$2 = $2-$8; newcigar=$10$11$12} print $0,newcigar}' _tmp2 | cut -d ' ' -f 1-7,13 >  _tmpOK2

		# Corrected SMH, with left size and SMH type separated for later processing. 
		cat _tmpOK _tmpOK2 | tr ' ' '\t' | sed 's:\([0-9]\+\)\([SMH]\)\([0-9A-Z]\+\)$:\1\t\2\t\3:' > _sa.SMHc

##==== 6. calculate query start, end
	    awk '{if ($9 == "M") print}' _sa.SMHc | tr ' ' '\t' | cut -f 1-9 > _M
	    awk '{if ($9 != "M") print}' _sa.SMHc | tr ' ' '\t' | cut -f 1-9 > _HS

	## 4 positions on a query read (after left and right alignments are merged by ReadID):
		##  left.query.start....left.query.end-[breakpoint]-right.query.start....right.query.end
		##         1                   2      -[breakpoint]-        3                  4
			## 2: left.query.end [corresponding mapping position in the genome is breakpoint 1]
			## 3: right.query.start [corresponding mapping position in the genome is breakpoint 2]

	## re-calculate read query start and end for '-' alignment
	## add read part ('left' or 'right' at field $12), breadkpoint chr ($13) and pos ($14)

	## if left of query is mapped
		    ## Left: $7+ and M; $7- and HS
		    ## Right: $7- and M; $7+ and HS
	    awk '{if ($7 =="+") {qstart=1; qend=$5-$4+1; print $0,qstart,qend,"left",$3,$5}}' _M | tr ' ' '\t' > _left
	    awk '{if ($7 =="-") {qstart=1; qend=$5-$4+1; mstart=$5; mend=$4; $4=mstart; $5=mend; print $0,qstart,qend,"left",$3,$5}}' _HS | tr ' ' '\t' >> _left

	## if right of query is mapped
	    awk '{if ($7 =="+") {qstart=$8+1; qend=$8+$5-$4; print $0,qstart,qend,"right",$3,$4}}' _HS | tr ' ' '\t' > _right
	    awk '{if ($7 =="-") {qstart=$2-$5+$4; qend=$2; mstart=$5; mend=$4; $4=mstart; $5=mend; print $0,qstart,qend,"right",$3,$4}}' _M | tr ' ' '\t' >> _right

##==== 7.1 Aggregate (+/- 5 base window) breakpoints and center on nearest peak
	grep 'left'  _left > _lefts
	grep 'right' _right > _rights

	cut -f 12-14 _lefts > _bk 
	cut -f 12-14 _rights >> _bk
	sort -k1,1b -k2,2b -k3,3n _bk > _bks
	uniq -c _bks | sed -e 's/^ \+//' -e 's/ /\t/g' > _peak1

	awk '{OFS="\t"; chr = $2"_"$3;
		print chr,$0
	}' _peak1 > _peak2

	awk '{OFS="\t";
		if ($1 == pre1){
			if ($5 - pre5 <=5){
				if (preCnt > $2){
					pos = prePos;
					cnt = preCnt
				} else {
					pos = $5; cnt = $2; prePos=$5; preCnt=$2
				};
			} else {
				pos = $5; cnt = $2; prePos=$5; preCnt=$2
			}
		} else {
			cnt = $2; pos = $5; prePos=$5; preCnt=0
		}
		pre1=$1; pre5=$5;
		print $0,pos,cnt
	}' _peak2 > _peak2a

	tac _peak2a > _peak2at
	awk '{OFS="\t";
		if ($1 == pre1){
			if (pre6 - $6 <=5){
				if (preCnt > $7){
					pos = prePos;
					cnt = preCnt
				} else {
					pos = $6; cnt = $7; prePos=$6; preCnt=$7
				};
			} else {
				pos = $6; cnt = $7; prePos=$6; preCnt=$7
			}
		} else {
			pos = $6; cnt = $7; prePos=$6; preCnt=$7
		}
		pre1=$1; pre6=$6;
		print $0,pos,cnt
	}' _peak2at > _peak2b
	tac _peak2b > loacl_breakpoint_peak

	# merge back
	awk '{OFS="\t"; print $1"_"$5,$4"_"$8}' loacl_breakpoint_peak > _loacl_breakpoint_peak
	sort -k1,1b _loacl_breakpoint_peak > _loacl_breakpoint_peak.s

	sed -e 's/\t/_/12' -e 's/\t/_/12' _lefts > _lefts.p
	sort -k12,12b _lefts.p > _lefts.ps
	join -1 12 -2 1 _lefts.ps _loacl_breakpoint_peak.s | tr ' ' '\t' | cut -f 2- > _lefts.psj
	
	sed -e 's/\t/_/12' -e 's/\t/_/12' _rights > _rights.p
	sort -k12,12b _rights.p > _rights.ps
	join -1 12 -2 1 _rights.ps _loacl_breakpoint_peak.s | tr ' ' '\t' | cut -f 2- > _rights.psj

##==== 7.2 join left and right by ReadID to form a query read
	# sort by Query position too
	sort -k1,1b -k10,10n _lefts.psj > _leftsrt
	sort -k1,1b -k10,10n _rights.psj > _rightsrt

	    ## left: $2 - $15; right: $16 - $29
	join _leftsrt _rightsrt > _left_right.j
	
##==== 8. update mapping start.site.umi
	sed 's/:umi:[^-]\+/ /' _left_right.j | awk '{
		pos2 = $5 + 100000000;
		ssu = ":umi:C"$4"P"pos2;
		print $1""ssu""$2,$0
	}' | cut -d ' ' -f 1,4- > _left_right.ssu

##==== 9. breakpoint candidates
	## left: 2-12; right 13-23
	awk '{if ($12 < $23){
			pp=$12"__"$23
		} else {pp=$23"__"$12
		}; 
		print $0,pp
	}' _left_right.ssu > _breakpoint.noFilter
	sort -k1,1b _breakpoint.noFilter > breakpoint.noFilter

# DONE breakpoint candidate 
