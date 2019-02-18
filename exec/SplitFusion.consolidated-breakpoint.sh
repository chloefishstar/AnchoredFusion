
#!/bin/bash
. $1

subii=$( pwd | sed "s:.*/::")
# get threads
	head -n 1 $sampleInfo > _head.1
        gawk '{for (i=1; i<=NF; i++) {if ($i ~ /cpuBWA/) {print $i,i}}}' _head.1 | sed 's/.* /cpuField=/' > _t.sh
        . _t.sh
        cpuBWA=$(grep $subii $sampleInfo | cut -f $cpuField | sed 's: .*::')


	## determine Panel field
        head -n 1 $sampleInfo > _head.1
        gawk '{for (i=1; i<=NF; i++) {if ($i ~ /Panel/) {print $i,i}}}' _head.1 | sed 's/.* /PanelField=/' > _t.sh
        . _t.sh
    	panel=$(grep $subii $sampleInfo | cut -f $PanelField | sed 's: .*::' | sed 's:[vV][0-9]::')
	

##==== 1.1. get reads with SA from both Read 1 and 2

#	$samtools view -@ $cpuBWA $bam_path/$subii.consolidated.bam | grep 'SA:' > _sa.sam
#	$samtools view -@ $cpuBWA -T $hgRef -bS _sa.sam > _sa.bam

#	$bedtools bamtobed -cigar -i _sa.bam > _sa.bed0
	$bedtools bamtobed -cigar -i $bam_path/$subii.consolidated.bam > _sa.bed0

	## bedtools uses 0-base, change to 1-base:
	# awk '{start = $2+1; $2=start; print}' _sa.bed0 | tr ' ' '\t' > _sa.bed
	# add compatbility to processing bed with '/1' '/2' appended to readID from Bedtools
	awk '{start = $2+1; $2=start; print}' _sa.bed0 | tr ' ' '\t' > _sa.bed 

	# re-format read ID if not already in the umi:C.P format
	goodformat=$(head -n 1 _sa.bed | cut -f 4 | grep umi: | sed 's:.*umi:umi:' | grep C | grep P | wc -l)
	if [ $goodformat -eq 0 ]; then
		mv _sa.bed _sa.bed1
		sort -k4,4 _sa.bed1 > _sa.bed2
		sed 's/::umi/:umi/' _sa.bed2 | sed 's/:umi:/\tumi\t/' |\
		gawk '{OFS="\t"; if ($4 != preID){
					if ($8 == "+"){posC = 100000001 + $2} else {posC = 100000001 + $3};
                                	umi="C"$1"P"posC"-"$6
                        	} else {umi=preUmi};
                        preUmi=umi;
                        preID=$4;
			$6=umi;
			print $0 > "_sa.bed3"
                }' 2>/dev/null
		cat _sa.bed3
		sed 's/\tumi\t/:umi:/' _sa.bed3 > _sa.bed
	fi

##==== 1.2. get read length
	# By default, number and order of reads in files _sa.sam, _sa.bam and _sa.bed0 are identical. So, paste.
#        awk '{n=length($10); print $1,n}' _sa.sam > _sa.len
	$samtools view -@ $cpuBWA $bam_path/$subii.consolidated.bam | awk '{n=length($10); print $1,n}' > _sa.len

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

##==== 6. calculate query start, end
	# corrected read length for small size clipping which is either
	# transform CIGAR strings
	cut -f 2,7,8 -d ' ' _sa.len.bed.mq > _sa.SMH1
	sed -e 's/M/ M /g' -e 's/S/ S /g' -e 's/H/ H /g' _sa.SMH1 > _sa.SMH2
	awk '{if ($4=="M"){
			if ($2=="+"){
				start=1; end=$3
			}else if ($2=="-"){
				start=$1-$3+1; end=$1
			}
		} else if ($6=="M"){
			if ($2=="+"){
				start=$3+1; end=$3+$5
			}else if ($2=="-"){
				start=$1-$3-$5+1; end=$1-$3
			}
		} else {start=0; end=0}; print start,end}' _sa.SMH2 > _sa.SMH3
	paste _sa.len.bed.mq _sa.SMH3 | tr ' ' '\t' > _sa.SMH4
	
	sort -k1,1b -k9,9n _sa.SMH4 > _sa.SMH4s

##==== sep left/right

	awk '{if ($1==pre1){
			n=n+1
		} else {n=1};
		if (n==1) {
			print $0 > "_left";
			if (pren !=1) {print pre0 > "_right"}
		} else if (n>2){
			print pre0 > "_split.mid"
		};
		pre1=$1; pren=n; pre0=$0;
		print n,$0 > "_sa.SMH4sn"
	}' _sa.SMH4s


#=====================================================================
#=====================================================================


	## 4 positions on a query read (after left and right alignments are merged by ReadID):
		##  left.query.start....left.query.end-[breakpoint]-right.query.start....right.query.end
		##         1                   2      -[breakpoint]-        3                  4
			## 2: left.query.end [corresponding mapping position in the genome is breakpoint 1]
			## 3: right.query.start [corresponding mapping position in the genome is breakpoint 2]

	## add read part ('left' or 'right'), breadkpoint (chr and pos)

	    awk '{if ($7 =="+") {
				print $0,"left",$3,$5
			} else {
				mstart=$5; mend=$4; $4=mstart; $5=mend; 
				print $0,"left",$3,$5
			}
		}' _left | tr ' ' '\t' > _lefts

	    awk '{if ($7 =="+") {
				print $0,"right",$3,$4
			} else {
	    			mstart=$5; mend=$4; $4=mstart; $5=mend; 
				print $0,"right",$3,$4
			}
		}' _right | tr ' ' '\t' | sed 1d > _rights

##==== 7.1 Aggregate (+/- 5 base window) breakpoints and center on nearest peak
	cut -f 11-13 _lefts > _bk 
	cut -f 11-13 _rights >> _bk
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
	tac _peak2b > local_breakpoint_peak

	# merge back
	awk '{OFS="\t"; print $1"_"$5,$4"_"$8}' local_breakpoint_peak > _local_breakpoint_peak
	sort -k1,1b _local_breakpoint_peak > _local_breakpoint_peak.s

	sed -e 's/\t/_/11' -e 's/\t/_/11' _lefts > _lefts.pk
	sed -e 's/\t/_/11' -e 's/\t/_/11' _rights > _rights.pk
		sort -k11,11b _lefts.pk > _lefts.pks
		sort -k11,11b _rights.pk > _rights.pks

	join -1 11 -2 1 _lefts.pks _local_breakpoint_peak.s > _lefts.pk1	
	join -1 11 -2 1 _rights.pks _local_breakpoint_peak.s > _rights.pk1	

##=== join leftmost and rightmost
	cut -f2- -d ' ' _lefts.pk1 | sort -k1,1b > _lefts.pk2
	cut -f2- -d ' ' _rights.pk1 | sort -k1,1b > _rights.pk2
	join _lefts.pk2 _rights.pk2 > _left_right.j

##==== 9. breakpoint candidates
	## left: 2-11; right 12-21
	awk '{if ($11 < $21){
			pp=$11"__"$21
		} else {pp=$21"__"$11
		}; 
		print $0,pp
	}' _left_right.j > _breakpoint.noFilter1
	sort -k1,1b _breakpoint.noFilter1 > _breakpoint.noFilter2

##====	correct breakpoint for those contain mid
	cut -f1 _split.mid | sort -u > _mid.id
	sort -k2,2b _sa.SMH4sn > breakpoint.reads # _sa.SMH4sns
	join -1 1 -2 2 _mid.id breakpoint.reads > split.mid

	awk '{if ($1==pre1){
		diff = $5 - pre5;
		if ($4 != pre4 || diff > 100000 || diff < -100000){
			if (pre8=="+"){bkp1=pre4"_"pre6} else {bkp1=pre4"_"pre5};
			if ($8=="+"){bkp2=$4"_"$5} else {bkp2=$4"_"$6}
                   };
        	};
	if (bkp1 < bkp2){bkp = bkp1"__"bkp2};
	if (bkp1 > bkp2){bkp = bkp2"__"bkp1};
        if (bkp != ""){print $1,bkp > "_sa.mid.bkp"}
        diff=""; bkp=""; bkp1=""; bkp2=""
        pre1=$1; pre4=$4; pre5=$5; pre6=$6; pre8=$8;
	}' split.mid

	sort -k1,1b _sa.mid.bkp > _sa.mid.bkps
	join -a1 _breakpoint.noFilter2 _sa.mid.bkps > _breakpoint.noFilter3

	awk '{if (NF==23) {$22=$23; $23=""}; print}' _breakpoint.noFilter3 > _breakpoint.noFilter.bkp.corrected

	join _breakpoint.noFilter.bkp.corrected _mid.id > breakpoint.noFilter.w.mid
	join -v 1 _breakpoint.noFilter.bkp.corrected _mid.id > breakpoint.noFilter.wo.mid

# DONE breakpoint candidate
