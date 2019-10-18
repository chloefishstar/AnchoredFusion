#!/bin/bash
. $1

subii=$( pwd | sed "s:.*/::")

##==== 1.1. get reads with SA from both Read 1 and 2
	$samtools view -@ $cpuBWA $bam_path/$subii.consolidated.bam \
		| awk '{if (!($2==0 || $2==16 || $2==147 || $2==163 || $2==99 || $2==83)){
				if ($0 ~ "SA:") {
					print $0 > "_sa.sam"
				} 
			}}'

	$samtools view -@ $cpuBWA -T $refGenome -bS _sa.sam > _sa.bam
	$bedtools bamtobed -cigar -i _sa.bam > _sa.bed0

	## bedtools uses 0-base, change to 1-base:
	awk '{start = $2+1; $2=start; print}' _sa.bed0 | tr ' ' '\t' > _sa.bed 

if [ -s _sa.bed ]; then
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
        awk '{n=length($10); print $1,n}' _sa.sam > _sa.len

##==== 2. Join bed and read length, and keep the longest
	# By default, number and order of reads in files _sa.sam, _sa.bam and _sa.bed0 are identical. So, can use paste.
    	paste _sa.len _sa.bed | tr ' ' '\t' | cut -f2- > _sa.len.bed
		
	cut -f1,5 _sa.len.bed > _len.readID1
	sort -k2,2b -k1,1nr _len.readID1 > _len.readID2
	sort -k2,2b -u _len.readID2 > _len.readID

	sort -k2,2b _len.readID > _len.readID.s
	sort -k4,4b _sa.bed > _sa.bed.s
	join -1 2 -2 4 _len.readID.s _sa.bed.s > _sa.len.bed

##==== 3. remove MQ low
	echo | awk -v minMQ=$minMQ '{OFS="\t"; if ($6 >= minMQ) print $0}' _sa.len.bed > _sa.len.bed.mq

##==== 4. calculate query start, end
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

##==== 5. separate left and right split alignments
	awk '{if ($1==pre1){
			n=n+1
		} else {n=1};

		if (n==1) {
			print $0 > "_left";
			if (pren !=1) {print pre0 > "_right"}
		} else if (n>2){
			print pre0 > "split.mid"
		};
		pre1=$1; pren=n; pre0=$0;
		print n,$0 > "_sa.SMH4sn"
	}' _sa.SMH4s

##==== 6. get the 4 positions on a SA query read (after left and right alignments are merged by ReadID):
		## skematic drawing:
		##  left.query.start....left.query.end-[breakpoint]-right.query.start....right.query.end
		##         1                   2      -[breakpoint]-        3                  4
			## 2: left.query.end [corresponding mapping position in the genome is breakpoint 1]
			## 3: right.query.start [corresponding mapping position in the genome is breakpoint 2]
	        ## add breadkpoint (chr and pos)
            awk '{if ($7 =="+") {
                                print $0,$3"_"$5
                        } else {
                                mstart=$5; mend=$4; $4=mstart; $5=mend;
                                print $0,$3"_"$5
                        }
                }' _left | tr ' ' '\t' > _lefts

            awk '{if ($7 =="+") {
                                print $0,$3"_"$4
                        } else {
                                mstart=$5; mend=$4; $4=mstart; $5=mend;
                                print $0,$3"_"$4
                        }
                }' _right | tr ' ' '\t' | sed 1d > _rights

	##=== join leftmost and rightmost
        join _lefts _rights > _left_right.j

##==== 7. breakpoint.candidates.preFilter
	## left: 2-11; right 12-21
	awk '{if ($11 < $21){
			pp=$11"__"$21
		} else {pp=$21"__"$11
		}; 
		print $0,pp
	}' _left_right.j > _breakpoint.noFilter1
	sort -k1,1b _breakpoint.noFilter1 > _breakpoint.noFilter2
else
	touch breakpoint.candidates.preFilter
fi
##==== 8. correct breakpoint for those containing mid
if [ -f split.mid ]; then
	cut -f1 split.mid | sort -u > _mid.id
	sort -k2,2b _sa.SMH4sn > _sa.SMH4sns
	join -1 1 -2 2 _mid.id _sa.SMH4sns > split.mid.expanded

	if [ -s split.mid.expanded ]; then
		awk '{if ($1==pre1){
			diff = $5 - pre5;
			if ($4 != pre4 || diff > 100000 || diff < -100000){
				if (pre8=="+"){bkp1=pre4"_"pre6} else {bkp1=pre4"_"pre5};
				if ($8=="+"){bkp2=$4"_"$5} else {bkp2=$4"_"$6}
				q1=pre10; q2=pre11
			   };
			};
		if (bkp1 < bkp2){bkp = bkp1"__"bkp2};
		if (bkp1 > bkp2){bkp = bkp2"__"bkp1};
		if (bkp != ""){print $1,bkp,q1,q2,q3,q4 > "_sa.mid.bkp"}

		diff=""; bkp=""; bkp1=""; bkp2=""; q1=""; q2=""
		pre1=$1; pre4=$4; pre5=$5; pre6=$6; pre8=$8; pre10=$10; pre11=$11
		}' split.mid.expanded

		sort -k1,1b _sa.mid.bkp > _sa.mid.bkps
		join -a1 _breakpoint.noFilter2 _sa.mid.bkps > _breakpoint.noFilter3

		awk '{if (NF==25) {$9=$24; $10=$25; $22=$23}; print}' _breakpoint.noFilter3 | cut -d ' ' -f 1-22 > _breakpoint.noFilter.bkp.corrected

		join _breakpoint.noFilter.bkp.corrected _mid.id > breakpoint.noFilter.w.mid
		join -v 1 _breakpoint.noFilter.bkp.corrected _mid.id > breakpoint.candidates.preFilter
	else 
        	cp _breakpoint.noFilter2 breakpoint.candidates.preFilter
	fi
else
        cp _breakpoint.noFilter2 breakpoint.candidates.preFilter
fi

# DONE breakpoint candidate
