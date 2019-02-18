#!/bin/bash
. $1
subii=$( pwd | sed "s:.*/::")
SA_flag=$2

#== consolidated bam ==
##get threads
head -n 1 $sampleInfo > _head.1
        gawk '{for (i=1; i<=NF; i++) {if ($i ~ /cpuBWA/) {print $i,i}}}' _head.1 | sed 's/.* /cpuField=/' > _t.sh
        . _t.sh
        cpuBWA=$(grep $subii $sampleInfo | cut -f $cpuField | sed 's: .*::')

if [ ! -s $bam_path/$subii.consolidated.bam ]; then

	#== fasq2bam ==
	if [ ! -s $bam_path/$subii.bam ]; then
		$bwa mem -T 20 -5a -t 5 $hgRef $bam_path/$subii.R1.fq $bam_path/$subii.R2.fq > $bam_path/$subii.sam 2> bwa.log
		$samtools view -@ $cpuBWA -T $hgRef -bS $bam_path/$subii.sam > $bam_path/$subii.bam
	fi

	if [ "$SA_flag"x = "SA"x ]; then 
		$samtools view -@ $cpuBWA $bam_path/$subii.bam | grep 'SA:' > $bam_path/$subii.sam
		$samtools view -@ $cpuBWA -T $hgRef -bS $bam_path/$subii.sam > $bam_path/$subii.bam
	fi
	
	#== processed UMI ==
	if [ ! -s $bam_path/$subii.UMI.bam ]; then

		#== raw bed ==
#        	$samtools view -@ $cpuBWA -T $hgRef -bS $bam_path/$subii.sam > $bam_path/$subii.bam
                $samtools flagstat $bam_path/$subii.bam > flag.stat
                $samtools view -H $bam_path/$subii.bam > header
		
		if [ ! -s $bam_path/$subii.sam ]; then
			$samtools view -@ $cpuBWA $bam_path/$subii.bam > $bam_path/$subii.sam
		fi

                $bedtools bamtobed -cigar -i $bam_path/$subii.bam > _raw.bed

		#=========================================
		#=== add chr.pos (at ligation site) to umi
		#=========================================
        	sed -e 's/:umi:/\tumi\t/' _raw.bed |\
                gawk '{OFS="\t";
                        if ($4 != preID){
                                if ($8 == "+"){
                                        posC = 100000001 + $2
                                } else {
                                        posC = 100000001 + $3
                                }
                                umi="C"$1"P"posC
                        } else {
                                umi=preUmi
                        };
 
                        preUmi=umi;
                        preID=$4;
 
                        if ($6 ~ /N/){
                                print $4":umi",umi > "_bed.umiN"
                        } else {
                                print $4":umi",umi > "_bed.umi"
                        }
                }' 2>/dev/null
 
	        uniq _bed.umi > _bed.umi.u
        	wc -l _bed.umi* > umiN
        	sort --parallel=$cpuBWA -k1,1b _bed.umi.u > _bed.umi.us
 
		#=== UMI BAM ===
		#=== join Cleaned umi with raw sam ===
        	sed -e 's:\t\t:\t*\t:g' -e 's/:umi:/:umi\t/' $bam_path/$subii.sam > _raw.samC
                rm $bam_path/$subii.sam
        	sort --parallel=$cpuBWA -k1,1b _raw.samC > _raw.sam.s
                rm _raw.samC
 
        	cp header _umi.sam
        	join _bed.umi.us _raw.sam.s | sed -e 's/ /:/' -e 's/ /-/' -e 's/ /\t/g' >> $bam_path/$subii.UMI.sam
                rm _raw.sam.s
 
	        $samtools view -@ $cpuBWA -T $hgRef -bS $bam_path/$subii.UMI.sam > _umi.bam
        	$samtools sort -@ $cpuBWA _umi.bam -o $bam_path/$subii.UMI.bam
		rm _umi.bam
#	        $samtools index $bam_path/$subii.UMI.bam $bam_path/$subii.UMI.bam.bai
	fi
	
	#== de-duplication ==
	#== consolidated UMI BAM==
        # prep data
        oriDir="$(pwd)"
        ramtmp="$(mktemp -dp .)"
        cd $oriDir/$ramtmp
        
	#$REPPATH/samtools/bin/samtools view --threads $cpuBWA $oriDir/../$subii.raw.bam | sed 's/::umi:/\t/' | awk '{OFS="\t"; if ($4 != "") {print $0 > $4}}'
        if [ ! -s $bam_path/$subii.UMI.sam ]; then
                sed -e 's/\t/_flag_/'  -e 's/::umi:/\t/' $bam_path/$subii.UMI.sam | awk '{OFS="\t"; if ($3 != "") {print $0 > $3}}'
        else
                $samtools view -@ $cpuBWA $bam_path/$subii.UMI.bam |sed 's/\t/_flag_/' | sed 's/::umi:/\t/' | awk '{OFS="\t"; if ($3 != "") {print $0 > $3}}'
        fi

        # keep 1st umi-flag
        # sort chr.pos_umi_flag in each chr
        find ./ -type f | parallel "sort -k2,2b -u {} | sed -e 's/\t/:umi:/' -e 's/_flag_/\t/' > {}.consolidated"
 
        cd $oriDir
        cp header consolidated.sam
        cat $ramtmp/*.consolidated >> consolidated.sam
        rm -rf $oriDir/$ramtmp
 
        $samtools view -@ $cpuBWA -T $hgRef -bS consolidated.sam > _consolidated.bam
        rm consolidated.sam
        $samtools sort -@ $cpuBWA _consolidated.bam -o $bam_path/$subii.consolidated.bam
        rm _consolidated.bam
	rm $bam_path/$subii.UMI.sam

#	$samtools index $bam_path/$subii.consolidated.bam $bam_path/$subii.consolidated.bam.bai
fi

		
