#!/bin/bash

. $1
subii=$( pwd | sed "s:.*/::")

	#== If specify fastq_dir, do fasq to bam ==
	if [ "$fastq_dir" != "" ]; then
		if [ "$r1filename" != "" ]; then
			$bwa mem -T 20 -5a -t $thread $refGenome $fastq_dir/$r1filename $fastq_dir/$r2filename > _raw.sam 2> bwa.log
		else	
			$bwa mem -T 20 -5a -t $thread $refGenome $fastq_dir/$subii.R1.fq $fastq_dir/$subii.R2.fq > _raw.sam 2> bwa.log
		fi
	elif [ "$bam_dir" != "" ]; then
	    #== If specify bam_dir, start from bam ==
	    # SplitFusion requires that the 5' end of Read_1 (ligation site) is integrated into UMI
	    #   so to calculate the number of unique ligation site supporting a fusion in the final step.
	    # If existing bam file was not mapped using the '-5a' option during BWA MEM, 
	    #   we recommend to specify fastq-dir only and leave bam_dir unspecified so that the pipeline
	    #   starts from fastq mapping using BWA MEM with the -5a option.

	    $samtools view -@ $thread $bam_dir/$subii.bam > _raw.sam
	else 
		echo "Must specify fastq_dir or bam_dir"
		exit
	fi

	head -n 1000 _raw.sam | grep -v ^@ > header

#== Consolidate reads based on unique UMI

	#==== if no umi, add umi:A for compatability
        hasUmi=$(grep -v ^@ _raw.sam | head -n 1 | cut -f 4 | grep umi: | wc -l)
		if [ $hasUmi -eq 0 ]; then
			grep -v ^@ _raw.sam | sed 's/\t/:umi:A\t/' > _raw.sam2
			mv _raw.sam2 _raw.sam
		fi

	$samtools view -@ $thread -T $refGenome -bS _raw.sam > _raw.bam

	#==== Prepare UMI: add chr.pos (at ligation site) to umi
                $bedtools bamtobed -cigar -i _raw.bam > _raw.bed

		#==== Reformat if not already in the required UMI format, which has ligation site integrated in UMI
	        goodformat=$(head -n 1 _raw.bed | cut -f 4 | grep umi: | sed 's:.*umi:umi:' | grep C | grep P | wc -l)
       		if [ $goodformat -eq 0 ]; then
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
		    sort --parallel=$thread -k1,1b _bed.umi.u > _bed.umi.us
     
		    #=== join Cleaned umi with raw sam by Read ID ===
		    sed -e 's:\t\t:\t*\t:g' -e 's/:umi:/:umi\t/' _raw.sam > _raw.samC
		    sort --parallel=$thread -k1,1b _raw.samC > _raw.sam.s
			    rm _raw.samC
     
		    join _bed.umi.us _raw.sam.s > _umi.sam0
		    sed -e 's/ /:/' -e 's/ /-/' -e 's/ /\t/g' _umi.sam0 > _umi.sam
			    rm _raw.sam.s
		fi
	
	#== Parallel de-duplication / consolidate umi bam ==
	        oriDir="$(pwd)"
       		ramtmp="$(mktemp -dp .)"
        	cd $oriDir/$ramtmp
        
		sed -e 's/\t/_flag_/'  -e 's/:umi:/\t/' $oriDir/_umi.sam | grep -ve '^@SQ' -ve '^@PG' | awk '{OFS="\t"; if ($3 != "") {print $0 > $3}}'

	        # keep 1st umi-flag
       		# sort chr.pos_umi_flag in each chr
	        find ./ -type f | parallel "sort -k2,2b -u {} | sed -e 's/\t/:umi:/' -e 's/_flag_/\t/' > {}.consolidated"
 
	        cd $oriDir
	        cat $ramtmp/*.consolidated > consolidated.sam
	        rm -rf $oriDir/$ramtmp
		rm _*
		grep 'SA:' consolidated.sam > _sa.sam
 
        $samtools view -@ $thread -T $refGenome -bS consolidated.sam > _consolidated.bam
        $samtools sort -@ $thread _consolidated.bam -o $subii.consolidated.bam
	        rm consolidated.sam _consolidated.bam
	$samtools index $subii.consolidated.bam 

