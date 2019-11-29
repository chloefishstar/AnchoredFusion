#!/bin/bash
. $1

SampleId=$( pwd | sed "s:.*/::")

##==== Filter 1: minMapLength, minExclusive, maxQueryGap, maxOverlap ====:
		## 1. at least minMapLength 
		## 2. at least minExclusive bp mapped exculsively to alignments 1 and 2
		## 3. no larger than maxQueryGap
		## 4. no larger than max overlapping length

		## reads without middle split
		    ## left:  $9-------$10
		    ## right:       $19------$20
		echo | awk -v minMapLength=$minMapLength -v minExclusive=$minExclusive -v maxQueryGap=$maxQueryGap -v maxOverlap=$maxOverlap \
		    '{ gap = $19-$10-1; 
			   overlap = $10-$19+1;
				if (  ($10-$9 >= minMapLength && $20-$19 >= minMapLength) \
					&& ($19-$9 >= minExclusive && $20-$10 >= minExclusive && gap <= maxQueryGap && overlap <= maxOverlap) \
			    ) {print $0,overlap} else {print $1 > "_filter1"}
			}' breakpoint.candidates.preFilter > _sa.fu01

		## reads with middle split, turn off maxQueryGap by let maxQueryGap=100
		    ## left:  $9-------$10
		    ## right:		       $19------$20
		if [ -f breakpoint.noFilter.w.mid ]; then
		 echo | awk -v minMapLength=$minMapLength -v minExclusive=$minExclusive -v maxQueryGap=1000 \
		    '{ gap = $19-$10-1; 
			   overlap = $10-$19+1;
				if (  ($10-$9 >= minMapLength && $20-$19 >= minMapLength) \
					&& ($19-$9 >= minExclusive && $20-$10 >= minExclusive && gap <= maxQueryGap) \
			    ) {print $0,overlap} else {print $1 > "_filter2"}
			}' breakpoint.noFilter.w.mid > _sa.fu02

			join -v 1 _sa.fu02 _filter2 > _sa.fu02f

			# filter gap and overlap
			echo | awk -v maxQueryGap=$maxQueryGap -v maxOverlap=$maxOverlap \
			'{ if ($1==pre1){
					gap = $10 - pre11 -1; overlap = pre11 - $10 + 1
					if (gap > maxQueryGap || overlap > maxOverlap) {print $1}
				} else {gap=0; overlap=0};
				pre1=$1; pre11=$11;
			}' split.mid.expanded > _mid.gap.id
			sort --parallel=$thread -u _mid.gap.id > _mid.gap.id.u

			join -v 1 _sa.fu02f _mid.gap.id.u > _sa.fu02ff 
			cat _sa.fu01 _sa.fu02ff > _sa.fu0
		else cp _sa.fu01 _sa.fu0
		fi

##==== Filter 2: FusionMinStartSite
	## breakpoint ($22 now, later $23) and separate start site (chr+pos)
	sed 's/:umi:/\t/' _sa.fu0 | tr ' ' '\t' | awk '{OFS="\t"; print $23,$2,$0}' | sed -e 's/C\([^\t]\+\)P\([0-9]\+\)-/\1\t\2\t/' -e 's:/[12]::' > _sa.fu2

	# sort by breakpoint and  start.site.umi
	sort --parallel=$thread -k1,1b -k28,28n -k6,6b _sa.fu2 > _sa.fu3

	## breakpoint stats: num_unique_molecule (numi), num_start_site (nss), num_start_site2 (diff by at least 2, nss2), average MQ (avgMQ1, avgMQ2)
        awk '{OFS="\t";
		if ($1 == pre1 && $2 == pre2 && $28 == pre28){
			i += 1;
			diff = $3-pre3;
			if (diff < 750000){
				siteID = preSiteID
				if (diff==0){
					if ($4 != pre4){
						numi += 1
					}
				} else {
					numi += 1;
					nss += 1;
					if (diff >1){nss2 += 1}
				}
			} else {
				siteID = NR
				numi=1; nss=1; nss2=1
			};

			mq1 += $11;
			mq2 += $21;
		} else {
				i=1; numi=1; nss=1; nss2=1; siteID=NR; mq1=$11; mq2=$21
		};

		avgMQ1 = mq1/i; avgMQ2 = mq2/i;
		print siteID,numi,nss,nss2,avgMQ1,avgMQ2,$0 > "breakpoint.stats";
     		pre1=$1; pre2=$2; pre3=$3; pre4=$4; pre28=$28; preSiteID=siteID
	}' _sa.fu3

	##====  Apply Filter2, and
	## min start site step size of 2 (Deprecated: set to 1)
		minStartStepSize=1
	cut -f1-6 breakpoint.stats | tac > _breakpoint.stats6
	sort --parallel=$thread -k1,1b -u _breakpoint.stats6 > _breakpoint.stats6.u
	echo | awk -v FusionMinStartSite=$FusionMinStartSite -v minStartStepSize=$minStartStepSize -v minAvgMQ=$minAvgMQ \
		'{OFS="\t"; if ($3 >= FusionMinStartSite && $4 >= minStartStepSize && $5 >= minAvgMQ && $6 >= minAvgMQ){
			print $1,$2,$3,$4 > "breakpoint.siteID.stats.filtered"
			}
		}' _breakpoint.stats6.u

	sort --parallel=$thread -k1,1b breakpoint.stats > _breakpoint.stats.s
	join -1 1 -2 1 breakpoint.siteID.stats.filtered  _breakpoint.stats.s \
		| sed 's/ /:umi:/14' \
		| tr ' ' '\t' | cut -f 2,3,4,14- > breakpoint.candidates

touch _0; rm _*

## End: breakpoint candidates