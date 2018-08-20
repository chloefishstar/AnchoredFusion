#!/bin/bash
. ../../run.info.sh

##==== Filter 1: minMapLength, minExclusive, maxQueryGap, maxOverlap ====:

		## reads with middle split
		cut -f1 split.mid | sort -u | awk '{print $0,"_mid"}' > _mid.id
		
		join -a1 breakpoint.noFilter _mid.id > breakpoint.noFilter2
		grep '_mid$' breakpoint.noFilter2 | sed 's:_mid$::' > breakpoint.noFilter.w.mid
		grep -v '_mid$' breakpoint.noFilter2 > breakpoint.noFilter.wo.mid

		## 1. at least minMapLength 
		## 2. at least minExclusive bp exculsive between alignments 1 and 2
		## 3. no larger than maxQueryGap
		## 4. less than max overlapping length

		## reads without middle split
		    ## left:  $10-------$11
		    ## right:       $21------$22
		echo | awk -v minMapLength=$minMapLength -v minExclusive=$minExclusive -v maxQueryGap=$maxQueryGap -v maxOverlap=$maxOverlap \
		    '{ gap = $21-$11-1; 
		       overlap = $11-$21+1;
			    if (  ($11-$10 >= minMapLength && $22-$21 >= minMapLength) \
				    && ($21-$10 >= minExclusive && $22-$11 >= minExclusive && gap <= maxQueryGap && overlap <= maxOverlap) \
			) {print $0,overlap}}' \
			     breakpoint.noFilter.wo.mid > _sa.fu0

		## reads with middle split, turn off maxQueryGap by let = 100
		    ## left:  $10-------$11
		    ## right:		       $21------$22
		 echo | awk -v minMapLength=$minMapLength -v minExclusive=$minExclusive -v maxQueryGap=1000 \
		    '{ gap = $21-$11-1; 
		       overlap = $11-$21+1;
			    if (  ($11-$10 >= minMapLength && $22-$21 >= minMapLength) \
				    && ($21-$10 >= minExclusive && $22-$11 >= minExclusive && gap <= maxQueryGap && overlap <= maxOverlap) \
			) {print $0,overlap}}' \
			     breakpoint.noFilter.w.mid >> _sa.fu0

##==== Filter 2: StrVarMinStartSite
	## breakpoint ($24) and separate start site (chr+pos)
	sed 's/:umi:/\t/' _sa.fu0 | tr ' ' '\t' | awk '{OFS="\t"; print $25,$2,$0}' | sed -e 's/C\([^\t]\+\)P\([0-9]\+\)-/\1\t\2\t/' -e 's:/[12]::' > _sa.fu2

	# sort by breakpoint and  start.site.umi
	sort -k1,1b -k6,6b _sa.fu2 > _sa.fu3

	## breakpoint stats: num_start_site (nss), num_unique_molecule (numi), num_start_site2 (diff by at least 2, nss2)
        awk '{OFS="\t";
		if ($1 == pre1 && $2 == pre2){
			diff = $3-pre3;
			if (diff < 100){
				siteID = preSiteID
				if (diff==0){
					if ($4 != pre4){
						numi = numi+1
					}
				} else {
					numi = numi+1;
					nss = nss+1;
					if (diff >1){
						nss2 = nss2+1
					}
				}
			} else {
				siteID = NR
				numi=1; nss=1; nss2=1
			};
		} else {
				numi=1; nss=1; nss2=1; siteID=NR
		}
		print siteID,numi,nss,nss2,$0 > "breakpoint.stats";
     		pre1=$1; pre2=$2; pre3=$3; pre4=$4; preSiteID=siteID
	}' _sa.fu3 

	## Apply Filter2
	cut -f1-4 breakpoint.stats | tac > _breakpoint.stats4
	sort -k1,1b -u _breakpoint.stats4 > _breakpoint.stats4.u
	echo | awk -v StrVarMinStartSite=$StrVarMinStartSite \
		'{if ($3 >= StrVarMinStartSite){
			print $0 > "breakpoint.siteID.MinStartSite"
			}
		}' _breakpoint.stats4.u

	sort -k1,1b breakpoint.stats > _breakpoint.stats.s
	join breakpoint.siteID.MinStartSite  _breakpoint.stats.s \
		| sed 's/ /:umi:/12' \
		| tr ' ' '\t' | cut -f 2,3,4,12- > breakpoint.candidates

## End: breakpoint candidates
