#!/bin/bash
. ../../run.info.sh

##==== Filter 1: minMapLength, minExclusive, maxQueryGap, maxOverlap ====:

		## reads with middle split
		cut -f1 split.mid | sort -u | awk '{print $0,"_mid"}' > _mid.id
		
		join -a1 breakpoint _mid.id > breakpoint2
		grep '_mid$' breakpoint2 | sed 's:_mid$::' > breakpoint.w.mid
		grep -v '_mid$' breakpoint2 > breakpoint.wo.mid

		## 1. at least minMapLength 
		## 2. at least minExclusive bp exculsive between hit1 and hit2
		## 3. no larger than maxQueryGap
		## 4. less than max overlapping length

		## reads without middle split
		    ## left:  $10-------$11
		    ## right:       $24------$25
		 echo | awk -v minMapLength=$minMapLength -v minExclusive=$minExclusive -v maxQueryGap=$maxQueryGap -v maxOverlap=$maxOverlap \
		    '{ gap = $24-$11-1; 
		       overlap = $11-$24+1;
			    if (  ($11-$10 >= minMapLength && $25-$24 >= minMapLength) \
				    && ($24-$10 >= minExclusive && $25-$11 >= minExclusive && gap <= maxQueryGap && overlap <= maxOverlap) \
			) {print $0,overlap}}' \
			     breakpoint.wo.mid > _sa.fu0

		## reads with middle split, turn off maxQueryGap by let = 100
		    ## left:  $10-------$11
		    ## right:		       $24------$25
		 echo | awk -v minMapLength=$minMapLength -v minExclusive=$minExclusive -v maxQueryGap=100 \
		    '{ gap = $24-$11-1; 
		       overlap = $11-$24+1;
			    if (  ($11-$10 >= minMapLength && $25-$24 >= minMapLength) \
				    && ($24-$10 >= minExclusive && $25-$11 >= minExclusive && gap <= maxQueryGap && overlap <= maxOverlap) \
			) {print $0,overlap}}' \
			     breakpoint.w.mid >> _sa.fu0

			## keep fusionID, gap, overlap
			cut -d ' ' -f 32,33 _sa.fu0 > fusionID.overlap


##==== Filter 2: StrVarMinUniqMap
## Default in run.info.sh: StrVarMinUniqMap=2
    ## or uncomment below and define the value:
    # StrVarMinUniqMap=new.value

	## stagger sites on left $4
	# sort by point-point, left-breakpoint, left-stagger.point
	sort -k32,32b  -k30,30b -k4,4n _sa.fu0 > _sa.fu0s
	sort -k32,32b  -k30,30b -k4,4n -u  _sa.fu0s > _fu_1u	

		## if stagger sites >= minUniqMap
                echo | awk -v minUniqMap=$StrVarMinUniqMap '{
			if ($32==pre32 && $30==pre30){
				if ($4 - pre4 >1){
					cnt=cnt+1
				}
			} else { cnt = 1};
                        pre32=$32; pre30=$30; pre4=$4;
			print $32,cnt,"stgLeft"}' _fu_1u | sort -k1,1b -k2,2nr | sort -k1,1b -u | awk '{if ($2 != 1) print $0}' > fu.breakpointID.stgLeft
		join -2 32 fu.breakpointID.stgLeft _sa.fu0s > fu.info.stgLeft	
	
	## stagger sites on right $19
	sort -k32,32b -k31,31b -k19,19n -u  _sa.fu0s > _fu_2u	

		## if stagger sites >= minUniqMap
                echo | awk -v minUniqMap=$StrVarMinUniqMap '{
			if ($32==pre32 && $31==pre31) {
				if ($19 - pre19 >1){
                            		cnt=cnt+1
				}
			} else {cnt = 1};
                        pre32=$32; pre31=$31; pre19=$19;
			print $32,cnt,"stgRight"}' _fu_2u | sort -k1,1b -k2,2nr | sort -k1,1b -u | awk '{if ($2 != 1) print $0}' > fu.breakpointID.stgRight
		join -2 32 fu.breakpointID.stgRight _sa.fu0s > fu.info.stgRight	


	##=== combine fusionID
	cat fu.info.stgLeft fu.info.stgRight | sort -k1,1b  > fusion.candidates
