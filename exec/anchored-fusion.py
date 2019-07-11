#!/usr/bin/python
from __future__ import print_function
import sys
import re
import os
import argparse

def parseArgs():

    parser = argparse.ArgumentParser(description='Anchored-Fusion is a \
                    fast data analysis pipeline detects gene fusion based \
                    on split reads and/or paired-end reads.')
    parser.add_argument('--AnchoredFusionPath', required=True
                        #, default=''
                        , help="the path where Anchored-Fusion pipeline is installed [required]")
    parser.add_argument('--R', required=True
                        #, default=''
                        , help="the path of R [required]")
    parser.add_argument('--hgRef', required=True
                        #, default=''
                        , help="the path where human genome reference is stored [required]")
    parser.add_argument('--bam_path', required=True
                        #, default=''
                        , help="the path where bam or fastq file is stored [required]")
    parser.add_argument('--sample_id', required=True
                        #, default=''
                        , help="the sample name of running [required]")
    parser.add_argument('--output', required=True
                        #, default=''
                        , help="the path where output is stored [required]")
    
    parser.add_argument('--panel', required=False
                        , default='NA'
                        , help="the path where target genes panel file is stored")
    parser.add_argument('--fusion_library', required=False
                        #, default=''
                        , help="the path where fusion library file is stored")
    parser.add_argument('--step', required=False
                        , default='bam-consolidate,breakpoint-consolidate,breakpoint-filter,breakpoint-anno,breakpoint-anno-post'
                        , help="the step of running")

    parser.add_argument('--samtools', required=False
                        #, default=''
                        , help="the path of samtools")
    parser.add_argument('--bedtools', required=False
                        #, default=''
                        , help="the path of bedtools")
    parser.add_argument('--java', required=False
                        #, default=''
                        , help="the path of java")
    parser.add_argument('--bwa', required=False
                        #, default=''
                        , help="the path of bwa")
    parser.add_argument('--snpEff', required=False
                        #, default=''
                        , help="the path of snpEff")
    parser.add_argument('--snpEff_ref', required=False
                        , default='hg19'
                        , help="the version of snpEff reference")

    parser.add_argument('--cpuBWA', type=int
                        , default=2
                        , help="threads of BWA")
    parser.add_argument('--strVarMinStartSite', type=int
                        , default=3
                        , help="minimum start site")
    parser.add_argument('--maxQueryGap', type=int
                        , default=0
                        , help="maximum gap length")
    parser.add_argument('--minMQ', type=int
                        , default=13
                        , help="minimum mapping quality")
    parser.add_argument('--minMapLength', type=int
			, default=25
			, help="minimum read mapping length")
    parser.add_argument('--maxOverlap', type=int
                        , default=9
                        , help="maximum overlap length")
    parser.add_argument('--minExclusive', type=int
                        , default=25
                        , help="minimum exclusive length")
    
    
	

    args = vars(parser.parse_args())

    return args


def mkdir(path):
 
	folder = os.path.exists(path)
 
	if not folder:
		os.makedirs(path)

#cleanup = "rm config.txt"
#os.system(cleanup)
if __name__ == '__main__':
    
    args = parseArgs()
    output=[str(k) + "=" + "\"" + str(v)+ "\"" for k,v in args.iteritems() if v != None]
    config_p = args['output'] + "/" + args['sample_id'] + "/" 
    mkdir(config_p)
    config_o = config_p + args['sample_id'] + ".config.txt"
    config = open(config_o, "w+")
    #print(args, file = config, sep="\n")
    print("\n".join(output), file = config, sep="\n")

    config.close()

#config_run = "library(AnchoredFusion);runAnchoredFusion(runInfo = \"config.txt\")"
#design =  args['AnchoredFusionPath'] + "/data/Database/R" +  " -e " + config_run
design =  args['R'] +  " -e " + "'library(AnchoredFusion);runAnchoredFusion(runInfo = " + "\"" + config_o + "\"" + ")'"
os.system(design)

## END
