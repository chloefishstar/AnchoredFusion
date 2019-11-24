#!/usr/bin/python
from __future__ import print_function
import sys
import re
import os
import argparse

def parseArgs():

    parser = argparse.ArgumentParser(description='Split-Fusion is a \
                    fast data analysis pipeline detects gene fusion based \
                    on split reads and/or paired-end reads.')
    parser.add_argument('--SplitFusionPath', required=True
                        , help="the path where Split-Fusion pipeline is installed [required]")
    parser.add_argument('--R', required=True
                        , help="the path of R [required]")
    parser.add_argument('--perl', required=True
                        , help="the path of perl [required]")
    parser.add_argument('--refGenome', required=True
                        , help="the path where human genome reference is stored [required]")
    parser.add_argument('--sample_id', required=True
                        , help="the sample name of running [required]")
    parser.add_argument('--bam_dir'
                        , help="the path where bam or fastq file is stored. [Kickstart] The bam file of the sameple_id (xxx.bam or xxx.consolidated.bam) will be used. Either fastq_dir or bam_dir should be specified")
    parser.add_argument('--fastq_dir'
                        , help="the path where fastq file is stored. Either fastq_dir or bam_dir should be specified")
    parser.add_argument('--r1filename'
                        , help="Read 1 fastq filename. Can be in gzipped format. If not specified, $fastq_dir/$sample_id.R1.fq will be used.")
    parser.add_argument('--r2filename'
                        , help="Read 2 fastq filename. Can be in gzipped format. If not specified, $fastq_dir/$sample_id.R2.fq will be used.")
    parser.add_argument('--output', required=True
                        , help="the path where output is stored [required]")
    parser.add_argument('--panel', required=False
                        , default='NA'
                        , help="the path where target genes panel file is stored")
    parser.add_argument('--steps', required=False
                        , default='1_fastq-bam,2_bam-breakpoint,3_breakpoint-filter,4_breakpoint-anno,5_breakpoint-anno-post'
			, help="specify steps to run")
    parser.add_argument('--AnnotationMethod', required=False
                        , default = 'annovar'
                        , help="the name of annotation tools (annovar or snpEff)")
    parser.add_argument('--snpEff_ref', required=False
                        , default='GRCh37.75'
                        , help="the version of snpEff reference")
    parser.add_argument('--thread', type=int
                        , default=4
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
                        , default=18
                        , help="minimum read mapping length")
    parser.add_argument('--maxOverlap', type=int
                        , default=12
                        , help="maximum overlap length")
    parser.add_argument('--minExclusive', type=int
                        , default=20
                        , help="minimum exclusive length")
    parser.add_argument('--StrVarMinStartSite', type=int
                        , default=1
                        , help="minimum number of Adaptor Ligation Read Starting Sites to call Structure Variation/Fusion")
    parser.add_argument('--minFusionUniqReads', type=int
                        , default=2
                        , help="minimum number of reads to call Structure Variation/Fusion")

    args = vars(parser.parse_args())


    parser.add_argument('--snpEff', required=False
                        , default = args['SplitFusionPath'] + '/data/Database/snpEff/'
                        , help="the path of snpEff")
    parser.add_argument('--annovar', required=False
                        , default = args['SplitFusionPath'] + '/data/Database/annovar/'
                        , help="the path of snpEff")
    parser.add_argument('--fusion_library', required=False
			, default = args['SplitFusionPath'] + '/data/'
                        , help="the path where fusion library file is stored")
    parser.add_argument('--samtools', required=False
			, default = args['SplitFusionPath'] + '/data/Database/samtools'
                        , help="the path of samtools")
    parser.add_argument('--bedtools', required=False
			, default = args['SplitFusionPath'] + '/data/Database/bedtools'
                        , help="the path of bedtools")
    parser.add_argument('--java', required=False
			, default = args['SplitFusionPath'] + '/data/Database/jre1.8.0_201/bin/java'
                        , help="the path of java")
    parser.add_argument('--bwa', required=False
			, default = args['SplitFusionPath'] + '/data/Database/bwa-0.7.17/bwa'
                        , help="the path of bwa")
    
	
    args = vars(parser.parse_args())

    return args

def mkdir(path):
 
	folder = os.path.exists(path)
 
	if not folder:
		os.makedirs(path)

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

design =  args['R'] +  " -e " + "'suppressMessages(library(SplitFusion)); runSplitFusion(configFile = " + "\"" + config_o + "\"" + ")'" #+ "> /dev/null 2>&1"
os.system(design)

## END
