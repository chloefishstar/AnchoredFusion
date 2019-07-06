#!/usr/bin/python
from __future__ import print_function
import sys
import re
import os
import argparse

def parseArgs():

    parser = argparse.ArgumentParser(description='Anchored-Fusion is a fast data analysis 
                    pipeline detects gene fusion based on split reads \
                    and/or paired-end reads.')
    parser.add_argument('--AnchoredFusionPath', required=True
                        , default='/tools/repo/anchored-fusion'
                        , help="the path where Anchored-Fusion pipeline is installed.")
    parser.add_argument('--sampleInfo', required=True
                        , default=''
                        , help="...")
    parser.add_argument('--minMapLength', type=int, default=25, help="minimum read mapping length")

    args = vars(parser.parse_args())

    return args

cleanup = "rm config.txt"
os.system(cleanup)
if __name__ == '__main__':
    
    args = parseArgs()
    config = open("config.txt", "w+")
    print(args, file = config, sep="\n")
    config.close()

design =  args['AnchoredFusionPath'] + "/main.R"
os.system(design)

## END
