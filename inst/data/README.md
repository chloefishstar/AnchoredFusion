The dependency data (e.g. in 'data') should contain:

| Filename                   | Content                                                        |
|:--------------------------:|:--------------------------------------------------------------:|
| Homo_sapiens_assembly19.fasta | Contains a list of human genome reference, please mannually downloaded from ucsc or other official site. |
| panel-name.target.genes.txt    | Contains a list of targets (gene name, e.g. ALK, ROS1, etc.), e.g: ITFTNA.target.genes.txt   |
| fusion.gene-exon.filter.txt | Contains recurrent breakpoints identified as data accumulates, but are not of interest. |
| fusion.gene-exon.txt          | By default, SplitFusion only outputs fusions that are *in-frame* fusion of two different genes or when number of breakpoint-supporting reads exceed predefined threashold. This file contains known breakpoints that do not belong to the above two kinds, but are clinically relevant, e.g. "MET_exon13---MET_exon15" an exon-skipping event forms an important theraputic target. Many exon skipping/alternative splicing events are normal or of unknown clinical relevance and are thus not output by default. |
| fusion.partners.txt | Contains a list of known fusion partners of targets. |
| ENSEMBL.orientation.txt | Due to the lack of transcript orientation in snpEff annotation, so this file include two columns, Orientation (+ or –) and transcript ID (ENST*). |


The above files could be updated periodically as a backend supporting database that facilitates automatc filtering and outputing of fusion candidates.
