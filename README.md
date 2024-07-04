# bamdam: A toolkit for ancient environmental DNA shotgun sequencing data
Written by Bianca De Sanctis in July 2024. Please contact me with any issues or requests! I want to improve this software! bddesanctis@gmail.com

Welcome to bamdam! The goal of this toolkit is to provide functionality after shotgun sequencing aeDNA reads have been mapped to a reference database and run through ngsLCA. The main program is BamDam.py, with additional plotting functionality in extra R scripts. 

There is no installation required - just download this Github repo, make sure you have the python and R libraries installed (they should be pretty standard) and run. 

/ Main Script BamDam.py / 

The main script BamDam.py will:
(1) Shorten the LCA file to only nodes which are equal to or below your tax threshold, meet your minimum read count (or are below nodes which meet your minimum read count). Also, optionally include all lines with a user-defined list of strings (such as "Eukaryota"), and optionally exclude all lines with a user-defined list of strings (such as "Homo sapiens").
(2) Shorten the bam file to include only reads which appear in the newly shortened LCA file. This can make 100-1000x or even more difference in the size of a bam file after typical ancient eDNA mapping steps. It will also annotate the new bam file with PMD scores as in PMDTools (in the DS:Z field).
(3) Gather various statistics per node and write to a stats file, sorted by total reads assigned. More information detailing this file is below.
(4) Gather substitution patterns per node and write to a subs file, sorted in the same order as the stats file (so line n of the stats file should correspond to line n of the subs file).
You can then look through your stats file and decide what seems real and ancient to you. You may also plot damage plots for individual tax IDs using an R script which is detailed in the next section.

Please ensure your input bam is sorted in the same order as your input lca file, and that every read in your LCA file is also in your bam file (the opposite need not be true). This will be true by default if you use fresh ngsLCA output. Also, for now, please provide a ngsLCA lca file, not a metaDMG lca file (they have different formats).

Usage: BamDam.py [-h] --in_lca IN_LCA --in_bam IN_BAM --out_lca OUT_LCA --out_bam OUT_BAM --out_stats OUT_STATS --out_subs OUT_SUBS --stranded STRANDED
                 [--mincount MINCOUNT] [--k K] [--upto UPTO] [--lcaheaderlines LCAHEADERLINES] [--minsim MINSIM]
                 [--keep_keywords KEEP_KEYWORDS [KEEP_KEYWORDS ...]] [--exclude_keywords EXCLUDE_KEYWORDS [EXCLUDE_KEYWORDS ...]]
                 [--exclude_keyword_file EXCLUDE_KEYWORD_FILE]

options:
  -h, --help            show this help message and exit
  --in_lca IN_LCA       Path to the original (sorted) LCA file (required)
  --in_bam IN_BAM       Path to the original (sorted) BAM file (rquired)
  --out_lca OUT_LCA     Path to the short output LCA file (required)
  --out_bam OUT_BAM     Path to the short output BAM file (required)
  --out_stats OUT_STATS
                        Path to the output stats file (required)
  --out_subs OUT_SUBS   Path to the output subs file (required)
  --stranded STRANDED   Either ss for single stranded or ds for double stranded (required)
  
  --mincount MINCOUNT   Minimum read count to keep a node (default: 5)
  --k K                 Value of k for kmer complexity calculations (default: 5)
  --upto UPTO           Keep nodes up to and including this tax threshold, use root to disable (default: family)
  --lcaheaderlines LCAHEADERLINES
                        Number of header lines in input LCA file (default: 0)
  --minsim MINSIM       Minimum similarity to reference to keep a read (default: 0.95)
  --keep_keywords KEEP_KEYWORDS [KEEP_KEYWORDS ...]
                        Other keyword(s) in LCA file for filtering to keep, e.g. Eukaryota (default: none)
  --exclude_keywords EXCLUDE_KEYWORDS [EXCLUDE_KEYWORDS ...]
                        Other keyword(s) in LCA file for filtering to delete (default Hominidae)- NOT IMPLEMENTED YET
  --exclude_keyword_file EXCLUDE_KEYWORD_FILE
                        Text file containing list of keywords or tax paths to remove in LCA file, one on each line- NOT IMPLEMENTED YET

BamDam was written by Bianca De Sanctis in 2024. For assistance please contact bddesanctis@gmail.com.

