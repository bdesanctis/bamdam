# bamdam: A toolkit for ancient environmental DNA shotgun sequencing data

Written by Bianca De Sanctis in July 2024. Please contact me with any issues or requests! I want to improve this software! bddesanctis@gmail.com

Welcome to bamdam! The goal of this toolkit is to provide functionality after shotgun sequencing aeDNA reads have been mapped to a reference database and run through ngsLCA. The main program is BamDam.py, with additional plotting functionality in extra R scripts.

## Table of Contents
- [Installation](#installation)
- [Main Script: BamDam.py](#main-script)
- [Explanation of the Output Stats File Columns](#explanation-of-the-output-stats-file-columns)
- [Plotting Script](#plotting-script)
- [Coming Soon](#coming-soon)
- [License](#license)

## Installation

Just download this Github repo, make sure you have the python and R libraries installed (they should be pretty standard) and run it.
```sh
git clone https://github.com/bdesanctis/bamdam
```

## Main Script

Input: Bam and LCA files. Output: Stats and subs files.

The main script BamDam.py will:
1. Shorten the LCA file to only nodes which are equal to or below your tax threshold, meet your minimum read count (or are below nodes which meet your minimum read count). This is not implemented just yet, but soon I will add the option to include all lines with a user-defined list of strings (such as "Eukaryota"), and optionally exclude all lines with a user-defined list of strings (such as "Homo sapiens, Galus Galus"), or with a tax node ID.
2. Shorten the bam file to include only reads which appear in the newly shortened LCA file. This can make 100-1000x or even more difference in the size of a bam file after typical ancient eDNA mapping steps. It will also annotate the new bam file with PMD scores as in PMDTools (in the DS:Z field).
3. Gather various statistics per node and write to a stats file, sorted by total reads assigned. More information detailing this file is below.
4. Gather substitution patterns per node and write to a subs file, sorted in the same order as the stats file (so line n of the stats file should correspond to line n of the subs file).

You can then look through your stats file and decide what seems real and ancient to you. You may also plot damage plots for individual tax IDs using an R script and the subs file which is detailed in the next section. The subs file should be fairly self-explanatory in formatting, if you are interested in looking at it.

Please ensure your input bam is sorted in the same order as your input lca file, and that every read in your LCA file is also in your bam file, which will be true by default if you use fresh ngsLCA output. However, the opposite need not be true, so you can subset your lca file before running bamdam if you'd like, and bamdam will feed only those lines in the lca file through to the shortened bam. Also, for now, please provide a ngsLCA lca file, not a metaDMG lca file (they have different formats).

Usage and options: 
```sh
BamDam.py [-h] --in_lca IN_LCA --in_bam IN_BAM --out_lca OUT_LCA --out_bam OUT_BAM --out_stats OUT_STATS --out_subs OUT_SUBS --stranded STRANDED
          [--mincount MINCOUNT] [--k K] [--upto UPTO] [--lcaheaderlines LCAHEADERLINES] [--minsim MINSIM]
          [--keep_keywords KEEP_KEYWORDS [KEEP_KEYWORDS ...]] [--exclude_keywords EXCLUDE_KEYWORDS [KEEP_KEYWORDS ...]]
          [--exclude_keyword_file EXCLUDE_KEYWORD_FILE]
          
options:

-h, --help: show this help message and exit
--in_lca IN_LCA: Path to the original (sorted) LCA file (required)
--in_bam IN_BAM: Path to the original (sorted) BAM file (required)
--out_lca OUT_LCA: Path to the short output LCA file (required)
--out_bam OUT_BAM: Path to the short output BAM file (required)
--out_stats OUT_STATS: Path to the output stats file (required)
--out_subs OUT_SUBS: Path to the output subs file (required)
--stranded STRANDED: Either ss for single stranded or ds for double stranded (required)
--mincount MINCOUNT: Minimum read count to keep a node (default: 5)
--k K: Value of k for kmer complexity calculations (default: 5)
--upto UPTO: Keep nodes up to and including this tax threshold, use root to disable (default: family)
--lcaheaderlines LCAHEADERLINES: Number of header lines in input LCA file (default: 0)
--minsim MINSIM: Minimum similarity to reference to keep a read. This should be the same as your ngslca minimum similarity! (default: 0.95)
--keep_keywords KEEP_KEYWORDS [KEEP_KEYWORDS ...]: Other keyword(s) in LCA file for filtering to keep, e.g. Eukaryota (default: none)
--exclude_keywords EXCLUDE_KEYWORDS [KEEP_KEYWORDS ...]: Other keyword(s) in LCA file for filtering to delete (default Hominidae) - NOT IMPLEMENTED YET
--exclude_keyword_file EXCLUDE_KEYWORD_FILE: Text file containing list of keywords or tax paths to remove in LCA file, one on each line - NOT IMPLEMENTED YET
```

Example usage:
```sh
python BamDam.py \
    --in_lca "S32_minsim95.lca" \
    --in_bam "S32.sorted.bam" \
    --out_lca "S32_shortened.lca" \
    --out_bam "S32_shortened.bam" \
    --out_stats "S32_allstats.txt" \
    --out_subs "S32_allsubs.txt" \
    --stranded ds \
    --mincount 5 \
    --k 5 \
    --upto "family" \
    --lcaheaderlines 0 \
    --minsim 0.95 \
    --keep_keywords Eukaryota
```

## Explanation of the output stats file columns
- **TaxNodeID**: The tax node ID from the lca file.
- **"TaxName"**: The tax name rom the lca file.
- **"TotalReads"**: The number of reads assigned to that node or underneath.
- **"PMDsover2"**: The number of read assigned that node or underneath with mean PMD scores over 2 (where the mean is over the alignments of that read, since a PMD score is computed per alignment).
- **"PMDSover4"**: The number of read assigned that node or underneath with mean PMD scores over 4.
- **"Damaged+1"**: The number of reads assigned to that node or underneath where every alignment of that read had a C->T on the 5' (+1) position.
- **"Damaged-1"**: The number of reads assigned to that node or underneath where every alignment of that read had a C->T if single stranded, or a G->A if double stranded, on the 3' (-1) position.
- **"TotalAlignments"**: Sum of the number of alignments for all the reads assigned to that node or underneath.
- **"MeanLength"**: The mean length of the reads assigned to that node or underneath.
- **"ANI"**: Average nucleotide identity of the reads assigned to that node or underneath. 
- **"PerReadKmerGI"**: Mean k-mer based Gini index, where the mean is over the reads assigned to that node or underneath. This is like a per-read complexity metric (are there many of the same k-mers in each read?). Ranges from 0 to 1. Lower numbers are better.
- **"ReadSetKmerGI"**: Read-set k-mer based Gini Index. This is like a per-read-set complexity metric (do the k-mers in all of these reads look like each other?). Ranges from 0 to 1. Lower numbers are better.
- **"AvgGC"**: Average GC content of the reads assigned to that node or underneath. 
- **"taxpath"**: The taxonomic path from the lca file.

In all cases where it is unclear, each read (not each alignment) is weighted equally.

You will notice your new shorter bam file is also annotated with PMD scores. A PMD score is calculated per alignment, and it tells you the log likelihood ratio of your read alignment being ancient over not. So, PMD=2 is like "this read is double as likely to be ancient than it is modern". PMD scores are from [Skoglund et al. 2014](https://doi.org/10.1073/pnas.131893411) - have a look at the first figure to understand them a bit better. The parameters of the PMD score calculation are in the Python script, not user changeable other than editing the code, and I haven't played around with them yet. So it might change a bit at some point, but I think they are sufficiently reliable already for evaluating damage when combined with the Damaged+1 and Damaged-1 columns.

## Plotting script

If you want to look at the damage smiley plot of a specific taxa, there is also an R script PlotDamage.R that works from the command line. It reads in the subs file written by BamDam.py, and needs to know a tax id (ideally) or a tax name (less ideally). You can find the relevant tax id by grep-ing the tax name in the lca file. 

Usage: 
```sh
./PlotDamage.R [options]

Options:
        -f CHARACTER, --subs_file=CHARACTER
                Path to the substitutions file (required)

        -t CHARACTER, --tax=CHARACTER
                Tax name or id. Accepts (unique) tax names in theory, but it's much safer to use tax ids. (required)

        -s CHARACTER, --stranded=CHARACTER
                ds for double stranded or ss for single stranded. (required)

        -o CHARACTER, --output_plot=CHARACTER
                Output file for the plot ending in .png or .pdf (required)

        -h, --help
                Show this help message and exit
```

Example plot (this one has a lot of reads contributing to it, they will generally be noisier):
![Damage plot of Dryas](https://github.com/bdesanctis/bamdam/blob/main/dryas.png)

## Coming soon

Aka my to do list :). Please let me know if you have requests, or really like any of these ideas.

- write a plotting function in R that runs on the lca + bam to plot pmd distribution
- write another plotting function for read length distribution
- [priority] write a python function, separately from this, to quickly extract reads from a bam or fq file assigned to a specific genus (see below) 
- consider writing a function to add some metadmg parameters to the output stats file 
- consider writing an option to NOT include nodes under ones with mincount if those underneath nodes don't meet mincount themselves. (you could just pass an option keepunder = False or something like this)
- implement progress bars w tqdm
- [priority] write a function to combine multiple stats and subs files from multiple different files, like from ellesmere . then you can combine the "read-pulling-out function" with this and make a full-on damage distribution for different taxa in ellesmere! 
- consider outputting a metadmg-like Dfit once you understand what the hell it is doing 
- add a logger (maybe this is more pain than it's worth)
- [priority] implement remove keywords from lca like Homonidae , maybe take in a whole txt file list of keywords from controls lcas?
- [priority] clean up the keyword handling in general
- consider "damtools" as a name - other suggestions also welcome haha
- [priority] consider rewriting bamfilter evenness of coverage calculation - this would go before ngslca
- rory: are more strands ending in purines (A vs G) than anything else? could check this! https://www.pnas.org/doi/abs/10.1073/pnas.0704665104 
- output ani without the damage substitutions . could make it comparable b/w ss and ds by cutting the read in half 
- eventually - null damage distribution?

## License
This project is licensed under the MIT License - see the LICENSE file for details.



BamDam was written by Bianca De Sanctis in 2024. For assistance please contact bddesanctis@gmail.com.

