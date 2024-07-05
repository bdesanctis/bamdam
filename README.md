# bamdam: A post-mapping toolkit for ancient environmental DNA capture or shotgun sequencing data

## Important note - this software is still in development!

Please contact me with any issues, unexpected behavior, or requests! bddesanctis@gmail.com

Welcome to bamdam! The goal of this toolkit is to provide functionality after shotgun sequencing aeDNA reads have been mapped to a reference database and run through ngsLCA using the reference database taxonomy to obtain a .lca file. If you have a lot of data, long term storing these often giant bams post-mapping can be annoying, so the first point of bamdam is to write (much) smaller versions that still include all the relevant information for the project. The second point of bamdam is to compute a ton of read set metrics to determine if a taxonomic node looks like a real ancient taxa. There are also associated extra scripts for plotting functionality and quickly extracting reads from a fq or bam belonging to a specific tax id(s). More to come, probably.

This pipeline and software is not appropriate for metabarcoding data. 

Bamdam was heavily inspired by [metaDMG](https://github.com/metaDMG-dev/metaDMG-cpp) and [filterBAM](https://github.com/genomewalker/bam-filter). It is not particularly optimized for speed, and doesn't thread yet. On the other hand, it doesn't ever read a whole file into memory, so it shouldn't need that much RAM (please tell me if you manage to crash it and how). It can subset a 40GB bam to family level in an hour on my laptop and this should scale close to linearly. That feels sufficient to me for now - let me know if it's not for you. I could implement threading.

## Table of Contents
- [Installation](#installation)
- [BamDam.py: Main Script](#main-script)
- [BamDam.py: Explanation of the Output Stats File Columns](#explanation-of-the-output-stats-file-columns)
- [PlotDamage.R: Plotting Script](#plotting-script)
- [extractreads.py: Script to extract reads](#script-to-extract-reads)
- [Steps before bamdam](#steps-before-bamdam)
- [Coming Soon](#coming-soon)
- [License](#license)

## Installation

Hopefully pretty straightforward. 

```sh
git clone https://github.com/bdesanctis/bamdam.git
cd bamdam
chmod +x PlotDamage.R
```

Should complain when you try to run scripts if you don't have the appropriate python or R libraries installed, in which case install then run again.

## Main Script

Input: Bam and LCA files. Output: Shortened bam, shortened lca, stats and subs files.

The main script BamDam.py will:
1. Shorten the LCA file to only nodes which are equal to or below your tax threshold, meet your minimum read count (or are below nodes which meet your minimum read count). You may also give it a list of keywords or phrases and it will exclude all lines including any of those keywords or phrases (e.g. Homonidae). 
2. Shorten the bam file to include only reads which appear in the newly shortened LCA file. This can make 100-1000x or even more difference in the size of a bam file after typical ancient eDNA mapping steps. It will also annotate the new bam file with PMD scores as in PMDTools (in the DS:Z field).
3. Gather various statistics per node and write to a stats file, sorted by total reads assigned. More information detailing this file is below.
4. Gather substitution patterns per node and write to a subs file, sorted in the same order as the stats file (so line n of the stats file should correspond to line n of the subs file). 

You can then look through your stats file and decide what seems real and ancient to you. You may also plot damage plots for individual tax IDs using an R script and the subs file which is detailed in the next section. The subs file should be fairly self-explanatory in formatting, if you are interested in looking at it.

Please ensure your input bam is sorted in the same order as your input lca file, and that every read in your LCA file is also in your bam file, both of which will be true by default if you use fresh ngsLCA output. On the other hand, bamdam will be fine if there are reads in your bam that are not in your LCA file, which means you can pre-subset your LCA file if you'd like. Pre-subsetting your lca file may give a speed-up and will end you up with an even smaller bam, which you should consider doing if you only wanted to keep eukaryotes, for example (e.g.: grep "Eukaryot" in.lca > out.lca). Also, for now, please provide a ngsLCA-formatted lca file, not a metaDMG-formatted lca file (they have slightly different formats).

Usage and options: 
```sh
BamDam.py [-h] --in_lca IN_LCA --in_bam IN_BAM --out_lca OUT_LCA --out_bam OUT_BAM --out_stats OUT_STATS --out_subs OUT_SUBS --stranded STRANDED
          [--mincount MINCOUNT] [--k K] [--upto UPTO] [--lcaheaderlines LCAHEADERLINES] [--minsim MINSIM]
          [--exclude_keywords KEEP_KEYWORDS [KEEP_KEYWORDS ...]]
          
options:

-h, --help : show this help message and exit
--in_lca : Path to the original (sorted) LCA file (required)
--in_bam : Path to the original (sorted) BAM file (required)
--out_lca : Path to the short output LCA file (required)
--out_bam : Path to the short output BAM file (required)
--out_stats : Path to the output stats file (required)
--out_subs : Path to the output subs file (required)
--stranded : Either ss for single stranded or ds for double stranded (required)
--mincount : Minimum read count to keep a node (default: 5)
--k : Value of k for kmer complexity calculations (default: 5)
--upto : Keep nodes up to and including this tax threshold, use root to disable (default: family)
--lcaheaderlines : Number of header lines in input LCA file - just go look at the file real quick and find out (default: 0)
--minsim : Minimum similarity to reference to keep a read. This needs to be the same as your ngslca minimum similarity. (default: 0.95)
--exclude_keywords : Keyword(s) in LCA file for filtering to exclude lines containing (default: Hominidae)
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
    --exclude_keywords "Homonidae" "Felidae" \
    > S32_bamdam.log
```

## Explanation of the output stats file columns
- **TaxNodeID**: The tax node ID from the lca file.
- **"TaxName"**: The tax name rom the lca file.
- **"TotalReads"**: The number of reads assigned to that node or underneath.
- **"PMDsover2"**: The number of reads assigned to that node or underneath with mean PMD scores over 2 (where the mean is over the alignments of that read, since a PMD score is computed per alignment).
- **"PMDSover4"**: The number of reads assigned to that node or underneath with mean PMD scores over 4.
- **"Damaged+1"**: The number of reads assigned to that node or underneath where every alignment of that read had a C->T on the 5' (+1) position.
- **"Damaged-1"**: The number of reads assigned to that node or underneath where every alignment of that read had a C->T if single stranded, or a G->A if double stranded, on the 3' (-1) position.
- **"TotalAlignments"**: Sum of the number of alignments for all the reads assigned to that node or underneath.
- **"MeanLength"**: The mean length of the reads assigned to that node or underneath.
- **"ANI"**: Average nucleotide identity of the reads assigned to that node or underneath. 
- **"PerReadKmerGI"**: Mean k-mer based Gini index, where the mean is over the reads assigned to that node or underneath. This is like a per-read complexity metric (are there many of the same k-mers in each read?). Ranges from 0 to 1. Lower numbers are better. If this number is substantially bigger for a node of interest than others, consider looking more carefully at the actual reads mapping to that node.
- **"ReadSetKmerGI"**: Read-set k-mer based Gini Index. This is like a per-read-set complexity metric (do the k-mers in all of these reads look like each other?). Ranges from 0 to 1. Lower numbers are better. If this number is substantially bigger for a node of interest than others, consider looking more carefully at the actual reads mapping to that node.
- **"AvgGC"**: Average GC content of the reads assigned to that node or underneath. 
- **"taxpath"**: The taxonomic path from the lca file.

In all cases where it is unclear, each read (not each alignment) is weighted equally.

You will notice your new shorter bam file is also annotated with PMD scores. A PMD score is calculated per alignment, and it tells you the log likelihood ratio of your read alignment being ancient over not. So, PMD = 2 is like "this read is double as likely to be ancient than it is modern". PMD scores are from [Skoglund et al. 2014](https://doi.org/10.1073/pnas.131893411) - have a look at the first figure to understand them a bit better. The parameters of the PMD score calculation are in the Python script, not user changeable other than editing the code, and I haven't played around with them yet. So it might change a bit at some point, but I think they are sufficiently reliable already for evaluating damage when combined with the Damaged+1 and Damaged-1 columns.

## Plotting script

If you want to look at the damage smiley plot of a specific taxa, there is also an R script PlotDamage.R that works from the command line. It reads in the subs file written by BamDam.py, and needs to know a tax id (ideally) or a unique tax name (less ideally). You can find the relevant tax id by grep-ing the tax name in the lca file. 

Usage: 
```sh
./PlotDamage.R [options]

Options:
        -f CHARACTER, --subs_file=CHARACTER
                Path to the substitutions file (required)
        -t CHARACTER, --tax=CHARACTER
                Tax name or id. Accepts (unique) tax names in theory, but it is much safer to use tax ids. (required)
        -s CHARACTER, --stranded=CHARACTER
                ds for double stranded or ss for single stranded. (required)
        -o CHARACTER, --output_plot=CHARACTER
                Output file for the plot ending in .png or .pdf (required)
        -h, --help
                Show this help message and exit
```

Usage example:
```sh
./PlotDamage.R -f subs.txt -t 48230 -s ds -o dryas.png
```

Example double-stranded plot (this one has a lot of reads contributing to it, they will generally be noisier):

![Damage plot of Dryas](https://github.com/bdesanctis/bamdam/blob/main/dryas.png)

## Script to extract reads

This is just like a fancy grep command. You could write a bash command to do the same thing, but it would probably be slower. This script leverages the fact that the lca and bam files are in the same order. Just like grep, you can use the full tax identifier, e.g. "48230:Dryas:genus", to be totally sure things will work correctly (perhaps there is a bacteria with "Dryas" in the name, I don't know).

Usage:
```sh
extractreads.py [-h] --in_lca IN_LCA --in_bam IN_BAM --out_bam OUT_BAM --keyword KEYWORD

options:
  -h, --help         show this help message and exit
  --in_lca IN_LCA    Path to the (sorted) LCA file
  --in_bam IN_BAM    Path to the (sorted) BAM file
  --out_bam OUT_BAM  Path to the filtered BAM file
  --keyword KEYWORD  Keyword or phrase to filter for
```

Usage example:
```sh
python extractreads.py --in_lca in.lca --in_bam in.bam --keyword "Salix" --out_bam onlysalix.bam
```

## Steps before bamdam

To do. Gonna write a paragraph here about mapping to reference databases, where to find taxonomy files, ngslca, prefiltering lca files etc.

## Coming soon

Aka my to do list :). Please let me know if you have requests, or really like any of these ideas.

- [priority!!!] make it so the new bam also has a shorter header, right now it is just copying the old header and that is mildly ridiculous
- [priority!!!] output the difference between other 5' and c>t 5' (same with 3') as a difference in proportions too i think

- [priority] consider rewriting bamfilter evenness of coverage calculation - this would go before ngslca
- write a plotting function in R that runs on the lca + bam to plot pmd distribution
- write another plotting function for read length distribution
- consider writing a function to add some metadmg parameters to the output stats file 
- consider writing an option to NOT include nodes under ones with mincount if those underneath nodes don't meet mincount themselves. (you could just pass an option keepunder = False or something like this)
-  write a function to combine multiple stats and subs files from multiple different files, like from ellesmere . then you can combine the "read-pulling-out function" with this and make a full-on damage distribution for different taxa in ellesmere! 
- add a logger (maybe this is more pain than it's worth)
- consider "damtools" as a name - other suggestions also welcome haha
- rory: are more strands ending in purines (A vs G) than anything else? could check this! https://www.pnas.org/doi/abs/10.1073/pnas.0704665104 
- output ani without the damage substitutions . could make it comparable b/w ss and ds by cutting the read in half
- yucheng: derive a kmer distribution from the modern reads, and remove reads with many of those kmers from ancient samples
- eventually - null damage distribution?
- add some example data

## License
This project is licensed under the MIT License - see the LICENSE file for details.



BamDam was written by Bianca De Sanctis in 2024. For assistance please contact bddesanctis@gmail.com.

