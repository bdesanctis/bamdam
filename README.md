# bamdam: A post-mapping toolkit for ancient environmental DNA capture or shotgun sequencing data

## Table of Contents
- [Description](#description)
- [Installation](#installation)
- [Main Script](#main-script)
- [Accessory scripts](#damage-plotting-script)
- [Coming Soon](#coming-soon)

## Description

WARNING - this software is still in development. Things are changing often here!

Please contact me with any issues, unexpected behavior, or requests. bddesanctis@gmail.com

Welcome to bamdam! The goal of this toolkit is to provide functionality after capture or shotgun sequencing aeDNA reads have been mapped to a reference database and run through a least common ancestor algorithm like [ngsLCA](https://github.com/miwipe/ngsLCA) using the reference database taxonomy to obtain an ngsLCA-style .lca file. If you have a lot of data, long term storing these often giant bams post-mapping can be annoying, so the first point of bamdam is to write (much) smaller versions that still include all the relevant information for the project. The second point of bamdam is to compute a ton of read set metrics to determine if a taxonomic node looks like a real ancient taxa, rather than a modern environmental, lab, or database contaminant. 

There are also the following accessory scripts:
- PlotDamage.R for smiley plots,
- extractreads.py for quickly extracting reads from a bam belonging to a specific tax id(s), and
- PlotPMD.R for plotting the PMD score distribution of a taxa. 

Bamdam was heavily inspired by [metaDMG](https://github.com/metaDMG-dev/metaDMG-cpp) and [filterBAM](https://github.com/genomewalker/bam-filter). It is not particularly optimized for speed, and doesn't thread yet. On the other hand, it doesn't ever read a whole file into memory, so it shouldn't need that much RAM (please tell me if you manage to crash it and how). It can subset a 40GB bam to family level in an hour on my laptop and this should scale close to linearly, which feels sufficient for now.

## Installation

Hopefully pretty straightforward. It's got a couple python and R packages as dependencies. 

```sh
git clone https://github.com/bdesanctis/bamdam.git
cd bamdam
chmod +x PlotDamage.R
chmod +x PlotPMD.R
```

## Main Script

Input: Bam and LCA files. Output: Shortened bam, shortened lca, stats and subs files.

The main script BamDam.py will:
1. Shorten the LCA file to only nodes which are equal to or below your tax threshold, which meet your minimum read count, or are below nodes which meet your minimum read count. You may also give it a list or file of tax identifiers and it will exclude all reads assigned to that node, and optionally all reads assigned to nodes below that node as well (--exclude_under). For exclusions, you can give it tax IDs, full tax names, or full tax entries; e.g. Homonidae, "Homo sapiens", 4919, etc, but ideally use full tax paths like "4919:Homo sapiens:species".
2. Shorten the bam file to include only reads which appear in the newly shortened LCA file. This can make 100-1000x or even more difference in the size of a bam file after typical ancient eDNA mapping steps. It will also annotate the new bam file with PMD scores as in PMDTools (in the DS:Z field).
3. Gather various statistics per node and write to a stats file, sorted by total reads assigned. More information detailing this file is below.
4. Gather substitution patterns per node and write to a subs file, sorted in the same order as the stats file (so line n of the stats file should correspond to line n of the subs file). 

You may also plot damage plots for individual tax IDs by inputting the subs file into an R script, which is detailed in the next section. 

Please ensure your input bam is sorted in the same order as your input lca file, and that every read in your LCA file is also in your bam file, both of which will be true by default if you use fresh ngsLCA output. On the other hand, bamdam will be fine if there are reads in your bam that are not in your LCA file, which means you can pre-subset your LCA file if you'd like. Pre-subsetting your lca file may give a speed-up and will end you up with an even smaller bam, which you should consider doing if you only wanted to keep eukaryotes, for example (e.g.: grep "Eukaryot" in.lca > out.lca). Also, for now, please provide a ngsLCA-formatted lca file, not a metaDMG-formatted lca file (they have slightly different formats).

Usage and options: 
```sh
usage: BamDam.py [-h] --in_lca IN_LCA --in_bam IN_BAM --out_lca OUT_LCA --out_bam OUT_BAM --out_stats OUT_STATS --out_subs OUT_SUBS --stranded STRANDED
                 [--mincount MINCOUNT] [--k K] [--upto UPTO] [--minsim MINSIM] [--exclude_keywords EXCLUDE_KEYWORDS [EXCLUDE_KEYWORDS ...]] [--exclude_keyword_file EXCLUDE_KEYWORD_FILE]

options:
  -h, --help            show this help message and exit
  --in_lca              Path to the original (sorted) LCA fil (required)
  --in_bam              Path to the original (sorted) BAM file (required)
  --out_lca             Path to the short output LCA file (required)
  --out_bam             Path to the short output BAM file (required)
  --out_stats           Path to the output stats file (required)
  --out_subs            Path to the output subs file (required)
  --stranded            Either ss for single stranded or ds for double stranded (required)
  --mincount            Minimum read count to keep a node (default: 5)
  --k                   Value of k for kmer complexity calculations (default: 5)
  --upto                Keep nodes up to and including this tax threshold, use root to disable (default: family)
  --minsim              Minimum similarity to reference to keep a read; must match ngslca min similarity (default: 0.95)
  --exclude_keywords    Keyword(s) to exclude when filtering (default: none)
  --exclude_keyword_file  File of keywords to exclude when filtering, one per line (default: none)
  --exclude_under       Set this flag if you also want to exclude all nodes underneath the ones you've specified

```
Example usage:
```sh
python BamDam.py \
    --in_lca "S32_minsim95.lca" \
    --in_bam "S32.sorted.bam" \
    --out_lca "S32_shortened.lca" \
    --out_bam "S32_shortened.bam" \
    --out_stats "S32_stats.txt" \
    --out_subs "S32_subs.txt" \
    --stranded ds \
    --mincount 5 \
    --k 5 \
    --upto "family" \
    --minsim 0.95 \
    --exclude_keywords "Homonidae" "Felidae" 
```

### Explanation of the output stats file columns
- **TaxNodeID**: The tax node ID from the lca file.
- **TaxName**: The tax name from the lca file.
- **TotalReads**: The number of reads assigned to that node or underneath.
- **PMDsover2**: The proportion of reads assigned to that node or underneath with mean PMD scores over 2 (where the mean is over the alignments of that read, since a PMD score is computed per alignment).
- **PMDSover4**: The proportion of reads assigned to that node or underneath with mean PMD scores over 4.
- **ND+1**: Normalized damage +1: The proportion of reads assigned to that node or underneath with a C->T on the 5' (+1) position, minus the mean (non C>T or G>A) divergence for that node.
- **ND-1**: Normalized damage -1: The proportion of reads assigned to that node or underneath with a C->T if single stranded, or a G->A if double stranded, on the 3' (-1) position, minus the mean (non C>T or G>A) divergence for that node.
- **TotalAlignments**: Sum of the number of alignments for all the reads assigned to that node or underneath.
- **MeanLength**: The mean length of the reads assigned to that node or underneath.
- **Div**: The mean divergence for that node, not including any C>T or G>A transitions.
- **ANI**: Average nucleotide identity of the reads assigned to that node or underneath. 
- **PerReadKmerGI**: Mean k-mer based Gini index, where the mean is over the reads assigned to that node or underneath. This is like a per-read complexity metric (are there many of the same k-mers in each read?). Ranges from 0 to 1. Lower numbers are better. If this number is substantially bigger for a node of interest than others, consider looking more carefully at the actual reads mapping to that node.
- **ReadSetKmerGI**: Read-set k-mer based Gini Index. This is like a per-read-set complexity metric (do the k-mers in all of these reads look like each other?). Ranges from 0 to 1. Lower numbers are better. If this number is substantially bigger for a node of interest than others, consider looking more carefully at the actual reads mapping to that node.
- **AvgGC**: Average GC content of the reads assigned to that node or underneath.
- **Damaged+1**: The raw proportion of reads assigned to that node or underneath where every alignment of that read had a C->T on the 5' (+1) position.
- **Damaged-1**: The raw proportion of reads assigned to that node or underneath where every alignment of that read had a C->T if single stranded, or a G->A if double stranded, on the 3' (-1) position.
- **taxpath**: The taxonomic path from the lca file.

In all cases where it is unclear, each read (not each alignment) is weighted equally.

You will notice your new shorter bam file is also annotated with PMD scores. A PMD score is calculated per alignment, and it tells you the log likelihood ratio of your read alignment being ancient over not. So, PMD = 2 is like "this read is double as likely to be ancient than it is modern". PMD scores are from [Skoglund et al. 2014](https://doi.org/10.1073/pnas.131893411) - have a look at the first figure to understand them a bit better. The parameters of the PMD score calculation are in the Python script, not user changeable other than editing the code, and I haven't played around with them yet. So it might change a bit at some point, but I think they are sufficiently reliable already for evaluating damage when combined with the other damage metrics.

## Accessory scripts

### Damage plotting script

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

<img src="https://github.com/bdesanctis/bamdam/blob/main/dryas.png" width="600">

### Script to extract reads

This is just like a fancy grep command. You could write a bash command to do the same thing, but it would probably be slower. This script leverages the fact that the lca and bam files are in the same order. Just like grep, you can use the full tax identifier, e.g. "48230:Dryas:genus", to be totally sure things will work correctly.

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

### PMD distribution plotting script

This is a very basic histogram plotting wrapper. Once you've extracted a bam file of reads associated to a specific node, you can plot its PMD score distribution. In theory you could also use this script to plot PMD scores for the entire bam if you'd like, but it would be slow and might crash if the bam is too big. The following will extract all PMD scores from a bam file, feed them into a basic R histogram plotting script, and save the histogram to "salixPMDs.png" with the title "Salix PMD scores":

```sh
samtools view onlysalix.bam | grep -o 'DS:Z:[^ ]*' | sed 's/DS:Z://' | ./PlotPMD.R -o salixpmds.png -t "Salix"
```

## Coming soon
Aka my to do list, suggestions welcome. I will probably do most if not all of the below in the next two months:
- [priority] make it so the new bam also has a shorter header, right now it is just copying the old header and that is mildly ridiculous
- evenness of coverage, probably immediately after shortening the bam + lca
- split functions to shorten and compute stats into two i think
- write another simple plotting function for read length distribution for a specific taxa
- are more strands ending in purines (A vs G) than anything else? could check this eg https://www.pnas.org/doi/abs/10.1073/pnas.0704665104 
- derive a kmer distribution from the modern reads, and remove reads with many of those kmers from ancient samples, maybe once you're already on the bams so you can see what you're removing, could do this simultaneously with evenness of coverage
- eventually - null damage distribution?
- add some example data
- implement some functionality for combining samples , eg. multi-sample damage plots
- clarify kmer metrics
- clarify y axis on damage plots is prop of reads not prop of subs
- make minsim params match ngslca ones , -editdistmin -editdistmax -simscorelow and -simscorehigh
- make extractreads.py work on fqs too
- larger-picture workflow including remove unnecessary headers and a bwa modification/little script for multimapping (maybe a separate repo for a bwa fork though)
- add a root line to the end of the subs file for quick whole-library damage assessment
- normalized c>t on +1, -1 minus c>t in the middle eg [+10,-10]
- conditional damage metric? prop of reads with c>t on 5' that also have g>a on 3'. this is a leipzig thing
- write a function to combine multi sample subs+stats files
- multi sample damage plots

## License
This project is licensed under the MIT License - see the LICENSE file for details.

BamDam was written by Bianca De Sanctis in 2024. For assistance please contact bddesanctis@gmail.com.
