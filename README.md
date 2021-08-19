# Random Fragment Generation

A critical part of any analysis of integration site preferences/biases is a comparison to random integration. Perhaps the most unbiased way to generate this random dataset is to produce it <i>in silico</i>. Below I will describe the pipeline that I have adapted to generate random fragments of a given genome. This pipeline outputs paired-end files (either fasta or fastq) containing genome fragments up to 150 bp in length. The generated files can then be aligned to the human genome using any user-preferred aligner.

All scripts referred to in this overview can be found on [GitHub](https://github.com/gbedwell/random_fragment_generation).

## Pipeline Overview

In the current pipeline, the only required starting file is a multi-fasta file containing entries for each of the chromosomes in the genome of interest. Scripts are included to operate on this file and:

1. Split the multi-fasta file into multiple single fasta files (one file for each chromosome).
2. Generate the required <code>base_counts</code> file.
3. Generate random fragments from the genome of interest given the outputs of (1) and (2).



### Splitting the Multi-Fasta File

The script provided to split the genome multi-fasta file into individual chromosome files is <code>fasta_split.sh</code>. This file takes three inputs, simply denoted as <code>-a</code>, <code>-b</code>, and <code>-c</code>.

+ Parameter <code>-a</code> defines the path to the genome multi-fasta file.

+ Parameter <code>-b</code> defines the name of a temporary file that converts the multi-fasta file to a two-column TSV file containing header information and sequence. This file will be removed before the script is finished running. <b>Note:</b> if the multi-fasta file being used lists individual chromosome sequences on multiple lines, the sequences are converted to a single line before being put in TSV format.

+ Parameter <code>-c</code> defines the path to the desired output location for the newly generated single-fasta files.



The script can be run in the following way:

```bash
./fasta_split.sh \
-a /path/to/multi-fasta \
-b single_line.txt \
-c /path/to/output/directory
```



### Generating a <code>base_counts</code> file

In addition to individual fasta files for each chromosome, fragment generation depends on a <code>base_counts</code> file. This file is a simple TSV file containing:

1. Each chromosome header.
2. The length of each chromosome.
3. The global start position of each chromosome.
4. The global end position of each chromosome.



In this file, chromosome order doesn't matter. All that matters is that the first entry has a global start position of 0 and the last entry has a global end position equal to the total sum of chromosome lengths. The purpose of this file is to define the chromosome and sequence that corresponds to a given random location in the human genome.



To generate this file from the starting multi-fasta file, run the script <code>base_counts_generation.sh</code>. <b>Note:</b> in order to run this script, the additional script <code>base_counts.R</code> must be in the same directory. <code>base_counts.R</code> requires [Tidyverse](https://www.tidyverse.org/).



<code>base_counts_generation.sh</code> requires two input terms.

+ <code>-a</code> defines the path to the genome multi-fasta file.
+ <code>-b</code> defines the suffix to append to the output files. This is useful for explicitly naming the <code>base_counts</code> file after specific genome builds.

The script can be run as follows:

```bash
./base_counts_generation.sh \
-a /path/to/multi-fasta \
-b CHM13_v1.1
```

When run, the script will print each chromosome sequence on a single line, convert the multi-fasta file a two-column TSV format, and count the length of each chromosome sequence. Then, the generated file will be fed to <code>base_counts.R</code> where the basic arithmetic required to calculate the global starts and global ends of each chromosome is performed. The final output is a tab-delimited text file named <code>base_counts_<SUFFIX>.txt</code> that is ready to be fed into the final script for fragment generation.



### Fragment Generation

The final step in generating random genome fragments is running the script <code>fragment_generation.py</code>. I have adapted this script from earlier versions written by Peter Cherepanov (peter.cherepanov@crick.ac.uk), now a [professor](https://www.crick.ac.uk/research/labs/peter-cherepanov) at The Francis Crick Institute.



In short, my version of the  script is written to:

1. Generate a random number and strand orientation corresponding to a random location in the human genome.
2. Parse the <code>base_counts</code> file to define the chromosome and relative location of the random site.
3. Parse the appropriate chromosome fasta file to obtain the sequence of the chromosome upstream of the random site.
4. "Fragment" the larger upstream fragment according to either specific restriction enzyme combinations or randomly (e.g. sonication or fragmentase).
   + Currently the following restriction enzyme combinations are explicitly defined: MseI, MseI/BglII, NheI/AvrII/SpeI/BamHI.
   + If using random fragmentation, the mean fragment length and standard deviation are user-defined. Fragments are then generated according to a hypothetical normal distribution.
5. Generate paired-end reads up to 150 bp in length corresponding to the 5' and 3' ends of the fragment.
6. Output either fasta or fastq files (user's preference) containing the paired-end reads. If fastq format is desired, dummy quality scores are generated with deteriorating quality approaching the 3' end of the read.



There are several required flags to define when calling <code>fragment_generation.py</code>.

The required flags are:

+ <code>-build</code>: The genome build being used. Currently human genome builds hg19, hg38, and CHM13_v1.1 are implemented.
+ <code>-path</code>: The path to the individual chromosome fasta files.
+ <code>-bc</code>: The full name of the <code>base_counts</code> file.
+ <code>-N</code>: The number of random sites to generate.
+ <code>-frag</code>: The desired fragmentation method.
+ <code>-dir</code>: The desired location of the output files.
+ <code>-distance</code>: The maximum distance allowed to the upstream cleavage site.
+ <code>-file_num</code>: The output file number. This is useful for generating multiple random fragment files in parallel.
+ <code>-output</code>: Defines the output file type (fasta or fastq).



If using random fragmentation, the following parameters should also be defined:

+ <code>-mu</code>: The median fragment size for random fragmentation.
+ <code>-sigma</code>: The standard deviation of fragment size for random fragmentation.



The script can be run as follows:

```bash
python2 fragment_generation.py \
-build CHM13_v1.1 \
-path /path/to/chromosome/files \
-bc base_counts_CHM13_v1.1.txt \
-N 100000 \
-frag random \
-dir /path/to/output/directory \
-distance 50000 \
-file_num 1 \
-output fasta \
-mu 500 \
-sigma 250
```


Descriptions of the script and the parameters themselves can be accessed as follows:

```
python2 fragment_generation.py --help
```
