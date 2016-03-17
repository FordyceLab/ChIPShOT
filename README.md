![ChIPShOT Logo](./logo.png "Gunga galunga!")

# The Chromatin Immunoprecipitation (and Sequencing) Shape/Occupancy Toolset

Recent developments in the study of transcription factor/DNA interactions have shown that, in addition to primary DNA sequence recognition, the local shape of DNA at and around the binding site can influence the the strength of the interaction [1, 2]. ChIPShOT is a Python 3 package and command line tool that makes it easy to download yeast ChIP-seq data from Gene Expression Omnibus, filter and align these data, and analyze the resulting binding sites and flanking regions for patterns in DNA shape that correlate strongly with occupancy.

## Installation

ChIPShOT relies on several external tools to do much of the heavy lifting. The external dependency list includes:

+ [BWA](http://bio-bwa.sourceforge.net)
+ [SRA Toolkit](http://www.ncbi.nlm.nih.gov/Traces/sra/sra.cgi?view=software)
+ [SAMtools](http://samtools.sourceforge.net)
+ [MACS2](https://github.com/taoliu/MACS)

These tools can be downloaded and installed using the instructions found on their respective websites. Alternatively, BWA, SRA Toolkit, and SAMtools are all available through [Homebrew](http://brew.sh/) (in the `science` tap) and MACS2 is available to be installed through pip (Python 2 only).

Some of ChIPShOT's Python dependecies can be difficult to install automatically. Therefore, it is suggested that users install these dependencies manually using `pip3`, as below:

```
pip3 install numpy
pip3 install scipy
pip3 install biopython
```

The remaining Python dependencies are installed automatically when the ChIPShOT package itself is installed.

ChIPShOT installation can be accomplished by:
1. Downloading or cloning this repository onto your local machine
2. Changing directory into `ChIPShOT` and running `pip3 install .` (don't forget the ".")

## Usage

ChIPShOT provides a command line tool (accessed as `chipshot`) with four useful subcommands to make basic ChIP-seq anaysis and peak recentering easy. Each subcommand has its own help statement that can be accessed with the command `chipshot <subcommand> -h`, filling in your subcommand of interest.

### `download`

This command takes as input an ID for a sample in GEO and a prefix for the output files. The sample ID can be taken diretly from the GEO accession viewer for your experiment of interest. For [this](http://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE29506) particular experiment, the sample IDs can be found under the "Samples" subsection and are numbered GSM730517-GSM730535. This command also requires an email adress (required by GEO to search their database programmatically). This can be specified with the `--email` argument. Two optional arguments, `--fastq` and `--sam`, can be passed to tell ChIPShOT to dump the respective file type from the Sequence Read Archive file (`.sra`).

### `cleanref`

This command takes as input the reference genome fasta file and an output prefix. It is suggested that you use the latest yeast genome release, which can be downloaded [here](http://downloads.yeastgenome.org/sequence/S288C_reference/genome_releases/). The `.fsa` file from this archive will be used for all alignments and peak recentering. However, the chromosome names from this version of the FASTA file are not easily parsed. Running the reference genome through the `cleanref` subcommand will rename the chromosomes so that they follow the yeast convention. 

### `align`

This command filters and trims the raw fastq reads and aligns the reads to the genome. Users can specify trimming parameters using the `--five`, `--three` or `--both` parameters, which specify the minimum quality score to trim to from the five prime end, three prime end, or both ends of each read, respectively. Using the `--quality` parameter, the user can specify a minimum mean quality score, below which the read will be filtered out before alignment. Similarly, the `--minlength` parameter, the user can specify a minimum length for each read after trimming is complete. The user must supply a reference genome with the `--reference` parameter in order for alignment to proceed after read trimming and filtering is complete. Additionally, the user must supply the path to fastq file as well as a prefix for the output file.

### `callpeaks`

This command calls and recenters peaks based on read depth at each locus. The command takes as input a bam file and a prefix for the output files. Users can optionally input an associated control bam file using the `--control` parameter. Additionally, users must specify a the reference genome with the `--reference` parameter, a position frequency matrix from [JASPAR](http://jaspar.genereg.net) with the `--pfm` parameter, and the distance to extend the centered peak on eaither side of the end of the putative binding site with the `--extend` parameter.

### Bonus: `tellmeastory`

ChIPShOT will even tell you a story if you ask nicely.

## References:

1. Abe, N., Dror, I., Yang, L., Slattery, M., Zhou, T., Bussemaker, H. J., et al. (2015). Deconvolving the recognition of DNA shape from sequence. Cell, 161(2), 307–318. http://doi.org/10.1016/j.cell.2015.02.008

2. Zhou, T., Shen, N., Yang, L., Abe, N., Horton, J., Mann, R. S., et al. (2015). Quantitative modeling of transcription factor binding specificities using DNA shape. Proceedings of the National Academy of Sciences, 112(15), 4654–4659. http://doi.org/10.1073/pnas.1422023112
