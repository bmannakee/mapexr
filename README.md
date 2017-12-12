# Title
**mapexr** - a post-variant-calling filter based on **BLASTN**
## Description
This package is an implementation of the **MAPEX** algorithm described in `Mannakee et al. ...`.
The algorithm is a post-variant-calling filter for variants that uses **BLASTN** to determine the fraction of reads supporting a variant call that align best to the site of the variant.
It is intended to be used to filter spurious variant calls resulting from mouse DNA contamination in patient-derived xenograft(PDX) tumors, and to evaluate supporting reads for potential paralog risk in primary tumors.

## Dependencies and Installation
**mapexr** depends on several bioconductor packages for bam file manipulation as well as functions included in the *tidyverse*.
The dependencies can be installed prior to installation of **mapexr** using the following commands. 
```
### Not run
## install dependencies
install.packages("tidyverse")
install.packages("devtools")
source("https://bioconductor.org/biocLite.R")
biocLite(c("Rsamtools","GenomicRanges","GenomicAlignments"))

## install mapexr
devtools::install_github("bmannakee/mapexr",build_vignettes=TRUE)
```

## BLASTN
The package requires a local installation of **BLASTN** and a blast database.
The **BLASTN** installation instructions can be found at <https://www.ncbi.nlm.nih.gov/books/NBK279671/>.
When used to evaluate PDXs **mapexr** requires a blast database constructed from concatentated human and mouse reference fasta files.
It is important that the human reference genome used match the reference genome used to do the initial alignment of reads prior to variant calling.
A human GRCh37 reference file can be downloaded from <http://ftp.1000genomes.ebi.ac.uk/vol1/ftp/technical/reference/>, and the mouse mm9 reference genome can be downloaded from <http://hgdownload.cse.ucsc.edu/downloads.html#mouse>.
Once the genomes are downloaded, the following shell code demonstrates how to prepare a blast database from a concatenated human and mouse reference.
```
### Not run
## Create the combined reference file, 
## replace chr with mchr in the mouse reference, 
## and concatenate to the human reference in the combined reference file
cp /path/to/human_g1k_v37.fasta /path/to/combined_g1k_v37_mm9.fasta
cat /path/to/mouse_mm9.fa | sed 's/chr/mchr/' >> /path/to/combined_g1k_v37_mm9.fasta

## Blast+ needs to be installed
makeblastdb -in /path/to/combined_g1k_v37_mm9.fasta -parse_seqids -dbtype nucl
```
When using **mapexr** to look for potential paralogs in primary tumors, a combined human and mouse reference genome can be used, or a database constructed from just the relevant human reference.
In general any human and mouse reference file can be used, as long as "m" is added to the beginning of the chromosome id in the mouse sequence.
General instructions for using `makeblastdb` can be found at <https://www.ncbi.nlm.nih.gov/books/NBK279688/>.

## Usage
The algorithm takes as input a tumor alignment `.bam` file (with index `.bai`), and a set of variant calls derived from that alignment.
Four types of variant call file are supported, with built in parsers for `VCF >=ver4.0`, the **MuTect1** output format `call_stats.txt`, the **Varscan2** output format `.snp` and the MAFLITE format.
The MAFLITE format is a 5-column tab-separated format with columns `chr`, `start`, `end`, `ref_allele`, `alt_allele`.
A header is not required, but if provided should match these column names for use with **mapexr**.
The package exports one main function `run_mapex()` (See function documentation) as well as the `load_variants` (see function documentation) function which is exported to facilitate trouble-shooting of the parsing of variant file formats.
Example code for `run_mapex()`
```
### not run
## file paths and variables
bam <- "/path/to/tumor.bam"
bamidx <- "/path/to/tumor.bai"
variants <- "/path/to/variants.vcf"
blastout <- "/path/to/blastoutput.txt" # Read level blast results, not required
blastpath <- "/path/to/blastn" # if blast is in the users path, just "blastn" here
blastdb <- "/path/to/combined_db"
threads <- 1 # number of threads consumed by blastn
mapq <- 1 # this is the default for minimum mapq score

results <- run_mapex(bam_file=bam,
		     bam_idx=bamidx,
		     variant_file=variants,
		     variant_file_type='vcf',
		     blast_out_file=blastout,
		     blastn_path=blastpath,
		     blast_threads=threads,
		     min_mapq=mapq)
```
See the `Example` vignette included with the package for detailed description of the results from `run_mapex()` and how they can be used.
