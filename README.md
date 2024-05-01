# BCREval
BCREval is used to evaluate bisulfite conversion ratio in WGSBS experiment

Whole genome shotgun bisulfite sequencing (WGSBS) also known as BS-seq has been widely used to measure the methylation of whole genome at single-base resolution. One of the key steps in the assay is converting unmethylated cytosines into thymines (BS conversion). Incomplete conversion of unmethylated cytosines can introduces false positive methylation call. Developing a quick method to evaluate bisulfite conversion ratio (BCR) is benefit for both quality control and data analysis of WGSBS. Here BCReval is a small python script to estimate the BCR by using telomeric repetitive DNA as native spike-in control.

## References
Zhou J, Zhao M, Sun Z, Wu F, Liu Y, Liu X, He Z, He Q, He Q. BCREval: a computational method to estimate the bisulfite conversion ratio in WGBS. BMC Bioinformatics. 2020 Jan 31;21(1):38. doi: 10.1186/s12859-019-3334-z. PMID: 32005131; PMCID: PMC6995172.

## Usage information

usage: BCReval.py [-h] [-s {0,1}] [-n MIN_LENGTH] filename

Estimates bisulfite conversion rate from telomere reads.

positional arguments:
  filename              The input_file should be a FASTQ file.

options:
  -h, --help            show this help message and exit
  -s {0,1}, --strand {0,1}
                        Whether to run in stranded mode (options are 0 and 1)
  -n MIN_LENGTH, --min-length MIN_LENGTH
                        Minimum telomeric repeats (default is 5)


## Output format
Unlike the original BCReval, this script will only print the BCR to stdout with no other metrics or data.
