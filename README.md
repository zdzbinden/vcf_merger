VCF File Processing and Merging Script
This script is designed to process and merge two VCF (Variant Call Format) files that have been aligned using the same reference genome. The purpose of this script is to combine VCF files with different sets of individuals (i.e., non-overlapping samples, but similar loci).

The script filters the VCF files to retain only the lines that contain contig names present in a provided .fai file. The .fai file is an index file for the reference genome fasta file, and it contains the names, sizes, and other information about the contigs in the reference genome. Each line of the .fai file represents one contig and is formatted as follows: <contig name>\t<size>\t<offset>\t<linebases>\t<linewidth>.

After filtering, the script replaces the headers of the VCF files with combined headers, sorts, compresses, and indexes the VCF files, and finally merges the VCF files. The script checks for errors at each step and stops if any command fails. It also checks that the required programs are installed before running.

Requirements
The script requires the following programs to be installed:

bcftools
bgzip
tabix
These can be installed using conda or other package managers.

Usage
The script takes as input two VCF files and a .fai file. The .fai file should contain the contig names that are present in the VCF files. The script can be run from the command line as follows:

bash
Copy code
./script.sh file1.vcf file2.vcf file.fai
Replace file1.vcf, file2.vcf, and file.fai with the names of your actual files.

The script outputs a merged VCF file named merged_file.vcf.

Error Checking
The script checks for errors at each step and stops if any command fails. It reports an error message indicating what went wrong. It also checks that the required programs are installed before running.

At the end of the script, it checks that the number of individuals in the merged file is as expected (the sum of the individuals in the two input VCF files). If the number of individuals in the merged file is not as expected, the script reports an error and exits with a non-zero status.

Temporary Files
The script creates several temporary files during its execution. These files are removed as soon as they are no longer needed. The original input files are not modified.
