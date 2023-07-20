#!/bin/bash

# This script takes as input two VCF files and a .fai file, filters the VCF files to retain only the lines that contain contig names present in the .fai file, replaces the headers of the VCF files with combined headers, sorts, compresses, and indexes the VCF files, and finally merges the VCF files. It checks for errors at each step and stops if any command fails. It also checks that the required programs are installed before running.

# Usage: ./script.sh file1.vcf file2.vcf file.fai

# Check that the required programs are installed
command -v bcftools >/dev/null 2>&1 || { echo "bcftools is not installed. Aborting."; exit 1; }
command -v bgzip >/dev/null 2>&1 || { echo "bgzip is not installed. Aborting."; exit 1; }
command -v tabix >/dev/null 2>&1 || { echo "tabix is not installed. Aborting."; exit 1; }

# Input files
file1=$1
file2=$2
fai_file=$3

# Extract the contig names from the .fai file and save them to a temporary file
awk '{print $1}' $fai_file > contigs.txt || { echo "Failed to extract contig names from .fai file"; exit 1; }

# Filter the VCF files to retain only the header lines and the lines that contain contig names present in the .fai file
awk 'BEGIN{while(getline < "contigs.txt") contigs[$1]} /^#/ || ($1 in contigs)' $file1 > ${file1%.vcf}_filtered.vcf || { echo "Failed to filter VCF file 1"; exit 1; }
awk 'BEGIN{while(getline < "contigs.txt") contigs[$1]} /^#/ || ($1 in contigs)' $file2 > ${file2%.vcf}_filtered.vcf || { echo "Failed to filter VCF file 2"; exit 1; }

# Extract the original headers from the VCF files, separating the final line
grep "^##" $file1 > ${file1%.vcf}_original.hdr || { echo "Failed to extract original header from VCF file 1"; exit 1; }
grep "^#" $file1 | grep -v "^##" > ${file1%.vcf}_final.hdr || { echo "Failed to extract final header line from VCF file 1"; exit 1; }
grep "^##" $file2 > ${file2%.vcf}_original.hdr || { echo "Failed to extract original header from VCF file 2"; exit 1; }
grep "^#" $file2 | grep -v "^##" > ${file2%.vcf}_final.hdr || { echo "Failed to extract final header line from VCF file 2"; exit 1; }

# Create a header file from the .fai file, only including contigs present in the VCF files
awk 'BEGIN{while(getline < "contigs.txt") contigs[$1]} {if ($1 in contigs) print "##contig=<ID="$1",length="$2">"}' $fai_file > contigs.hdr || { echo "Failed to create header file from .fai file"; exit 1; }

# Append the contig information to the headers of the VCF files
cat contigs.hdr >> ${file1%.vcf}_original.hdr || { echo "Failed to append contig information to header of VCF file 1"; exit 1; }
cat contigs.hdr >> ${file2%.vcf}_original.hdr || { echo "Failed to append contig information to header of VCF file 2"; exit 1; }

# Append the final line of the original header to the new header
cat ${file1%.vcf}_final.hdr >> ${file1%.vcf}_original.hdr || { echo "Failed to append final header line to header of VCF file 1"; exit 1; }
cat ${file2%.vcf}_final.hdr >> ${file2%.vcf}_original.hdr || { echo "Failed to append final header line to header of VCF file 2"; exit 1; }

# Replace the headers of the VCF files with the combined headers
bcftools reheader -h ${file1%.vcf}_original.hdr -o ${file1%.vcf}_reheader.vcf ${file1%.vcf}_filtered.vcf || { echo "Failed to replace header of VCF file 1"; exit 1; }
bcftools reheader -h ${file2%.vcf}_original.hdr -o ${file2%.vcf}_reheader.vcf ${file2%.vcf}_filtered.vcf || { echo "Failed to replace header of VCF file 2"; exit 1; }

# Remove the temporary files
rm contigs.txt contigs.hdr ${file1%.vcf}_final.hdr ${file2%.vcf}_final.hdr ${file1%.vcf}_original.hdr ${file2%.vcf}_original.hdr ${file1%.vcf}_filtered.vcf ${file2%.vcf}_filtered.vcf || { echo "Failed to remove temporary files"; exit 1; }

# Sort the VCF files
bcftools sort ${file1%.vcf}_reheader.vcf -Ov -o ${file1%.vcf}_reheader_sorted.vcf || { echo "Failed to sort VCF file 1"; exit 1; }
bcftools sort ${file2%.vcf}_reheader.vcf -Ov -o ${file2%.vcf}_reheader_sorted.vcf || { echo "Failed to sort VCF file 2"; exit 1; }

# Compress and index the VCF files
bgzip -c ${file1%.vcf}_reheader_sorted.vcf > ${file1%.vcf}_reheader_sorted.vcf.gz || { echo "Failed to compress sorted VCF file 1"; exit 1; }
bgzip -c ${file2%.vcf}_reheader_sorted.vcf > ${file2%.vcf}_reheader_sorted.vcf.gz || { echo "Failed to compress sorted VCF file 2"; exit 1; }
tabix -p vcf ${file1%.vcf}_reheader_sorted.vcf.gz || { echo "Failed to index compressed VCF file 1"; exit 1; }
tabix -p vcf ${file2%.vcf}_reheader_sorted.vcf.gz || { echo "Failed to index compressed VCF file 2"; exit 1; }

# Merge the VCF files
bcftools merge ${file1%.vcf}_reheader_sorted.vcf.gz ${file2%.vcf}_reheader_sorted.vcf.gz -Ov -o merged_file.vcf || { echo "Failed to merge VCF files"; exit 1; }

# Remove the temporary files
rm ${file1%.vcf}_reheader.vcf ${file2%.vcf}_reheader.vcf ${file1%.vcf}_reheader_sorted.vcf ${file2%.vcf}_reheader_sorted.vcf ${file1%.vcf}_reheader_sorted.vcf.gz ${file2%.vcf}_reheader_sorted.vcf.gz ${file1%.vcf}_reheader_sorted.vcf.gz.tbi ${file2%.vcf}_reheader_sorted.vcf.gz.tbi || { echo "Failed to remove temporary files"; exit 1; }

# Count the number of individuals in each file
num_individuals_file1=$(bcftools query -l $file1 | wc -l) || { echo "Failed to count individuals in VCF file 1"; exit 1; }
num_individuals_file2=$(bcftools query -l $file2 | wc -l) || { echo "Failed to count individuals in VCF file 2"; exit 1; }
num_individuals_merged=$(bcftools query -l merged_file.vcf | wc -l) || { echo "Failed to count individuals in merged VCF file"; exit 1; }

# Check that the number of individuals in the merged file is as expected
expected_num_individuals=$((num_individuals_file1 + num_individuals_file2))
if [ $num_individuals_merged -eq $expected_num_individuals ]; then
    echo "The number of individuals in the merged file is as expected: $num_individuals_merged"
    echo "Number of individuals in file 1: $num_individuals_file1"
    echo "Number of individuals in file 2: $num_individuals_file2"
else
    echo "Error: The number of individuals in the merged file ($num_individuals_merged) is not as expected ($expected_num_individuals)"
    exit 1
fi
