#!/bin/bash

# This script takes as input two VCF files and a .fai file, filters the VCF files to retain only the lines that contain contig names present in the .fai file, replaces the headers of the VCF files with combined headers, sorts, compresses, and indexes the VCF files, and finally merges the VCF files. It checks for errors at each step and stops if any command fails. It also checks that the required programs are installed before running.

# Usage: ./script.sh file1.vcf file2.vcf file.fai

# Get start time
start=$(date +%s)

# Check that the required programs are installed
command -v bcftools >/dev/null 2>&1 || { echo -e "bcftools is not installed. Aborting.\n"; exit 1; }
command -v bgzip >/dev/null 2>&1 || { echo -e "bgzip is not installed. Aborting.\n"; exit 1; }
command -v tabix >/dev/null 2>&1 || { echo -e "tabix is not installed. Aborting.\n"; exit 1; }

echo -e "All required programs are installed.\n"

# Input files
file1=$1
file2=$2
fai_file=$3

# Extract the contig names from the .fai file and save them to a temporary file
awk '{print $1}' $fai_file > contigs.txt || { echo -e "Failed to extract contig names from .fai file\n"; exit 1; }

echo -e "Contig names extracted from .fai file.\n"

# Filter the VCF files to retain only the header lines and the lines that contain contig names present in the .fai file
awk 'BEGIN{while(getline < "contigs.txt") contigs[$1]} /^#/ || ($1 in contigs)' $file1 > ${file1%.vcf}_filtered.vcf || { echo -e "Failed to filter VCF file 1\n"; exit 1; }
awk 'BEGIN{while(getline < "contigs.txt") contigs[$1]} /^#/ || ($1 in contigs)' $file2 > ${file2%.vcf}_filtered.vcf || { echo -e "Failed to filter VCF file 2\n"; exit 1; }

echo -e "VCF files filtered to retain only the header lines and the lines that contain contig names present in the .fai file.\n"

# Extract the original headers from the VCF files, separating the final line
grep "^##" $file1 > ${file1%.vcf}_original.hdr || { echo -e "Failed to extract original header from VCF file 1\n"; exit 1; }
grep "^#" $file1 | grep -v "^##" > ${file1%.vcf}_final.hdr || { echo -e "Failed to extract final header line from VCF file 1\n"; exit 1; }
grep "^##" $file2 > ${file2%.vcf}_original.hdr || { echo -e "Failed to extract original header from VCF file 2\n"; exit 1; }
grep "^#" $file2 | grep -v "^##" > ${file2%.vcf}_final.hdr || { echo -e "Failed to extract final header line from VCF file 2\n"; exit 1; }

echo -e "Original headers extracted from the VCF files.\n"

# Create a header file from the .fai file, only including contigs present in the VCF files
awk 'BEGIN{while(getline < "contigs.txt") contigs[$1]} {if ($1 in contigs) print "##contig=<ID="$1",length="$2">"}' $fai_file > contigs.hdr || { echo -e "Failed to create header file from .fai file\n"; exit 1; }

echo -e "Header file created from the .fai file.\n"

# Append the contig information to the headers of the VCF files
cat contigs.hdr >> ${file1%.vcf}_original.hdr || { echo -e "Failed to append contig information to header of VCF file 1\n"; exit 1; }
cat contigs.hdr >> ${file2%.vcf}_original.hdr || { echo -e "Failed to append contig information to header of VCF file 2\n"; exit 1; }

echo -e "Contig information appended to the headers of the VCF files.\n"

# Append the final line of the original header to the new header
cat ${file1%.vcf}_final.hdr >> ${file1%.vcf}_original.hdr || { echo -e "Failed to append final header line to header of VCF file 1\n"; exit 1; }
cat ${file2%.vcf}_final.hdr >> ${file2%.vcf}_original.hdr || { echo -e "Failed to append final header line to header of VCF file 2\n"; exit 1; }

echo -e "Final line of the original header appended to the new header.\n"

# Replace the headers of the VCF files with the combined headers
bcftools reheader -h ${file1%.vcf}_original.hdr -o ${file1%.vcf}_reheader.vcf ${file1%.vcf}_filtered.vcf || { echo -e "Failed to replace header of VCF file 1\n"; exit 1; }
bcftools reheader -h ${file2%.vcf}_original.hdr -o ${file2%.vcf}_reheader.vcf ${file2%.vcf}_filtered.vcf || { echo -e "Failed to replace header of VCF file 2\n"; exit 1; }

echo -e "Headers of the VCF files replaced with the combined headers.\n"

# Remove filtered VCF files after use
rm ${file1%.vcf}_filtered.vcf ${file2%.vcf}_filtered.vcf || { echo -e "Failed to remove filtered VCF files\n"; exit 1; }

# Change the ID column information for the loci in each input VCF file to be a combination of the information in the CHROM column and the POS column
bcftools annotate -Ov -I '%CHROM\_%POS' ${file1%.vcf}_reheader.vcf > ${file1%.vcf}_annotated.vcf || { echo -e "Failed to annotate VCF file 1\n"; exit 1; }
bcftools annotate -Ov -I '%CHROM\_%POS' ${file2%.vcf}_reheader.vcf > ${file2%.vcf}_annotated.vcf || { echo -e "Failed to annotate VCF file 2\n"; exit 1; }

echo -e "VCF file IDs annotated.\n"

# Remove reheader VCF files after use
rm ${file1%.vcf}_reheader.vcf ${file2%.vcf}_reheader.vcf || { echo -e "Failed to remove reheader VCF files\n"; exit 1; }

# Normalize to combine multiple records created by multi-allelic loci
bcftools norm -Ov -m+any ${file1%.vcf}_annotated.vcf -o ${file1%.vcf}_normed.vcf || { echo -e "Failed to normalize VCF file 1\n"; exit 1; }
bcftools norm -Ov -m+any ${file2%.vcf}_annotated.vcf -o ${file2%.vcf}_normed.vcf|| { echo -e "Failed to normalize VCF file 2\n"; exit 1; }

echo -e "Input files are normalized to remove duplicate records.\n"

# Check for unique IDs in the VCF files
num_ids_file1=$(awk '!/^#/ {print $3}' ${file1%.vcf}_normed.vcf | sort | uniq | wc -l)
num_rows_file1=$(awk '!/^#/ {print $3}' ${file1%.vcf}_normed.vcf | wc -l)
if [ $num_ids_file1 -ne $num_rows_file1 ]; then
    echo -e "***************************************************"
    echo -e "Warning: IDs in VCF file 1 are not unique";
    echo -e "***************************************************"
fi

num_ids_file2=$(awk '!/^#/ {print $3}' ${file2%.vcf}_normed.vcf | sort | uniq | wc -l)
num_rows_file2=$(awk '!/^#/ {print $3}' ${file2%.vcf}_normed.vcf | wc -l)
if [ $num_ids_file2 -ne $num_rows_file2 ]; then
    echo -e "***************************************************"
    echo -e "Warning: IDs in VCF file 2 are not unique";
    echo -e "***************************************************"
fi

echo -e "Checked for unique IDs in the VCF files.\n"

# Sort the VCF files
bcftools sort ${file1%.vcf}_normed.vcf -Ov -o ${file1%.vcf}_annotated_sorted.vcf || { echo -e "Failed to sort VCF file 1\n"; exit 1; }
bcftools sort ${file2%.vcf}_normed.vcf -Ov -o ${file2%.vcf}_annotated_sorted.vcf || { echo -e "Failed to sort VCF file 2\n"; exit 1; }

echo -e "VCF files sorted.\n"



# Compress and index the VCF files
bgzip -c ${file1%.vcf}_annotated_sorted.vcf > ${file1%.vcf}_annotated_sorted.vcf.gz || { echo -e "Failed to compress sorted VCF file 1\n"; exit 1; }
bgzip -c ${file2%.vcf}_annotated_sorted.vcf > ${file2%.vcf}_annotated_sorted.vcf.gz || { echo -e "Failed to compress sorted VCF file 2\n"; exit 1; }
tabix -p vcf ${file1%.vcf}_annotated_sorted.vcf.gz || { echo -e "Failed to index compressed VCF file 1\n"; exit 1; }
tabix -p vcf ${file2%.vcf}_annotated_sorted.vcf.gz || { echo -e "Failed to index compressed VCF file 2\n"; exit 1; }

echo -e "VCF files compressed and indexed.\n"



# Merge the VCF files
bcftools merge -m id ${file1%.vcf}_annotated_sorted.vcf.gz ${file2%.vcf}_annotated_sorted.vcf.gz -Ov -o merged_file.vcf || { echo -e "Failed to merge VCF files\n"; exit 1; }

echo -e "VCF files merged.\n"

# Check for unique IDs in the merged VCF file
num_ids_merged=$(awk '!/^#/ {print $3}' merged_file.vcf | sort | uniq | wc -l)
num_rows_merged=$(awk '!/^#/ {print $3}' merged_file.vcf | wc -l)
if [ $num_ids_merged -ne $num_rows_merged ]; then
    echo -e "***************************************************"
    echo -e "Warning: IDs in the merged VCF file are not unique\n"
    echo -e "***************************************************"
fi

echo -e "Checked for unique IDs in the merged VCF file.\n"

# Count the number of loci in each input VCF file
num_loci_file1=$(awk '!/^#/ {print $3}' ${file1%.vcf}_annotated_sorted.vcf | wc -l) || { echo -e "Failed to count loci in VCF file 1\n"; exit 1; }
num_loci_file2=$(awk '!/^#/ {print $3}' ${file2%.vcf}_annotated_sorted.vcf | wc -l) || { echo -e "Failed to count loci in VCF file 2\n"; exit 1; }

echo -e "Counted the number of loci in each input VCF file: File 1 = $num_loci_file1; File 2 = $num_loci_file2.\n"

# Count the number of shared loci between the two input VCF files
bcftools isec -n=2 -Ov ${file1%.vcf}_annotated_sorted.vcf.gz ${file2%.vcf}_annotated_sorted.vcf.gz -w 1 > intersect.txt || { echo -e "Failed to count shared loci between VCF files\n"; exit 1; }

# Calculate num_shared_loci by piping the temporary file into grep and wc
num_shared_loci=$(grep -v "^#" intersect.txt | wc -l) || { echo -e "Failed to count shared loci between VCF files\n"; exit 1; }

echo -e "Counted the number of shared loci between the two input VCF files: $num_shared_loci.\n"

# Calculate the number of unique loci in each input VCF file
num_unique_loci_file1=$(bcftools isec -C -n-1 -w1 -Ov  ${file1%.vcf}_annotated_sorted.vcf.gz ${file2%.vcf}_annotated_sorted.vcf.gz | grep -v "^#" | wc -l) || { echo -e "Failed to count unique loci in File 1\n"; exit 1; }
num_unique_loci_file2=$(bcftools isec -C -n-1 -w1 -Ov  ${file2%.vcf}_annotated_sorted.vcf.gz ${file1%.vcf}_annotated_sorted.vcf.gz | grep -v "^#" | wc -l) || { echo -e "Failed to count unique loci in File 2\n"; exit 1; }

echo -e "Calculated the number of unique loci in each input VCF file: File 1 = $num_unique_loci_file1; File 2 = $num_unique_loci_file2.\n"

# Count the number of loci in the merged VCF file
num_loci_merged=$(cat merged_file.vcf | grep -v "^#" | wc -l) || { echo -e "Failed to count loci in merged VCF file\n"; exit 1; }

echo -e "Counted the number of loci in the merged VCF file: $num_loci_merged.\n"

# Remove annotated VCF files after use
rm ${file1%.vcf}_annotated.vcf ${file2%.vcf}_annotated.vcf ${file1%.vcf}_normed.vcf ${file2%.vcf}_normed.vcf || { echo -e "Failed to remove annotated VCF files\n"; exit 1; }
# Remove sorted VCF files after use
rm ${file1%.vcf}_annotated_sorted.vcf ${file2%.vcf}_annotated_sorted.vcf || { echo -e "Failed to remove sorted VCF files\n"; exit 1; }
# Clean up temporary files
rm contigs.txt ${file1%.vcf}_original.hdr ${file1%.vcf}_final.hdr ${file2%.vcf}_original.hdr ${file2%.vcf}_final.hdr contigs.hdr || { echo -e "Failed to clean up temporary files\n"; exit 1; }
rm ${file1%.vcf}_annotated_sorted.vcf.gz ${file1%.vcf}_annotated_sorted.vcf.gz.tbi ${file2%.vcf}_annotated_sorted.vcf.gz ${file2%.vcf}_annotated_sorted.vcf.gz.tbi || { echo -e "Failed to clean up temporary files\n"; exit 1; }

echo -e "Temporary files cleaned up.\n"

# Change the ID column information for the loci in output VCF file to be a combination of the information in the CHROM column and the POS columm
bcftools annotate -Ov -I '%CHROM\_%POS' merged_file.vcf > combined_alignment.vcf || { echo -e "Failed to annotate merged VCF file 1\n"; exit 1; }
rm merged_file.vcf

echo -e "Merged IDs annotated.\n"


# Create Intersection VCF
bgzip -c combined_alignment.vcf > combined_alignment.vcf.gz || { echo -e "Failed to compress final file 1\n"; exit 1; }
tabix -p vcf combined_alignment.vcf.gz || { echo -e "Failed to index final file 1\n"; exit 1; }
# Extract shared loci from the merged VCF file
bcftools view -R intersect.txt -Ov -o combined_intersection.vcf combined_alignment.vcf.gz


# Peep the fixed info of the vcf file
awk '!/^##/ {print $1, $2, $3, $4, $5}' combined_alignment.vcf  | head

echo

echo -e "Script completed successfully.\n"

# Get end time and report runtime
end=$(date +%s)
runtime=$((end-start))

echo "Job ran in: $runtime seconds"

# Convert the runtime to hours, minutes, and seconds
hours=$((runtime / 3600))
minutes=$(( (runtime % 3600) / 60))
seconds=$((runtime % 60))

echo "Job ran in: $hours hours $minutes minutes $seconds seconds"

