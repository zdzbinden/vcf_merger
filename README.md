# VCF File Processing and Merging Script

This script is designed to process and merge two VCF (Variant Call Format) files that have been aligned using the same reference genome. The purpose of this script is to combine VCF files with different sets of individuals (i.e., non-overlapping samples, but similar loci).

We do this with **bcftools merge**, but this requires a few things of our VCF files..

The VCF files need:
- SNPs ordered by contig/chromosome
- Contig names and lengths provided in the header
- Compressed with bgzip
- Indexed with tabix

This script takes care of all those things!

The script filters the VCF files to retain only the lines that contain contig names present in a provided .fai file. The .fai file is an index file for the reference genome fasta file, and it contains the names, sizes, and other information about the contigs in the reference genome. Each line of the .fai file represents one contig and is formatted as follows: `<contig name>\t<size>\t<offset>\t<linebases>\t<linewidth>`.

After filtering, the script replaces the headers of the VCF files with combined headers, sorts, compresses, and indexes the VCF files, and finally merges the VCF files. The script checks for errors at each step and stops if any command fails. It also checks that the required programs are installed before running.

## Requirements

The script requires the following programs to be installed:

- `bcftools`
- `bgzip`
- `tabix`

These can be installed using conda or other package managers.

## Usage

The script takes as input two VCF files and a .fai file. The .fai file should contain the contig names that are present in the VCF files. The script can be run from the command line as follows:

```bash
./script.sh file1.vcf file2.vcf file.fai
```

Replace `file1.vcf`, `file2.vcf`, and `file.fai` with the names of your actual files.

The script outputs a merged VCF file named `merged_file.vcf`.

Here's a brief guide on how to use the script:

1. **Clone the repository from GitHub to your local machine.** This can be done using the `git clone` command followed by the URL of the repository. For example:

    ```
    git clone https://github.com/zdzbinden/vcf_merger
    ```

This command will create a new directory in your current location with the same name as the repository. This directory will contain all the files from the repository.

2. **Navigate to the directory containing the script.** Use the `cd` command to change directories. For example:

    ```
    cd vcf_merger
    ```

3. **Make the script executable.** This can be done using the `chmod` command followed by `u+x` and the name of the script. For example:

    ```
    chmod u+x vcf_merger.sh
    ```

4. **Run the script.** This can be done using the `./` prefix followed by the name of the script and any necessary arguments. For example:

    ```
    ./vcf_merger.sh file1.vcf file2.vcf file.fai
    ```

    Replace `file1.vcf`, `file2.vcf`, and `file.fai` with the actual names of your input files (if they are in the same directory as vef_merger.sh) or the paths to those files. This command runs the script with the specified input files.

Remember to replace the placeholders in the commands with your actual file names or paths. If you encounter any errors, make sure that your input files are in the correct format and that you have the necessary permissions to read and write in the directory.

## Error Checking

The script checks for errors at each step and stops if any command fails. It reports an error message indicating what went wrong. It also checks that the required programs are installed before running.

At the end of the script, it checks that the number of individuals in the merged file is as expected (the sum of the individuals in the two input VCF files). If the number of individuals in the merged file is not as expected, the script reports an error and exits with a non-zero status.

## Temporary Files

The script creates several temporary files during its execution. These files are removed as soon as they are no longer needed. The original input files are not modified.

---
### License

This software is provided for free: you can redistribute it and/or modify it under the terms of the GNU Public License as published by the Free Software Foundation. You should have received a copy of the GNU Public License with the software. If not, see:  [www.gnu.org/licenses](http://www.gnu.org/licenses)

The author claims no liability nor resposibility for the functionality of this software.
