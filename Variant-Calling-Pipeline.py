#!/usr/bin/python
import sys
import argparse
import subprocess
import logging
import os

__author__ = "Josh Haddock"
__email__ = "jrhaddock1@gmail.com"

# Set up logging
logging.basicConfig(level=logging.INFO, format='%(asctime)s - %(levelname)s - %(message)s')

def run_command(command):
    """Runs a shell command and checks for errors."""
    logging.info(f"Running command: {command}")
    result = subprocess.run(command, shell=True, stdout=subprocess.PIPE, stderr=subprocess.PIPE)
    if result.returncode != 0:
        logging.error(f"Error running command: {command}\n{result.stderr.decode()}")
        raise RuntimeError(f"Error running command: {command}\n{result.stderr.decode()}")
    return result.stdout.decode()

def bwa_mem(bwa_ref, input_sample_name):
    """Run BWA-MEM alignment."""
    command = f"bwa mem {bwa_ref} {input_sample_name}_1.fq.gz {input_sample_name}_2.fq.gz > {input_sample_name}.sam"
    run_command(command)

def samtools_sort(input_sample_name):
    """Sort SAM file using Samtools."""
    command = f"samtools sort {input_sample_name}.sam -o {input_sample_name}_s.bam"
    run_command(command)

def samtools_mpileup_and_bcftools_call(sam_ref, input_sample_name):
    """Run Samtools mpileup and BCFtools call to generate VCF."""
    command = f"samtools mpileup -Ou -f {sam_ref} {input_sample_name}_s.bam | bcftools call -vmO v -o {input_sample_name}.vcf"
    run_command(command)

def filter_variants(input_sample_name):
    """Filter VCF file using BCFtools."""
    input_vcf = f"{input_sample_name}.vcf"
    output_vcf = f"{input_sample_name}.filtered.vcf"
    command = f"bcftools filter -O v -o {output_vcf} -s LOWQUAL -e '%QUAL<20 || DP<10' {input_vcf}"
    run_command(command)
    return output_vcf

def annotate_variants(input_sample_name, annovar_db):
    """Annotate variants using ANNOVAR."""
    input_vcf = f"{input_sample_name}.filtered.vcf"
    output_prefix = f"{input_sample_name}.annotated"
    command = f"table_annovar.pl {input_vcf} {annovar_db} -buildver hg19 -out {output_prefix} -remove -protocol refGene,cytoBand,gnomad30_genome -operation g,r,f -nastring . -vcfinput"
    run_command(command)
    return f"{output_prefix}.hg19_multianno.vcf"

def calculate_average_qual(vcf_file):
    """Calculate the average QUAL score from a VCF file."""
    total_qual = 0
    total_count = 0

    with open(vcf_file, 'r') as vcf_reader:
        for line in vcf_reader:
            if not line.startswith("#"):
                parts = line.strip().split("\t")
                qual = float(parts[5])
                total_qual += qual
                total_count += 1

    if total_count > 0:
        average_qual = total_qual / total_count
        return average_qual
    else:
        return 0  # Return 0 if no records with QUAL values are found

def parse_args():
    """Standard argument parsing."""
    parser = argparse.ArgumentParser(description="My Variant Calling Pipeline Script")
    parser.add_argument('-i', '--input_sample_name', type=str, required=True, help='Base sample name')
    parser.add_argument('-d', '--annovar_db', type=str, required=True, help='Path to ANNOVAR database')
    return parser.parse_args()

def main():
    """Parses args, then runs the variant calling pipeline."""
    args = parse_args()
    bwa_ref = "/home/jrhaddock1/tiny-test-data/genomes/Hsapiens/hg19/bwa/hg19.fa"
    sam_ref = "/home/jrhaddock1/tiny-test-data/genomes/Hsapiens/hg19/seq/hg19.fa"
    input_sample_name = args.input_sample_name
    annovar_db = args.annovar_db

    try:
        # Run BWA-MEM
        logging.info("Running BWA-MEM...")
        bwa_mem(bwa_ref, input_sample_name)
        
        # Sort the SAM file with Samtools
        logging.info("Sorting SAM file with Samtools...")
        samtools_sort(input_sample_name)
        
        # Run Samtools mpileup and BCFtools call
        logging.info("Running Samtools mpileup and BCFtools call...")
        samtools_mpileup_and_bcftools_call(sam_ref, input_sample_name)

        # Filter variants
        logging.info("Filtering variants...")
        filtered_vcf = filter_variants(input_sample_name)

        # Annotate variants
        logging.info("Annotating variants...")
        annotated_vcf = annotate_variants(input_sample_name, annovar_db)

        # Calculate the average QUAL score
        logging.info("Calculating average QUAL score...")
        average_qual = calculate_average_qual(annotated_vcf)
        logging.info(f"Average QUAL score: {average_qual:.2f}")

    except Exception as err:
        logging.error(f"An error occurred: {err}")
        sys.exit(1)

if __name__ == "__main__":
    main()
