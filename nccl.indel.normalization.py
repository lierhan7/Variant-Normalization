#!/usr/bin/env python

import argparse
import pysam

def query_reference_sequence(genome_fasta_file, chromosome, start, end):
    genome_fasta = pysam.FastaFile(genome_fasta_file)
    sequence = genome_fasta.fetch(chromosome, int(start), int(end))
    sequence = sequence.upper()

    return sequence

def nccl_indel_normalization(chromosome, position, ref_allele, alt_allele, strand, genome_fasta_file):

    if len(ref_allele) > 1 and len(alt_allele) == 1:
        # deletion
        deletion_sequence = ref_allele[1:]      # CCTA
        deletion_start = int(position) + 1      # 10183821
        deletion_end = int(position) + len(deletion_sequence)   # 10183824
        reference_sequence = query_reference_sequence(genome_fasta_file, chromosome, deletion_start - 1, deletion_start + len(deletion_sequence) * 20 - 1)      # CCTACCCA
        reference_sequence_deletion = reference_sequence[len(deletion_sequence):]
        deletion_sequence_normalize = deletion_sequence
        deletion_start_normalize = deletion_start
        deletion_end_normalize = deletion_end
        if strand == '+':
            idx = 0
            while True:
                idx += 1
                deletion_sequence_temp = reference_sequence[idx:idx+len(deletion_sequence)]
                deletion_start_temp = deletion_start + idx
                deletion_end_temp = deletion_start + idx + len(deletion_sequence_temp) - 1
                reference_sequence_deletion_temp = reference_sequence[:idx] + reference_sequence[idx+len(deletion_sequence_temp):]
                if reference_sequence_deletion_temp == reference_sequence_deletion:
                    deletion_sequence_normalize = deletion_sequence_temp
                    deletion_start_normalize = deletion_start_temp
                    deletion_end_normalize = deletion_end_temp
                else:
                    break

        return chromosome, deletion_start_normalize, deletion_end_normalize, deletion_sequence_normalize, '-'
    elif len(ref_allele) == 1 and len(alt_allele) > 1:
        # insertion
        insertion_sequence = alt_allele[1:]
        insertion_site = int(position)
        reference_sequence = query_reference_sequence(genome_fasta_file, chromosome, insertion_site, insertion_site + len(insertion_sequence) * 20 - 1)
        reference_sequence_insertion = insertion_sequence + reference_sequence 
        insertion_sequence_normalize = insertion_sequence
        insertion_site_normalize = insertion_site
        if strand == '+':
            idx = 0
            while True:
                idx += 1
                insertion_sequence_temp = reference_sequence_insertion[idx:idx+len(insertion_sequence)]
                insertion_site_temp = insertion_site + idx
                reference_sequence_insertion_temp = reference_sequence[:idx] + insertion_sequence_temp + reference_sequence[idx:]
                if reference_sequence_insertion_temp == reference_sequence_insertion:
                    insertion_sequence_normalize = insertion_sequence_temp
                    insertion_site_normalize = insertion_site_temp
                else:
                    break
                
        return chromosome, insertion_site_normalize, '-', '-', insertion_sequence_normalize
    elif len(ref_allele) == 1 and len(alt_allele) == 1:
        # snv
        return chromosome, position, position, ref_allele, alt_allele
    elif len(ref_allele) > 1 and len(alt_allele) > 1:
        # complex
        return chromosome, position, int(position) + len(ref_allele) - 1, ref_allele, alt_allele

# testing
#chromosome = 'chr4'
#position = '55593608'
#ref_allele = 'GGTTGTT'
#alt_allele = 'G'
#strand = '+'
#genome_fasta_file = 'hg19.fa'
#chromosome, start_normalize, end_normalize, ref_allele_normalize, alt_allele_normalize = nccl_indel_normalization(chromosome, position, ref_allele, alt_allele, strand, genome_fasta_file)
#print(chromosome, start_normalize, end_normalize, ref_allele_normalize, alt_allele_normalize)

if __name__ == "__main__":
    args_parser = argparse.ArgumentParser(description="normalize variant by nccl rule")
    args_parser.add_argument("--input", dest="input_file", required=True, help="input file to normalize")
    args_parser.add_argument("--output", dest="output_file", required=True, help="output normalized file")
    args_parser.add_argument("--ref", dest="reference_genome", required=True, help="reference genome fasta file")

    args = args_parser.parse_args()
    input_file = args.input_file
    output_file = args.output_file
    reference_genome = args.reference_genome

    with open(input_file, 'r') as readfile, open(output_file, 'w') as writefile:
        lines = readfile.readlines()
        for line in lines:
            line = line.strip()
            if line.startswith('#'):
                writefile.write(line+'\n')
            else:
                line_parts = line.split('\t')
                chromosome = line_parts[0]
                position = line_parts[1]
                ref_allele = line_parts[2]
                alt_allele = line_parts[3]
                strand = line_parts[5]
                genome_fasta_file = reference_genome
                chromosome, start_normalize, end_normalize, ref_allele_normalize, alt_allele_normalize = nccl_indel_normalization(chromosome, position, ref_allele, alt_allele, strand, genome_fasta_file)
                writefile.write(chromosome+'\t'+str(start_normalize)+'\t'+ref_allele_normalize+'\t'+alt_allele_normalize+'\t'+'\t'.join(line_parts[5:])+'\n')
