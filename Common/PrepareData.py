#!/usr/local/bin/python
#Created on 9/8/13

__author__ = 'Juan A. Ugalde'

def read_genome_list(input_file):
    """
    This function reads a genome list. The first column is the identifier of the genome depending on the annotation
    source used (IMG, JGI or Genbank), the second column is the genome name (Full species name) and the third column
    is the prefix used for this genome in the OrthoMCL analysis
    """
    import sys

    genome_count = 0
    genome_info = {}
    genome_name_list = []
    genome_prefix_list = []

    for line in open(input_file, 'r'):
        if line.strip():
            line = line.rstrip()
            element = line.split("\t")

            #Check for duplicates
        if element[0] in genome_info.keys():  # Duplicate genome ID
            print "Duplicate genome id found:  " + line
            sys.exit("Check for duplicates")

        elif element[1] in genome_name_list:  # Duplicate genome name
            print "Duplicate genome name found: " + line
            sys.exit("Check for duplicates")

        elif element[2] in genome_prefix_list:  # Duplicate prefix name
            print "Duplicate prefix found: " + line
            sys.exit("Check for duplicates")

        else:
            genome_info[element[2]] = element[0]
            genome_count += 1
            genome_name_list.append(element[1])
            genome_prefix_list.append(element[2])

    return genome_info, genome_count


def modify_fasta(input_fasta_file, new_fasta_prefix, output_file):
    """
The input is a fasta file. This will rename the file with the new ID,
and modify each ID in each entry of the fasta file
"""

    from Bio import SeqIO  # Import from tools to read fasta, from Biopython

    #Rename the fasta
    for record in SeqIO.parse(input_fasta_file, "fasta"):
        output_file.write(">" + new_fasta_prefix + "|" + record.id + "\n")
        output_file.write(str(record.seq) + "\n")