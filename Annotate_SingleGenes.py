#Created on 9/8/13
author = "Juan A. Ugalde"

from Common import AnnotationTools
import os
import argparse
from collections import defaultdict

program_description = "This script annotates a list of genes, where the first column is the genome name" \
                      "and the second is the gene ID. The annotation is taken from the folder generated by" \
                      "the import scripts"

parser = argparse.ArgumentParser(description=program_description)

#Arguments
parser.add_argument("-a", "--annotation_folder", type=str,
                    help="Folder with the annotation files", required=True)
parser.add_argument("-c", "--gene_file", type=str,
                    help="Gene file", required=True)
parser.add_argument("-o", "--output_file", type=str,
                    help="Output file", required=True)

args = parser.parse_args()

#Read the gene list

genome_gene_info = defaultdict(list)

for line in open(args.gene_file):
    if line.strip():
        line = line.rstrip()
        genome_gene_info[line.split("\t")[0]].append(line.split("\t")[1])

#Get the annotation information
protein_annotation, function_definitions = \
    AnnotationTools.parse_annotation_folder(genome_gene_info.keys(), args.annotation_folder)

#Print output table
output_file = open(args.output_file, 'w')

for genome in genome_gene_info:
    for protein in genome_gene_info[genome]:

        try:
            product = protein_annotation[protein]["Product"]
        except KeyError:
            product = None

        try:
            COG_number = protein_annotation[protein]["COG"]
        except KeyError:
            COG_number = None

        COG_description = None

        if not COG_number is None:
            COG_description = function_definitions[COG_number]

        try:
            PFAM_number = protein_annotation[protein]["PFAM"]
        except KeyError:
            PFAM_number = None

        PFAM_description = None
        if not PFAM_number is None:
            PFAM_description = function_definitions[PFAM_number]

        output_line = [genome, protein, product, COG_number, COG_description, PFAM_number, PFAM_description]

        output_file.write("\t".join(str(x) for x in output_line) + "\n")








#print protein_annotation
#print function_definitions

