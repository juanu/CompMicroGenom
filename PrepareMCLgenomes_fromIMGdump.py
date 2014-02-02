#!/Users/juan/anaconda/bin/python
#Created on 7/29/2013

__author__ = 'Juan A. Ugalde'


#I will rewrite some of the code later. This class should go into a separate file. Right now I'm not
#using it, but it will be useful later


class GffEntry():
    """
    This is used to read a gff line, and store the appropiate information into an object
    """
    def __init__(self, gff_line):
        entries = gff_line.split("\t")
        self.type = entries[2]
        self.start = entries[3]
        self.stop = entries[4]
        self.strand = entries[6]
        self.info = entries[8]

    def annot(self, annot_value):

        """
        Returns the value of the annotation entry on the gff line
        the usual entries are ID, locus_tag and product
        """
        annot_dic = {key: result for (key, result) in [entry.split("=") for entry in self.info]}

        try:
            return annot_dic[annot_value]
        except KeyError:
            return None



#Run the script
import os
import argparse
from collections import defaultdict
from Common import PrepareData

program_description = "This script takes a folder with gbk files from RAST, and output four files for" \
                      "each gbk: protein, nucleotide, genome and annotation\n"

parser = argparse.ArgumentParser(description=program_description)

parser.add_argument("-g", "--genome_list", type=str,
                    help="Tabular file with the list the genome to include and the unique prefix to use\n", required=True)
parser.add_argument("-i", "--input_folder", type=str,
                    help="location of the input folder", required=True)
parser.add_argument("-o", "--output_directory", type=str,
                    help="Output folder for the modified fasta files\n", required=True)

args = parser.parse_args()

#Make output directory
if not os.path.exists(args.output_directory):
    os.makedirs(args.output_directory)
    os.makedirs(args.output_directory + "/nucleotide")
    os.makedirs(args.output_directory + "/protein")
    os.makedirs(args.output_directory + "/genome")
    os.makedirs(args.output_directory + "/annotation")
    os.makedirs(args.output_directory + "/coords")

genome_dictionary, total_genome_count = PrepareData.read_genome_list(args.genome_list)

for prefix in genome_dictionary:
    print "Processing genome: %s" % prefix
    img_id = genome_dictionary[prefix]

    organism_folder = args.input_folder + "/" + img_id
    file_prefix = organism_folder + "/" + img_id

    #Create the output files
    nucleotide_file = open(args.output_directory + "/nucleotide/" + prefix + ".fna", 'w')
    aminoacid_file = open(args.output_directory + "/protein/" + prefix + ".fasta", 'w')
    genome_file = open(args.output_directory + "/genome/" + prefix + ".fna", 'w')
    annotation_file = open(args.output_directory + "/annotation/" + prefix + ".txt", 'w')
    coords_file = open(args.output_directory + "/coords/" + prefix + ".txt", 'w')

    #Generate the output files with the sequences
    input_nucleotide = file_prefix + ".genes.fna"
    input_aa = file_prefix + ".genes.faa"
    input_genome = file_prefix + ".fna"

    PrepareData.modify_fasta(input_nucleotide, prefix, nucleotide_file)
    PrepareData.modify_fasta(input_aa, prefix, aminoacid_file)
    PrepareData.modify_fasta(input_genome, prefix, genome_file)

    #Generate the output of the annotation
    input_cog = file_prefix + ".cog.tab.txt"
    input_pfam = file_prefix + ".pfam.tab.txt"
    input_ko = file_prefix + ".ko.tab.txt"
    input_gff = file_prefix + ".gff"

    #Store the COG, PFAM, KO, and product annotation
    #Format Gene_ID Description(COG, PFAM, KO) Annotation
    annotation_summary = defaultdict(lambda: defaultdict(str))
    annotation_defs = defaultdict()

    for line in open(input_cog, 'r'):
        if line.strip():
            if line.startswith("gene_oid"):
                continue
            else:
                line = line.rstrip()
                elements = line.split("\t")
                annotation_summary[elements[0]]["COG"] = elements[9]
                annotation_defs[elements[9]] = elements[10]

    if os.path.exists(input_pfam):

        for line in open(input_pfam, 'r'):
            if line.strip():
                if line.startswith("gene_oid"):
                    continue
                else:
                    line = line.rstrip()
                    elements = line.split("\t")
                    annotation_summary[elements[0]]["PFAM"] = elements[8]
                    annotation_defs[elements[8]] = elements[9]

    for line in open(input_ko, 'r'):
        if line.strip():
            if line.startswith("gene_oid"):
                continue
            else:
                line = line.rstrip()
                elements = line.split("\t")
                annotation_summary[elements[0]]["KO"] = elements[9]

                try:
                    annotation_defs[elements[9]] = elements[10]
                except IndexError:
                    continue

    for line in open(input_gff, 'r'):
        if line.strip():
            if line.startswith("#"):
                continue
            else:
                line = line.rstrip()
                elements = line.split("\t")

                feature_type = elements[2]

                #annotation_summary[protein_id]["feature_type"] = feature_type

                if feature_type == "CRISPR":
                    continue

                try:
                    features = elements[8]
                except IndexError:
                    continue

                protein_id = None
                product_desc = None


                for feature in features.split(";"):
                    if feature.startswith("ID"):
                        name, protein_id = feature.split("=")
                    if feature.startswith("product"):
                        info = feature.split("=")
                        product_desc = info[1]

                if protein_id is not None:
                    annotation_summary[protein_id]["Feature_type"] = elements[2]
                    annotation_summary[protein_id]["start"] = elements[3]
                    annotation_summary[protein_id]["end"] = elements[4]
                    annotation_summary[protein_id]["contig"] = elements[0]

                if product_desc is not None:
                    annotation_summary[protein_id]["Product"] = product_desc

    for protein in annotation_summary:
        for annotation_type in annotation_summary[protein]:

            function = annotation_summary[protein][annotation_type]

            output_line = None

            if annotation_type == "Product" or annotation_type == "Feature_type":
                output_line = protein + "\t" + annotation_type + "\t" + function
                annotation_file.write(output_line + "\n")

            else:
                if annotation_defs.has_key(function):
                    output_line = protein + "\t" + annotation_type + "\t" + function + "\t" + annotation_defs[function]
                else:
                    output_line = protein + "\t" + annotation_type + "\t" + function
                annotation_file.write(output_line + "\n")

        coords_file.write(annotation_summary[protein]["contig"] + "\t" + protein + "\t" + annotation_summary[protein]["start"] + "\t" + annotation_summary[protein]["end"] + "\n")


    nucleotide_file.close()
    aminoacid_file.close()
    genome_file.close()
    annotation_file.close()
    coords_file.close()