#!/Users/juan/anaconda/bin/python
#Created on 7/29/2013

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


#Run the script
if __name__ == '__main__':
    import os
    import argparse
    from collections import defaultdict

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

    genome_dictionary, total_genome_count = read_genome_list(args.genome_list)

    for prefix in genome_dictionary:
        img_id = genome_dictionary[prefix]

        organism_folder = args.input_folder + "/" + img_id
        file_prefix = organism_folder + "/" + img_id

        #Create the output files
        nucleotide_file = open(args.output_directory + "/nucleotide/" + prefix + ".fna", 'w')
        aminoacid_file = open(args.output_directory + "/protein/" + prefix + ".fasta", 'w')
        genome_file = open(args.output_directory + "/genome/" + prefix + ".fna", 'w')
        annotation_file = open(args.output_directory + "/annotation/" + prefix + ".txt", 'w')

        #Generate the output files with the sequences
        input_nucleotide = file_prefix + ".genes.fna"
        input_aa = file_prefix + ".genes.faa"
        input_genome = file_prefix + ".fna"

        modify_fasta(input_nucleotide, prefix, nucleotide_file)
        modify_fasta(input_aa, prefix, aminoacid_file)
        modify_fasta(input_genome, prefix, genome_file)

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
                    annotation_defs[elements[9]] = elements[10]

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

                    features = elements[8]
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
                    output_line = protein + "\t" + annotation_type + "\t" + function + "\t" + annotation_defs[function]
                    annotation_file.write(output_line + "\n")

        nucleotide_file.close()
        aminoacid_file.close()
        genome_file.close()
        annotation_file.close()