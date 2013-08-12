#!/Users/juan/anaconda/bin/python
#Created on 6/22/2013

__author__ = 'Juan A. Ugalde'


def read_gbk_list(input_file):

    """
    Reads a list of file names and prefix to be used for the further steps. The input is a tabular file, where
    the first column is the full path to the file, and the second column is the prefix to use later.

    """

    genome_count = 0
    genome_info = {}
    for line in open(input_file, 'r'):
        if line.strip():
            line = line.rstrip()
            element = line.split("\t")
            genome_info[element[1]] = element[0]
            genome_count += 1

    return genome_info, genome_count


#Run the script
if __name__ == '__main__':
    import os
    import argparse
    from Bio import SeqIO

    program_description = "This script takes a folder with gbk files from RAST, and output four files for" \
                          "each gbk: protein, nucleotide, genome and annotation\n"

    parser = argparse.ArgumentParser(description=program_description)

    parser.add_argument("-g", "--gbk_list", type=str,
                        help="Tabular file with the list of the gbk files and the unique prefix to use\n", required=True)
    #parser.add_argument("-i", "--input_folder", type=str,
                        #help="location of the gbk folder", required=True)
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

    #Read the genome list, and create the dictionary
    genome_dictionary, total_genome_count = read_gbk_list(args.gbk_list)

    for genome_prefix in genome_dictionary:
        gbk_file = open(genome_dictionary[genome_prefix], "r")

        #Create the output files
        nucleotide_file = open(args.output_directory + "/nucleotide/" + genome_prefix + ".fna", 'w')
        aminoacid_file = open(args.output_directory + "/protein/" + genome_prefix + ".fasta", 'w')
        genome_file = open(args.output_directory + "/genome/" + genome_prefix + ".fna", 'w')
        annotation_file = open(args.output_directory + "/annotation/" + genome_prefix + ".txt", 'w')

        protein_count = 1

        for record in SeqIO.parse(gbk_file, "genbank"):

            scaf_id = genome_prefix + "|" + record.id
            genome_file.write(">" + scaf_id + "\n" + str(record.seq) + "\n")

            for feature in record.features:
                    if feature.type == "CDS":

                        protein_id = genome_prefix + "_" + str(protein_count)
                        protein_count += 1

                        cds_id = genome_prefix + "|" + protein_id
                        nucleotide_sequence = feature.extract(record.seq)
                        aminoacid_sequence = feature.qualifiers["translation"][0]

                        nucleotide_file.write(">" + cds_id + "\n" + str(nucleotide_sequence) + "\n")
                        aminoacid_file.write(">" + cds_id + "\n" + str(aminoacid_sequence) + "\n")

                        #Create the annotation output
                        if "product" in feature.qualifiers:
                            annotation_file.write(protein_id + "\tProduct\t" + feature.qualifiers["product"][0] + "\n")
                        if "EC_number" in feature.qualifiers:
                            annotation_file.write(protein_id + "\tEC_number\t" + feature.qualifiers["EC_number"][0] + "\n")

                        #Missing the COG information

        gbk_file.close()

        nucleotide_file.close()
        aminoacid_file.close()
        genome_file.close()
        annotation_file.close()