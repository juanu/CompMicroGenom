#!/Users/juan/anaconda/bin/python
#Created on 6/20/2013

__author__ = 'Juan A. Ugalde'


def read_gbk_list(input_file):

    """
    Reads a list of file names and prefix to be used for the further steps. The input is a tabular file, where
    the first column is the organisms name in the NCBI format (downloaded directly from the FTP server), and the
    second column is a unique prefix to use for that genome

    """

    genome_count = 0
    genome_info = {}
    for line in open(input_file, 'r'):
        if line.strip():
            line = line.rstrip()
            element = line.split("\t")
            genome_info[element[0]] = element[1]
            genome_count += 1

    return genome_info, genome_count


def get_gbk_files(folder_name):
    import os
    gbk_files = []
    for input_gbk in os.listdir(folder_name):
        if input_gbk.endswith(".gbk"):
            gbk_files.append(input_gbk)

    return gbk_files


def get_cog_note(note):
    """

    :param note:
    :return:
    """
    import re
    cog_search = re.search('(COG\d+).*', note)

    cog_number = None

    if cog_search is not None:
        cog_number = cog_search.group(1)

    return cog_number

#Run the script
if __name__ == '__main__':
    import os
    import argparse
    from Bio import SeqIO

    program_description = "This script takes a folder with gbk files from the NCBI server, and output four files for" \
                          "each gbk: protein, nucleotide, genome and annotation"

    parser = argparse.ArgumentParser(description=program_description)

    parser.add_argument("-g", "--gbk_list", type=str,
                        help="Tabular file with the list of the gbk files and the unique prefix to use", required=True)
    parser.add_argument("-i", "--input_folder", type=str,
                        help="location of the gbk folder", required=True)
    parser.add_argument("-o", "--output_directory", type=str,
                        help="Output folder for the modified fasta files", required=True)

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

    for gbk_folder_name in genome_dictionary:
        prefix = genome_dictionary[gbk_folder_name]

        organism_folder = args.input_folder + "/" + gbk_folder_name

        #Create the output files
        nucleotide_file = open(args.output_directory + "/nucleotide/" + prefix + ".fna", 'w')
        aminoacid_file = open(args.output_directory + "/protein/" + prefix + ".fasta", 'w')
        genome_file = open(args.output_directory + "/genome/" + prefix + ".fna", 'w')
        annotation_file = open(args.output_directory + "/annotation/" + prefix + ".txt", 'w')

        #Get the list of gbk files in the folder
        gbk_files = get_gbk_files(organism_folder)

        for gbk in gbk_files:
            gbk_location = organism_folder + "/" + gbk
            input_handle = open(gbk_location, "r")

            for record in SeqIO.parse(input_handle, "genbank"):
                scaf_id = prefix + "|" + record.id
                genome_file.write(">" + scaf_id + "\n" + str(record.seq) + "\n")

                for feature in record.features:
                    if feature.type == "CDS":

                        protein_id = feature.qualifiers["protein_id"][0]

                        cds_id = prefix + "|" + protein_id
                        nucleotide_sequence = feature.extract(record.seq)
                        aminoacid_sequence = feature.qualifiers["translation"][0]

                        nucleotide_file.write(">" + cds_id + "\n" + str(nucleotide_sequence) + "\n")
                        aminoacid_file.write(">" + cds_id + "\n" + str(aminoacid_sequence) + "\n")

                        #Create the annotation output
                        if "product" in feature.qualifiers:
                            annotation_file.write(protein_id + "\tProduct\t" + feature.qualifiers["product"][0] + "\n")
                        if "EC_number" in feature.qualifiers:
                            annotation_file.write(protein_id + "\tEC_number\t" + feature.qualifiers["EC_number"][0] + "\n")
                        if "note" in feature.qualifiers:

                            cog_number = get_cog_note(feature.qualifiers["note"][0])

                            if cog_number is not None:
                                annotation_file.write(protein_id + "\tCOG\t" + cog_number + "\n")

            input_handle.close()

        nucleotide_file.close()
        aminoacid_file.close()
        genome_file.close()
        annotation_file.close()