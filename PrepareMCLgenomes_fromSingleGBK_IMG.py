#!/Users/juan/anaconda/bin/python
#Created on 8/7/2013

__author__ = 'Juan A. Ugalde'

#Run the script
if __name__ == '__main__':
    import os
    import argparse
    from Bio import SeqIO

    program_description = "This script takes a single gbk from IMG and its corresponding annotation file " \
                          "(or gene information) and generates the output needed for the next steps of the pipeline"

    parser = argparse.ArgumentParser(description=program_description)

    parser.add_argument("-g", "--gbk_file", type=str,
                        help="Tabular file with the list of the gbk files and the unique prefix to use", required=True)
    parser.add_argument("-i", "--annotation_file", type=str,
                        help="location of the gbk folder", required=True)
    parser.add_argument("-p", "--prefix", required=True, help="Genome prefix to use", type=str)
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

    #Create the output files
    nucleotide_file = open(args.output_directory + "/nucleotide/" + args.prefix + ".fna", 'w')
    aminoacid_file = open(args.output_directory + "/protein/" + args.prefix + ".fasta", 'w')
    genome_file = open(args.output_directory + "/genome/" + args.prefix + ".fna", 'w')
    annotation_file = open(args.output_directory + "/annotation/" + args.prefix + ".txt", 'w')

    #Read the gbk file
    input_gbk = open(args.gbk_file, 'r')

    locustag_to_id = dict()

    for line in open(args.annotation_file, 'r'):
        if line.strip():
            if line.startswith("gene_oid"):
                continue
            else:
                line = line.rstrip()

                elements = line.split("\t")
                locustag_to_id[elements[1]] = elements[0]

                if elements[2].startswith("pfam"):
                    annotation_file.write("%s\tPFAM\t%s\t%s\n" % (elements[0], elements[2], elements[3]))

                if elements[2].startswith("Product_name"):
                    annotation_file.write("%s\tProduct\t%s\n" % (elements[0], elements[4]))

                if elements[2].startswith("Locus_type"):
                    annotation_file.write("%s\tFeature_type\t%s\n" % (elements[0], elements[4]))

                if elements[2].startswith("COG") and elements[2] is not "COG_category":
                    annotation_file.write("%s\tCOG\t%s\t%s\n" % (elements[0], elements[2], elements[3]))

                if elements[2].startswith("KO"):
                    annotation_file.write("%s\tKO\t%s\t%s\n" % (elements[0], elements[2], elements[3]))

    for record in SeqIO.parse(input_gbk, "genbank"):
        scaf_id = args.prefix + "|" + record.id
        genome_file.write(">" + scaf_id + "\n" + str(record.seq) + "\n")

        for feature in record.features:
            if feature.type == "CDS":

                locus_id = feature.qualifiers["locus_tag"][0]

                protein_id = locustag_to_id[locus_id]
                cds_id = args.prefix + "|" + protein_id
                nucleotide_sequence = feature.extract(record.seq)
                aminoacid_sequence = feature.qualifiers["translation"][0]

                nucleotide_file.write(">" + cds_id + "\n" + str(nucleotide_sequence) + "\n")
                aminoacid_file.write(">" + cds_id + "\n" + str(aminoacid_sequence) + "\n")

    input_gbk.close()

    #Store the annotation

    nucleotide_file.close()
    aminoacid_file.close()
    genome_file.close()
    annotation_file.close()