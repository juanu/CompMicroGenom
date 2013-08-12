#!/usr/local/bin/python
#Created on 8/1/13

__author__ = 'Juan A. Ugalde'


def create_sequence_dictionary(fasta_folder):
    """
    This function creates a dictionary that where the key is the ID of the sequence, and the value is it sequence
    """
    import os
    from Bio import SeqIO
    fasta_list = os.listdir(fasta_folder)
    sequence_dictionary = {}

    for infile in fasta_list:
        input_file = fasta_folder + "/" + infile
        fasta_dict = SeqIO.index(input_file, "fasta")

        sequence_dictionary.update(fasta_dict)

    return sequence_dictionary

def get_cluster_information(input_cluster_file):
    cluster_file = open(input_cluster_file, 'r')

    clusters = {}

    for line in cluster_file:
        line = line.rstrip()

        cluster_id, gene_list = line.split("\t")
        clusters[cluster_id] = gene_list

    return clusters

if __name__ == '__main__':
    from cogent import LoadSeqs, DNA, PROTEIN
    from cogent.app.mafft import align_unaligned_seqs as mafft_align_unaligned_seqs
    from cogent.core.genetic_code import DEFAULT as standard_code
    from cogent.app.fasttree import build_tree_from_alignment as build_tree_fasttree
    import os
    import argparse

    program_description = "Script that takes a list of clusters, and the sequence information of the genes. " \
                          "The output is tha alignment of each cluster, and a tree (FastTree)"

    parser = argparse.ArgumentParser(description=program_description)

    parser.add_argument("-c", "--cluster_file", type=str, help="Cluster file", required=True)
    parser.add_argument("-n", "--fasta_nuc_directory", type=str, help="Output folder", required=True)
    parser.add_argument("-o", "--output_directory", type=str, help="Output folder", required=True)

    args = parser.parse_args()

    #Create the output folders
    protein_unaligned_folder = args.output_directory + "/protein_unaligned"
    protein_alignment_folder = args.output_directory + "/protein_alignment"
    dna_unaligned_folder = args.output_directory + "/dna_unaligned"
    dna_aligned_folder = args.output_directory + "/dna_aligned"
    tree_folder = args.output_directory + "/protein_trees"

    folder_list = [args.output_directory, protein_unaligned_folder, protein_alignment_folder, dna_unaligned_folder,
                   dna_aligned_folder, tree_folder]

    for folder in folder_list:
        if not os.path.exists(folder):
            os.makedirs(folder)

    #Get the cluster information
    cluster_information = get_cluster_information(args.cluster_file)

    #Create the sequence dictionary
    dna_sequence_dic = create_sequence_dictionary(args.fasta_nuc_directory)

    #Iterate over each cluster and generate the alignments

    frameshift_cases = []
    inframe_stops = []
    clusters_too_short = []
    nucleotide_not_found = []

    for cluster in cluster_information:

        protein_list = cluster_information[cluster].split(",")
        curated_protein_list = {}

        for protein in protein_list:
            genome_id, protein_id = protein.split("|")

            #Check if the sequence exists, if not store it and keep going

            if not protein in dna_sequence_dic:
                nucleotide_not_found.append(protein)
                continue

            sequence_with_stop_codons = DNA.makeSequence(dna_sequence_dic[protein])

            #Check if the sequence is the right one, and check for in frame stops
            #It seems that in JGI annotation, when scaffolds are joined, the resulted proteins do not match
            #the DNA sequence
            #Right now, I'll just remove those sequences, and deal with that later

            if len(sequence_with_stop_codons) % 3 == 0:
                seq_no_stop_codon = sequence_with_stop_codons.withoutTerminalStopCodon()

                #Chec for inframe stop codons
                stops_frame = standard_code.getStopIndices(seq_no_stop_codon, start=0)

                if len(stops_frame) > 0:
                    inframe_stops.append([cluster, genome_id, protein_id])

                else:
                    curated_protein_list[protein] = seq_no_stop_codon

            else:
                frameshift_cases.append([cluster, genome_id, protein_id])

        if len(curated_protein_list) < 3:  # Only take those clusters with 3 sequences or more
            clusters_too_short.append(cluster)
            continue

        #Alignments and output data

        unaligned_DNA = LoadSeqs(data=curated_protein_list, moltype=DNA, aligned=False)

        unaligned_AA = unaligned_DNA.getTranslation()

        #Generate alignments using muscle
        aligned_AA = mafft_align_unaligned_seqs(unaligned_AA, PROTEIN)

        #Replace the aminoacid sequences with the nucleotide sequence
        aligned_DNA = aligned_AA.replaceSeqs(unaligned_DNA)

        #Make protein trees using FastTree
        protein_tree = build_tree_fasttree(aligned_AA, PROTEIN)

        #Tree file
        tree_output = open(tree_folder + "/" + cluster + ".tre", 'w')

        tree_output.write(protein_tree.getNewick(with_distances=True))

        #Write the unaligned and aligned sequences
        unaligned_DNA.writeToFile(dna_unaligned_folder + "/" + cluster + ".fna", format="fasta")
        unaligned_AA.writeToFile(protein_unaligned_folder + "/" + cluster + ".faa", format="fasta")
        aligned_DNA.writeToFile(dna_aligned_folder + "/" + cluster + ".fna", format="fasta")
        aligned_AA.writeToFile(protein_alignment_folder + "/" + cluster + ".faa", format="fasta")

        tree_output.close()

    #Print log files

    logfile = open(args.output_directory + "/logfile.txt", 'w')
    file_frameshifts = open(args.output_directory + "/frameshifts.txt", 'w')
    file_inframe_stops = open(args.output_directory + "/inframe_stops.txt", 'w')
    file_short_clusters = open(args.output_directory + "/small_clusters.txt", 'w')
    file_nuc_not_found = open(args.output_directory + "/nucleotide_not_found.txt", 'w')

    logfile.write("Total number of analyzed cluster: %d\n" % len(cluster_information))
    logfile.write("Clusters with less than 3 sequences (after cleaning): %d\n" % len(clusters_too_short))
    logfile.write("Total sequences with frameshifts: %d\n" % len(frameshift_cases))
    logfile.write("Sequences with inframe stops: %d\n" % len(inframe_stops))
    logfile.write("Proteins not found in the nucleotide sequences: %d\n" % len(nucleotide_not_found))

    for item in frameshift_cases:
        file_frameshifts.write("\t".join(item) + "\n")

    for item in inframe_stops:
        file_inframe_stops.write("\t".join(item) + "\n")

    file_short_clusters.write("\n".join(clusters_too_short))

    file_nuc_not_found.write("\n".join(nucleotide_not_found))

    logfile.close()
    file_frameshifts.close()
    file_inframe_stops.close()
    file_nuc_not_found.close()