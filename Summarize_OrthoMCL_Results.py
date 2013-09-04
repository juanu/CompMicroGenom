#!/usr/local/bin/python
#Created on 8/1/13

__author__ = 'Juan Ugalde'


def read_genome_list(input_file):
    """
    This function reads a genome list. The first column is the identifier of the genome depending on the annotation
    source used (IMG, JGI or Genbank), the second column is the genome name (Full species name) and the third column
    is the prefix used for this genome in the OrthoMCL analysis
    """
    import sys

    genome_count = 0  # Numbers of genomes analyzed
    genome_info = {}  #Store the information for the genomes
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


def get_protein_info(genome_list, fasta_directory):
    """
    This function goes into a folder with fasta files (.fasta) of the proteins and get the information for the proteins
    in each genome, including the ID and the length of the protein
    """
    from Bio import SeqIO
    from collections import defaultdict

    proteins_in_genomes = defaultdict(list)  # output list of the proteins that are present in each genome
    protein_length = defaultdict(int)  # length of each protein. Assumes a unique ID for each protein

    files_read_counter = 0  # Count the number of files being read

    fasta_files = [fasta_directory + "/" + fasta + ".fasta" for fasta in genome_list]

    for fasta in fasta_files:
        for record in SeqIO.parse(fasta, "fasta"):
            proteins_in_genomes[record.id.split('|')[0]].append(record.id)
            protein_length[record.id] = int(len(record.seq))

        files_read_counter += 1

    return proteins_in_genomes, protein_length, files_read_counter


def get_orthomcl_results(cluster_file, genome_list):
    """
    This function reads the results of the OrthoMCL clustering, and returns a dictionary with the information.
    The format is:
    ClusterID: protein1 protein2 ....
    """

    from collections import defaultdict

    orthomcl_results = open(cluster_file, 'r')

    cluster_dictionary = defaultdict(list)  # Dictionary with the clusters. Stores only the clusters that will be used
    unique_proteins_genome_count = defaultdict(int)  # Proteins that are not in any cluster
    proteins_in_cluster = set()  # Proteins in cluster. Needed to search for proteins that are absent
    total_cluster_count = 0
    removed_clusters = 0

    for line in orthomcl_results:
        total_cluster_count += 1
        line = line.strip('\n')
        ids_proteins = line.split(": ")
        proteins = ids_proteins[1].split(" ")

        clean_protein_list = []  # This is used to remove proteins from genomes not in the list

        for genome in genome_list:  # Adding the proteins that we need
            [clean_protein_list.append(protein) for protein in [x for x in proteins if x.startswith(genome)]]

        #Now I need to evaluate those clusters that now are unique and zero
        if len(clean_protein_list) == 0:
            removed_clusters += 1
            continue

        if len(clean_protein_list) == 1:
            unique_proteins_genome_count[clean_protein_list[0].split("|")[0]] += 1
            removed_clusters += 1
            continue

        for protein in clean_protein_list:
            cluster_dictionary[ids_proteins[0]].append(protein)
            proteins_in_cluster.add(protein)

    return cluster_dictionary, proteins_in_cluster, unique_proteins_genome_count, total_cluster_count, removed_clusters


def read_group_files(group_file):
    """
    Reads a file with the group list, and returns a dictionary containing
    the name of the group and a list with the genomes that are in that group.
    The first column is the prefix of the genome, the second the associated group
    """

    from collections import defaultdict

    genome_groups = defaultdict(list)

    for line in open(group_file, 'r'):
        line = line.rstrip()
        element = line.split("\t")

        genome_groups[element[0]].append(element[1])

    return genome_groups

#Here we go into the main calculations


def get_unique_seqs_genome(sequences_in_genomes, total_set_sequences, sequence_lengths, min_length):
    """
    This module takes a dictionary with their genome and proteins, and a set of proteins
    and look for proteins that are not in the total set
    """

    from collections import defaultdict

    large_unique_sequences = defaultdict(list)
    short_unique_sequences = defaultdict(list)
    processed_sequences = 0

    for genome in sequences_in_genomes:
        for seq_name in sequences_in_genomes[genome]:
            processed_sequences += 1

            if seq_name in total_set_sequences:
                continue
            else:
                if sequence_lengths[seq_name] > min_length:
                    large_unique_sequences[genome].append(seq_name)
                else:
                    short_unique_sequences[genome].append(seq_name)

    return large_unique_sequences, short_unique_sequences, processed_sequences


def seqs_shared_clusters(cluster_dic, genome_dictionary):
    """
    This module takes the cluster dictionary and the genome dictionary and generates outputs which includes:
    Clusters that are unique to each genome
    Clusters that are shared across all genomes (single and multiple copy)

    """
    from collections import defaultdict

    unique_clusters = defaultdict(list)  # Dictionary with the unique clusters
    shared_single_clusters = []  # List with clusters that are shared and single copy
    shared_multiple_clusters = []  # List with clusters that are shared and in multiple copies

    genomes_in_matrix = sorted(genome_dictionary.keys())
    header = ["Cluster_ID"]
    header.extend(sorted(genome_dictionary.keys()))
    all_clusters_matrix = [header]

    for cluster in cluster_dic:

        genome_list = [protein.split("|")[0] for protein in cluster_dic[cluster]]  # Create a list with the genomes

        count = {x: genome_list.count(x) for x in genome_list}  # Count the occurences

        #Create the matrix
        cluster_matrix = [cluster]

        for genome in genomes_in_matrix:

            if genome in genome_list:
                cluster_matrix.append(count[genome])
            else:
                cluster_matrix.append(0)

        #print cluster_matrix
        all_clusters_matrix.append(cluster_matrix)

        if len(count) == 1:
            unique_clusters[genome_list[0]].append(cluster)

        elif len(count) == len(genome_dictionary.keys()):
            if sum(count.itervalues()) == len(genome_dictionary.keys()):
                shared_single_clusters.append(cluster)
            else:
                shared_multiple_clusters.append(cluster)

    return unique_clusters, shared_single_clusters, shared_multiple_clusters, all_clusters_matrix


def clusters_in_groups(clusters, groups):
    """

    """
    from collections import defaultdict

    import itertools
    unique_group_clusters = defaultdict(list)  # Count the unique clusters in each group

    combination_clusters = defaultdict(list)

    #Create inverted dictionary
    genome_group_info = dict()

    for group in groups:
        for genome in groups[group]:
            genome_group_info[genome] = group

    for cluster in clusters:

        group_count = defaultdict(lambda: defaultdict(int))

        for protein in clusters[cluster]:
            genome_id = protein.split("|")[0]
            group_for_genome = genome_group_info[genome_id]

            group_count[group_for_genome][genome_id] += 1

        ##Unique clusters for each group

        if len(group_count) == 1:
            for group in group_count:
                if len(group_count[group]) == len(groups[group]):
                    unique_group_clusters[group].append(cluster)

                else:  # I could add something here to count the number of proteins not unique
                    pass

        #Shared, all possible combinations
        else:
            for combination_count in range(2, len(group_count) + 1):
                for group_combinations in itertools.combinations(groups.keys(), combination_count):

                    # Check that all the genomes in the group are represented
                    group_check = 0

                    for i in range(0, len(group_combinations)):
                        if len(group_count[group_combinations[i]]) == len(groups[group_combinations[i]]):
                            continue
                        else:
                            group_check = 1

                    if group_check == 0:
                        combination_clusters[group_combinations].append(cluster)

    return unique_group_clusters, combination_clusters

if __name__ == '__main__':
    import os
    import sys
    import argparse

    #Create the options and program description
    program_description = "This script summarize the results of orthoMCL, and create several summary files." \
                          " The inputs are:" \
                          "-List of clusters, generated by orthoMCL" \
                          "-A genome list" \
                          "-A folder with fasta files" \
                          "- An optional group file, to group genomes. For example, all genomes from the same species "

    parser = argparse.ArgumentParser(description=program_description)

    parser.add_argument("-l", "--genome_list_index", type=str,
                        help="File with the genome list. Format GenomeID, FullName, ShortName", required=True)
    parser.add_argument("-c", "--cluster_file", type=str, help="Ortholog file, generated by OrthoMCL", required=True)
    parser.add_argument("-f", "--fasta_aa_directory", type=str, help="Directory with the fasta files", required=True)
    parser.add_argument("-g", "--group_information", type=str, help="Group file")
    parser.add_argument("-o", "--output_directory", type=str, help="Output directory", required=True)

    args = parser.parse_args()

    #Create the output directory
    if not os.path.exists(args.output_directory):
        os.makedirs(args.output_directory)

    #Create a log file
    run_summary = open(args.output_directory + "/logfile.txt", 'w')

    #####Read the genome list
    genome_id_dictionary, genome_count = read_genome_list(args.genome_list_index)

    run_summary.write("Genomes in the genome list: %d" % genome_count + "\n")

    ######Read the cluster information, and check that everything is ok
    #cluster_information, set_of_proteins_in_clusters, unique_cluster_count, total_clusters, removed_clusters = \
    #    get_orthomcl_results(args.cluster_file, [i for i in genome_id_dictionary.itervalues()])

    cluster_information, set_of_proteins_in_clusters, unique_cluster_count, total_clusters, removed_clusters = \
        get_orthomcl_results(args.cluster_file, genome_id_dictionary.keys())

    run_summary.write("Total number of clusters: %d" % len(cluster_information) + "\n")
    run_summary.write("Total number of protein in clusters: %d" % len(set_of_proteins_in_clusters) + "\n")
    run_summary.write("Total number of removed clusters (not present in the genome file): %d" % removed_clusters + "\n")

    #Check the counts, to see if everything is going ok
    if total_clusters - removed_clusters != len(cluster_information):
        sys.exit("The number of removed clusters clusters plus the retained clusters, "
                 "doesn't match the total of original clusters in the file")

    #####Read the fasta file
    #dic_protein_in_genomes, dic_protein_length, files_read_counter = \
    #    get_protein_info([i for i in genome_id_dictionary.itervalues()], args.fasta_aa_directory)

    dic_protein_in_genomes, dic_protein_length, files_read_counter = \
        get_protein_info(genome_id_dictionary.keys(), args.fasta_aa_directory)

    run_summary.write("Total fasta files: %d" % files_read_counter + "\n")
    run_summary.write("Total number of proteins in the fasta files: %d" % len(dic_protein_length) + "\n")

    #####Read the genome groups (if present)
    genome_groups = {}
    if args.group_information:
        genome_groups = read_group_files(args.group_information)
        run_summary.write("Total defined groups: %d" % len(genome_groups) + "\n")

    #####################################################
    #Look for unique protein in each genome
    selected_unique_proteins, removed_unique_proteins, total_number_proteins = \
        get_unique_seqs_genome(dic_protein_in_genomes, set_of_proteins_in_clusters, dic_protein_length, 50)

    #Check that everything looks ok
    check_number = 0
    for value in selected_unique_proteins.itervalues():
        check_number += len(value)
    for value in removed_unique_proteins.itervalues():
        check_number += len(value)

    if total_number_proteins - len(set_of_proteins_in_clusters) - check_number != 0:
        print "Total number of proteins:" + str(total_number_proteins)
        print "Total number of proteins in clusters:" + str(len(set_of_proteins_in_clusters))
        print "check number:" + str(check_number)
        sys.exit("Failed checkpoint. The number of unique proteins and proteins in "
                 "clusters does not match the total number of proteins")

    #Print the output files
    count_unique_proteins = open(args.output_directory + "/count_unique_sequences.txt", 'w')
    list_unique_proteins = open(args.output_directory + "/list_unique_sequences.txt", 'w')

    count_unique_proteins.write("Genome\tSelected\tTooShort\n")

    for genome in selected_unique_proteins:
        count_unique_proteins.write(genome + "\t" + str(len(selected_unique_proteins[genome])) + "\t" +
                                    str(len(removed_unique_proteins[genome])) + "\n")

        for protein in selected_unique_proteins[genome]:
            list_unique_proteins.write(genome + "\t" + protein.split("|")[1] + "\n")

    count_unique_proteins.close()
    list_unique_proteins.close()

    ############################
    ##Get the clusters shared between genomes and unique clusters to each genome
    matrix_output = open(args.output_directory + "/matrix_output.txt", 'w')
    list_unique_clusters = open(args.output_directory + "/list_unique_clusters.txt", 'w')
    count_unique_clusters = open(args.output_directory + "/count_unique_clusters.txt", 'w')
    list_shared_single_copy_clusters = open(args.output_directory + "/list_single_copy_clusters.txt", 'w')
    list_shared_multiple_copy_clusters = open(args.output_directory + "/list_shared_multiple_copy.txt", 'w')

    unique_clusters, shared_single_clusters, shared_multiple_clusters, all_clusters_matrix = \
        seqs_shared_clusters(cluster_information, genome_id_dictionary)

    #Print counters
    run_summary.write("Number of shared single copy clusters: %d" % len(shared_single_clusters) + "\n")
    run_summary.write("Number of shared multiple copy clusters: %d" % len(shared_multiple_clusters) + "\n")

    #Print the outputs
    matrix_output.write("\n".join(["\t".join(map(str, r)) for r in all_clusters_matrix]))  # Matrix output

    # Unique clusters per genome (duplicate or paralogs?)
    count_unique_clusters.write("Genome\tNumber of Clusters\n")
    for genome in unique_clusters:
        count_unique_clusters.write(genome + "\t" + str(len(unique_clusters[genome])) + "\n")

        for cluster in unique_clusters[genome]:
            list_unique_clusters.write(genome + "\t" + cluster + "\t"
                                       + ",".join(protein for protein in cluster_information[cluster]) + "\n")

    # Single copy shared clusters

    for cluster in shared_single_clusters:
        list_shared_single_copy_clusters.write(cluster + "\t" + ",".join(cluster_information[cluster]) + "\n")

    # Multiple copy shared clusters

    for cluster in shared_multiple_clusters:
        list_shared_multiple_copy_clusters.write(cluster + "\t" + ",".join(cluster_information[cluster]) + "\n")

    matrix_output.close()
    list_unique_clusters.close()
    count_unique_clusters.close()
    list_shared_single_copy_clusters.close()
    list_shared_multiple_copy_clusters.close()

    ###Save the cluster information
    list_all_clusters = open(args.output_directory + "/list_all_clusters.txt", 'w')
    for cluster in cluster_information:
        list_all_clusters.write(cluster + "\t" + ",".join(cluster_information[cluster]) + "\n")

    list_all_clusters.close()

    ###############
    ##Get clusters shared by groups
    if args.group_information:

        unique_group_clusters, combination_clusters = clusters_in_groups(cluster_information, genome_groups)

        list_unique_clusters_group = open(args.output_directory + "/list_unique_clusters_group.txt", 'w')
        list_all_group_combinations = open(args.output_directory + "/list_all_group_combinations.txt", 'w')
        count_group_results = open(args.output_directory + "/count_groups.txt", 'w')

        for group in unique_group_clusters:
            protein_count = sum(len(cluster_information[cluster]) for cluster in unique_group_clusters[group])

            count_group_results.write(group + "\t" +
                                      str(len(unique_group_clusters[group])) + "\t" + str(protein_count) + "\n")

            for cluster in unique_group_clusters[group]:
                list_unique_clusters_group.write(group + "\t" + cluster + "\t" + ",".join(cluster_information[cluster]))

        count_group_results.write("\n")

        for combination in combination_clusters:
            combination_name = "-".join(combination)
            protein_count = sum(len(cluster_information[cluster]) for cluster in combination_clusters[combination])

            count_group_results.write(combination_name + "\t" +
                                      str(len(combination_clusters[combination])) + "\t" + str(protein_count) + "\n")

            for cluster in combination_clusters[combination]:
                list_all_group_combinations.write(combination_name + "\t"
                                                  + cluster + "\t" + ",".join(cluster_information[cluster]) + "\n")

        count_group_results.close()
        list_unique_clusters_group.close()
        list_all_group_combinations.close()
        run_summary.close()