#!/usr/local/bin/python
#Created on 8/20/13

__author__ = 'Juan Ugalde'


def run_paml_per_group(groups, alignment, tree, output_dir, working_dir):
    """
    This function take the group, alignment, tree and folder information and runs a paml analysis
    on each defined group.
    The steps needed are to modify the tree to add the #1 that defines the branches in the tree for paml
    and then runs PAML on that tree, using the provided alignment.
    The working dir is important (different from the output dir), because different PAML runs at the same time may
    override each other.
    This is particularly important if running this script in more than one processor
    """
    from Bio import Phylo
    from SelectionAnalysis import paml_run

    cluster_tree = Phylo.read(tree, "newick")  # Read the input tree

    #Names have a pipe sign (|) with the organism|protein_id.
    clades_in_tree_by_gene_id = {str(clade).split("|")[1]: str(clade).split("|")[0]
                                 for clade in cluster_tree.get_terminals()}

    species_in_tree = set(str(clade).split("|")[0] for clade in cluster_tree.get_terminals())

    clade_results = dict()

    #Iterate on each group
    for group in groups:

        #Check that all the branches are present on the tree (and is not the only branch)
        if set(groups[group]).issubset(species_in_tree) and len(species_in_tree) > len(groups[group]):

            dict_new_clade_names = dict()

            for gene_id in clades_in_tree_by_gene_id:
                genome = clades_in_tree_by_gene_id[gene_id]

                if genome in groups[group]:
                    dict_new_clade_names[genome + "|" + gene_id] = genome + "|" + gene_id + " #1"
                else:
                    continue

            #Replace the names in the tree and save the tree
            old_tree_information = open(tree).read()

            new_tree_information = multiple_replace(dict_new_clade_names, old_tree_information)

            group_tree = working_dir + "/" + group + ".tre"

            new_tree_file = open(group_tree, 'w')
            new_tree_file.write(new_tree_information)
            new_tree_file.close()

            #Run model for the new tree
            paml_results = paml_run.ma_m1a(alignment, group_tree, output_dir, working_dir)

            clade_results[group] = paml_results

        else:
            clade_results[group] = None

    return clade_results


def multiple_replace(replace_dict, text):
    """
    Replace string based on dictionary. Taken from:
    http://stackoverflow.com/questions/15175142/how-can-i-do-multiple-substitutions-using-regex-in-python
    """
    import re

    #Create the regular expression from the dictionary keys
    regex = re.compile("(%s)" % "|".join(map(re.escape, replace_dict.keys())))

    #For each match, look-up corresponding value in dictionary
    return regex.sub(lambda mo: replace_dict[mo.string[mo.start():mo.end()]], text)


def cluster_analysis(cluster_list, cluster_folder, group_branches, output_folder, temporal_folder, results, no_data, not_found):
    """
    Function used to run the analysis on the cluster list. It will run PAML for each group, and then it will
    calculate the stats
    """

    from SelectionAnalysis import paml_stats
    from SelectionAnalysis import paml_prepare

    for cluster in cluster_list:
        cluster_file = cluster_folder + "/" + cluster + ".fna"  # Add fna extension

        #Check that the cluster file exists, if not continue
        if not os.path.exists(cluster_file):
            not_found.append(cluster)
            continue

        #Make a new tree, no confidence values in the branches
        new_tree = paml_prepare.run_fasttree(cluster_file, temporal_folder)

        #Make the new alignment, and get information about the alignment
        new_alignment_file, number_sequences, alignment_length = paml_prepare.adjust_alignment(cluster_file, temporal_folder)

        #Run PAML for each branch in the cluster with both models
        paml_site_branch_results = run_paml_per_group(group_branches, new_alignment_file, new_tree,
                                                      output_folder, temporal_folder)

        for group in paml_site_branch_results:

            #Store those clusters and groups that were not analyzed
            if paml_site_branch_results[group] is None:
                no_data.append([cluster, group])

            else:
                pvalue = paml_stats.lrt(paml_site_branch_results[group]["Ma"].get("lnL"),
                                        paml_site_branch_results[group]["M1a"].get("lnL"), 1)

                proportion_sites = float(paml_site_branch_results[group]["Ma"]["site_classes"][2]["proportion"]) + \
                                float(paml_site_branch_results[group]["Ma"]["site_classes"][3]["proportion"])

                average_omega = (float(paml_site_branch_results[group]["Ma"]["site_classes"][2]["branch types"]["foreground"]) +
                                  float(paml_site_branch_results[group]["Ma"]["site_classes"][3]["branch types"]["foreground"])) / 2

                #Store the final results
                #Group, Nseqs, Length, p-value, P1 in Ma, Omega in W
                results.append([cluster, group, number_sequences, alignment_length,
                                              round(pvalue, 3), proportion_sites, average_omega])


if __name__ == '__main__':
    import os
    import argparse
    from collections import defaultdict
    import multiprocessing as mp
    from SelectionAnalysis import paml_run, paml_stats

    program_description = "Script that takes a list of clusters and runs PAML (codeml). The model used is a branch-site" \
                          "with relaxed test (MA vs MA with omega fixed at 1). "

    parser = argparse.ArgumentParser(description=program_description)

    parser.add_argument("-c", "--cluster_list", type=str, help="Cluster file", required=True)
    parser.add_argument("-n", "--cluster_folder", type=str, help="Output folder", required=True)
    parser.add_argument("-g", "--groups", type=str, help="Group constrains", required=True)
    parser.add_argument("-o", "--output_directory", type=str, help="Output folder", required=True)
    parser.add_argument("-p", "--num_processors", type=int, help="Number of processors to use", required=True)
    parser.add_argument("-f", "--fdr", help="Perform false discovery rate")

    args = parser.parse_args()

    #Check for the output folder and also create the temporal folder

    if not os.path.exists(args.output_directory):
        os.makedirs(args.output_directory)

    #Read the cluster file and group file
    clusters_to_analyze = [line.rsplit()[0] for line in open(args.cluster_list) if line.strip()]

    group_constrains = defaultdict(list)  # Define the group of branches to analyze

    for line in open(args.groups):
        if line.strip():
            line = line.rstrip()
            group_constrains[line.split("\t")[0]].append(line.split("\t")[1])

    #Result and output files
    manager = mp.Manager()
    cluster_paml_results = manager.list([])
    groups_no_data = manager.list()
    clusters_not_found = manager.list()

    #Run in parallel, split the list
    num_proc = args.num_processors
    num_chunks = len(clusters_to_analyze) / num_proc
    clusters_chunks = [clusters_to_analyze[i:i+num_chunks] for i in range(0, len(clusters_to_analyze), num_chunks)]
    jobs = []
    i = 1

    #Create the jobs to run
    for chunk in clusters_chunks:
        temporal_folder = args.output_directory + "/temp_" + str(i)
        if not os.path.exists(temporal_folder):
            os.makedirs(temporal_folder)

        i += 1

        p = mp.Process(target=cluster_analysis, args=(chunk, args.cluster_folder,
        group_constrains, args.output_directory, temporal_folder, cluster_paml_results, groups_no_data, clusters_not_found))

        jobs.append(p)

    #Run the jobs
    [proc.start() for proc in jobs]
    [proc.join() for proc in jobs]

    #Print the results
    output_file = open(args.output_directory + "/paml_results.txt", 'w')
    no_results_file = open(args.output_directory + "/no_results.txt", 'w')
    not_found_file = open(args.output_directory + "/clusters_not_present.txt", 'w')

    for result in cluster_paml_results:
        output_file.write("\t".join(str(x) for x in result) + "\n")

    for entry in groups_no_data:
        no_results_file.write("\t".join(entry) + "\n")

    output_file.close()
    no_results_file.close()

    #Run False discovery rate analysis, if chosen
    if args.fdr:
        corrected_pvalue_results = paml_stats.fdr(cluster_paml_results, 4)

        corrected_results_file = open(args.output_directory + "/paml_results_corrected.txt", "w")

        for result in corrected_pvalue_results:
            corrected_results_file.write("\t".join(str(x) for x in result) + "\n")




