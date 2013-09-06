
def run_paml_per_group(groups, alignment, tree, output_dir, working_dir):
    """
    This is to take each defined group, modify the tree, and the run PAML
    """
    from Bio import Phylo
    import re

    cluster_tree = Phylo.read(tree, "newick")

    #Names have a pipe sign (|) with the organism|protein_id.
    # I need to keep track of everything to replace in the final tree
    clades_in_tree = {str(clade).split("|")[0]: str(clade).split("|")[1] for clade in cluster_tree.get_terminals()}

    clade_results = dict()

    #Iterate on each group
    for group in groups:

        #Check that all the branches are present on the tree (and is not the only branch
        if set(groups[group]).issubset(set(clades_in_tree.keys())) and len(clades_in_tree.keys()) > len(groups[group]):
            dict_new_clade_names = {name + "|" + clades_in_tree[name]: name + "|" + clades_in_tree[name] + " #1"
                                    for name in groups[group]}

            #Replace the names in the tree and save the tree

            old_tree_informtion = open(tree).read()

            new_tree_information = multiple_replace(dict_new_clade_names, old_tree_informtion)

            group_tree = working_dir + "/" + group + ".tre"

            new_tree_file = open(group_tree, 'w')
            new_tree_file.write(new_tree_information)
            new_tree_file.close()

            #Run model for the new tree

            paml_results = run_paml_site_branch_models(alignment, group_tree, output_dir, working_dir)


            clade_results[group] = paml_results


        else:
            clade_results[group] = None

    return clade_results

if __name__ == '__main__':
    import os
    import argparse
    from collections import defaultdict
    import multiprocessing
    from multiprocessing import Manager

    program_description = "Script that takes a list of clusters and runs PAML (codeml). The model used is a branch-site" \
                          "with relaxed test (MA vs M1a). "

    parser = argparse.ArgumentParser(description=program_description)

    parser.add_argument("-c", "--cluster_list", type=str, help="Cluster file", required=True)
    parser.add_argument("-n", "--cluster_folder", type=str, help="Output folder", required=True)
    parser.add_argument("-g", "--groups", type=str, help="Group constrains", required=True)
    parser.add_argument("-o", "--output_directory", type=str, help="Output folder", required=True)
    parser.add_argument("-p", "--num_processors", type=int, help="Number of processors to use", required=True)
    parser.add_argument("-f", "--fdr", help="Perform false discovery rate")

    args = parser.parse_args()

    #Check for the output folder and also create the temporal folder
    #I'm using the PID to create the temporary folder, which should allow multiple instances of the script to run

    if not os.path.exists(args.output_directory):
        os.makedirs(args.output_directory)

    #Read the cluster file and group file
    clusters_to_analyze = [line.rsplit()[0] for line in open(args.cluster_list) if line.strip()]

    group_constrains = defaultdict(list)

    for line in open(args.groups):
        if line.strip():
            line = line.rstrip()
            group_constrains[line.split("\t")[0]].append(line.split("\t")[1])

    #Result and output files
    manager = Manager()
    cluster_paml_results = manager.list([])
    groups_no_data = manager.list([])

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

        p = multiprocessing.Process(target=cluster_analysis, args=(chunk, args.cluster_folder,
        group_constrains, args.output_directory, temporal_folder, cluster_paml_results, groups_no_data))

        jobs.append(p)

    #Run the jobs
    [proc.start() for proc in jobs]
    [proc.join() for proc in jobs]

    #Print the results
    output_file = open(args.output_directory + "/paml_results.txt", 'w')
    no_results_file = open(args.output_directory + "/no_results.txt", 'w')

    for result in cluster_paml_results:
        output_file.write("\t".join(str(x) for x in result) + "\n")

    for entry in groups_no_data:
        print entry + "\n"

    output_file.close()
    no_results_file.close()

    #False discovery
    if args.fdr:
        from operator import itemgetter

        total_tests = len(cluster_paml_results)  # Total number of performed tests

        position = 1
        prev_adjusted_pvalue = 0

        sorted_results = sorted(cluster_paml_results, key=itemgetter(4))

        for entry in sorted_results:
            adjusted_pvalue = float(entry[4]) * (total_tests / position)

            #If the value is greater than 1, we set as one (0 < p < 1)
            adjusted_pvalue = min(adjusted_pvalue, 1)

            #Check that the value is not greater than the previous one
            adjusted_pvalue = max(adjusted_pvalue, prev_adjusted_pvalue)

            prev_adjusted_pvalue = adjusted_pvalue
            position += 1

            entry.insert(5, adjusted_pvalue)

        corrected_pvalue_file = open(args.output_directory + "/paml_results_corrected.txt", "w")

        for result in sorted_results:
            corrected_pvalue_file.write("\t".join(str(x) for x in result) + "\n")


