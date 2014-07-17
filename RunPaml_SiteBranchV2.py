#Created on 7/11/2014
__author__ = 'Juan Ugalde'

#TODO
#Clean and document the script


def run_site_branch(cluster_name, treefile, alignment, folder_temp, folder_plots):
    from ete2 import EvolTree
    from ete2.treeview.layouts import evol_clean_layout
    import os
    from collections import defaultdict
    import math
    from scipy.stats import chi2

    print "Processing cluster: " + cluster_name

    tree = EvolTree(treefile)
    tree.link_to_alignment(alignment, alg_format="fasta", nucleotides=True)

    #Create temporal folder
    temp_cluster_folder = folder_temp + "/" + cluster_name

    if not os.path.exists(temp_cluster_folder):
        os.makedirs(temp_cluster_folder)

    tree.workdir = temp_cluster_folder

    #Run M0 as the null model
    tree.run_model("M0")

    #Look at the site selection on each branch

    printed_tree = 0

    i = 0

    #Output list with the results
    output_list = []

    for node in tree.iter_descendants():

        #Mark the tree for the leaf under analysis
        tree.mark_tree([node.node_id], marks=["#1"])

        #Use the node id as folder name
        temp_leaf_name = str(node.node_id)

        print "Processing: " + cluster_name + " " + temp_leaf_name + " " + ",".join(node.get_leaf_names())

        #Run computation of each model.
        #From the notes on ETE:
        # to organize a bit, we name model with the name of the marked node
        # any character after the dot, in model name, is not taken into account
        # for computation. (have a look in /tmp/ete2.../bsA.. directory)

        tree.run_model("bsA." + temp_leaf_name)
        tree.run_model("bsA1." + temp_leaf_name)

        bsA = tree.get_evol_model("bsA." + temp_leaf_name)
        bsA1 = tree.get_evol_model("bsA1." + temp_leaf_name)

        ps_sites = defaultdict()
        total_sites = 0
        sites_over_95 = 0

        for s in range(len(bsA.sites['BEB']['aa'])):
            p_value_site = float(bsA.sites['BEB']['p2'][s])

            if p_value_site > 0.50:
                ps_sites[s] = [bsA.sites['BEB']['aa'][s], bsA.sites['BEB']['p2'][s]]
                total_sites += 1

                if p_value_site > 0.95:
                    sites_over_95 += 1

        #ps = float(tree.get_most_likely("bsA." + temp_leaf_name, "bsA1." + temp_leaf_name))
        rx = float(tree.get_most_likely("bsA1." + temp_leaf_name, "M0"))

        lrt_value = 2 * math.fabs(bsA1.lnL - bsA.lnL)  # LRT test value
        ps = 1 - chi2.cdf(lrt_value, 1)  # p-value based on chi-square

        test_status = None

        #Evidence of positive selection in the branch
        omega_value = float(bsA.classes['foreground w'][2])
        proportion_sites = float(bsA.classes['proportions'][2])

        #Plot file
        plot_file = folder_plots + "/" + cluster_name

        if ps < 0.05 and omega_value > 1:
            #Save plots, both in jpg and svg of the clusters with evidence of positive selection
            test_status = "Positive"

            if printed_tree == 0:
                #tree.render(plot_file + ".svg", layout=evol_clean_layout)
                #tree.render(plot_file + ".jpg", layout=evol_clean_layout)
                printed_tree = 1

            else:
                continue

        elif rx < 0.05 and ps >= 0.05:
            test_status = "Relaxed"

        else:
            #print "no signal"
            test_status = None

        #Remove marks on the tree
        tree.mark_tree(map(lambda x: x.node_id, tree.get_descendants()), marks=[''] * len(tree.get_descendants()),
                       verbose=False)

        result_entry = [cluster_name, node.node_id, omega_value, proportion_sites, ps, test_status,
                        total_sites, sites_over_95, ",".join(node.get_leaf_names())]

       # print result_entry
        #print ps_sites
        #node_results[node.node_id] = [result_entry, ps_sites]
        node_result = [result_entry, ps_sites]

        output_list.append(node_result)

    return output_list


if __name__ == '__main__':
    import os
    import argparse
    from collections import defaultdict
    import multiprocessing
    from ete2 import EvolTree
    from ete2.treeview.layouts import evol_clean_layout
    import shutil
    import sys

    program_description = "Script that takes a list of clusters, their trees and nucleotide alignments and run the" \
                          "site-branch test on them." \
                          "This script will test all the branches on the tree."

    parser = argparse.ArgumentParser(description=program_description)

    parser.add_argument("-c", "--cluster_list", type=str, help="Cluster file", required=True)
    parser.add_argument("-n", "--align_folder", type=str, help="Alignment folder", required=True)
    parser.add_argument("-t", "--tree_folder", type=str, help="Tree folder", required=True)
    parser.add_argument("-o", "--output_directory", type=str, help="Output folder", required=True)
    parser.add_argument("-p", "--num_processors", type=int, help="Number of processors to use (Default is 1)",
                        default=1)

    args = parser.parse_args()

     #Check for the output folder and also create the temporal folder

    if not os.path.exists(args.output_directory):
        os.makedirs(args.output_directory)

    temp_folder = args.output_directory + "/tmp"
    if not os.path.exists(temp_folder):
        os.makedirs(temp_folder)

    plot_folder = args.output_directory + "/plots"
    if not os.path.exists(plot_folder):
        os.makedirs(plot_folder)

    sites_folder = args.output_directory + "/sites"
    if not os.path.exists(sites_folder):
        os.makedirs(sites_folder)

    #Create output files

    output_file = open(args.output_directory + "/paml_results.txt", 'w')
    no_results_file = open(args.output_directory + "/no_results.txt", 'w')

    #Read the cluster file
    clusters_to_analyze = [line.rsplit()[0] for line in open(args.cluster_list) if line.strip()]

    results_list = []

    #Prepare for multiprocessing

    #Function to store the results
    def store_results(cluster_results):

        for entry in cluster_results:

            node_results, sites_results = entry

            output_file.write("\t".join(str(x) for x in node_results) + "\n")

            cluster_id = node_results[0]
            node = node_results[1]

            if sites_results:
                site_file = open(sites_folder + "/" + cluster_id + "_" + str(node) + ".txt", 'w')

                for position in sites_results:
                    aa, prob = sites_results[position]
                    site_file.write("\t".join(str(x) for x in [position, aa, prob]) + "\n")

                site_file.close()

            results_list.append(node_results)

    #Create the pool of processors
    pool = multiprocessing.Pool(args.num_processors)

    run_results = []

    for cluster in clusters_to_analyze:

        tree_file = args.tree_folder + "/" + cluster + ".tre"
        align_file = args.align_folder + "/" + cluster + ".fna"

        #Check that the files exists
        if not os.path.exists(tree_file):
            print "Tree file missing: " + tree_file
            no_results_file.write(cluster + "\n")
            continue

        if not os.path.exists(align_file):
            print "Alignment missing: " + align_file
            no_results_file.write(cluster + "\n")

        #Check alignment length. If only two sequences, move to the next one
        fasta_count = 0
        for line in open(align_file, 'r'):
            line = line.strip()
            if line.startswith(">"):
                fasta_count += 1

        if not fasta_count > 2:
            continue

        node_id_2_names = defaultdict()

        for entry in EvolTree(tree_file).iter_descendants():
            node_id_2_names[entry.node_id] = entry.get_leaf_names()

        #Results, the first element has:
        #The second is a dictionary with the positive selected sites

        #results_dict[cluster] = run_site_branch(cluster, tree_file, align_file, temp_folder, plot_folder)

        p = pool.apply_async(run_site_branch, args=(cluster, tree_file, align_file, temp_folder, plot_folder,),
                             callback=store_results)

        run_results.append(p)

    pool.close()
    pool.join()

    output_file.close()


