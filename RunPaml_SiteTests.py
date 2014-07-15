#Created on 7/14/14
__author__ = 'Juan Ugalde'

#Created on 7/11/2014
__author__ = 'Juan Ugalde'

#TODO
#Clean and document the script
#Check when files are not available

def run_site_tests(cluster_name, treefile, alignment, folder_temp, folder_plots):
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

    #Run M1 as the null model
    tree.run_model("M1")

    #Run M2 as the alternative model
    tree.run_model("M2")
    model1 = tree.get_evol_model("M1")
    model2 = tree.get_evol_model("M2")  # Get the results of the model

    #Run the LRT test, using ETE
    #pval = tree.get_most_likely("M2", "M1")

    #Get the positive selected sites
    ps_sites = defaultdict()
    total_sites = 0
    sites_over_95 = 0

    #Output file
    output_list = []

    for s in range(len(model2.sites['BEB']['aa'])):
        p_value_site = float(model2.sites['BEB']['p2'][s])

        if p_value_site > 0.50:
            ps_sites[s] = [model2.sites['BEB']['aa'][s], model2.sites['BEB']['p2'][s]]
            total_sites += 1

            if p_value_site > 0.95:
                sites_over_95 += 1


    #LRT Test

    lrt_value = 2 * math.fabs(model1.lnL - model2.lnL)  # LRT test value
    pval = 1 - chi2.cdf(lrt_value, 2)  # p-value based on chi-square

    test_status = None

    #Evidence of positive selection in the branch
    omega_value = float(model2.classes['w'][2])
    proportion_sites = float(model2.classes['proportions'][2])

    #Plot file
    plot_file = folder_plots + "/" + cluster_name

    if pval < 0.05 and omega_value > 1:
        #Save plots, both in jpg and svg of the clusters with evidence of positive selection
        test_status = "Positive"
        tree.render(plot_file + ".svg", layout=evol_clean_layout)
        tree.render(plot_file + ".jpg", layout=evol_clean_layout)
    else:
         #print "no signal"
        test_status = None


    result_entry = [cluster_name, omega_value, proportion_sites, pval, test_status,
                        total_sites, sites_over_95]

       # print result_entry
        #print ps_sites
        #node_results[node.node_id] = [result_entry, ps_sites]
    output_list = [result_entry, ps_sites]

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
                          "site-branch test on them. It can also perform a FDR analysis using the XXX approach." \
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
    def store_results(combined_results):

        entry_results, sites_results = combined_results

        output_file.write("\t".join(str(x) for x in entry_results) + "\n")

        cluster_id = entry_results[0]

        if not sites_results:

            site_file = open(sites_folder + "/" + cluster_id + ".txt", 'w')

            for position in sites_results:
                aa, prob = sites_results[position]
                site_file.write("\t".join(str(x) for x in [position, aa, prob]) + "\n")

            site_file.close()

        results_list.append(entry_results)


    #Create the pool of processors
    pool = multiprocessing.Pool(args.num_processors)

    run_results = []

    for cluster in clusters_to_analyze:

        tree_file = args.tree_folder + "/" + cluster + ".tre"
        align_file = args.align_folder + "/" + cluster + ".fna"

        node_id_2_names = defaultdict()

        for entry in EvolTree(tree_file).iter_descendants():
            node_id_2_names[entry.node_id] = entry.get_leaf_names()

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

        #Results, the first element has:
        #The second is a dictionary with the positive selected sites

        #results_dict[cluster] = run_site_branch(cluster, tree_file, align_file, temp_folder, plot_folder)

        p = pool.apply_async(run_site_tests, args=(cluster, tree_file, align_file, temp_folder, plot_folder,),
                             callback=store_results)

        run_results.append(p)

    pool.close()
    pool.join()

    #output_file.close()


