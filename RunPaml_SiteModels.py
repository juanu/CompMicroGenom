def single_cluster_analysis(cluster_id, cluster_folder, output_folder, temp_folder, outfile_notfound):
    """
    This function take the group, alignment, tree and folder information and runs a paml analysis using the M1a, M2a,
    M7 and M8 models
    The working dir is important (different from the output dir), because different PAML runs at the same time may
    override each other.
    This is particularly important if running this script in more than one processor
    """

    from SelectionAnalysis import paml_stats
    from SelectionAnalysis import paml_prepare
    from SelectionAnalysis import paml_run

    cluster_file = cluster_folder + "/" + cluster_id + ".fna"  # Add fna extension

    #Check that the cluster file exists, if not continue
    if not os.path.exists(cluster_file):
        outfile_notfound.write(cluster_id + "\n")

        #Make a new tree, no confidence values in the branches
    new_tree = paml_prepare.run_fasttree(cluster_file, temp_folder)

        #Make the new alignment, and get information about the alignment
    new_alignment_file, number_sequences, alignment_length = \
        paml_prepare.adjust_alignment(cluster_file, temp_folder)

    #Run PAML for each branch in the cluster with both models
    paml_sites_results = paml_run.paml_sites(new_alignment_file, new_tree, output_folder, temp_folder)

    #Calculate pvalue

    pvalue_m1_m2 = paml_stats.lrt(paml_sites_results[1].get("lnL"), paml_sites_results[2].get("lnL"), 2)
    pvalue_m7_m8 = paml_stats.lrt(paml_sites_results[7].get("lnL"), paml_sites_results[8].get("lnL"), 2)

    #Store the omega and proportion of sites,based on the M8 model
    try:
        proportion_sites = float(paml_sites_results[8]["site_classes"][10]["proportion"])
        omega_value = float(paml_sites_results[8]["site_classes"][10]["omega"])
    except TypeError:
        proportion_sites = 0
        omega_value = 0

    #Store final results

    summary_results = [cluster_id, number_sequences, alignment_length, round(pvalue_m1_m2, 3), round(pvalue_m7_m8, 3),
                        proportion_sites, omega_value]

    print summary_results



    return summary_results

if __name__ == '__main__':
    import os
    import argparse
    import multiprocessing
    from SelectionAnalysis import paml_stats

    program_description = "Script that takes a list of clusters and runs PAML (codeml), using a site approach." \
                          "The models used are M1a vs M2a and M7 vs M8. "

    parser = argparse.ArgumentParser(description=program_description)

    parser.add_argument("-c", "--cluster_list", type=str, help="Cluster file", required=True)
    parser.add_argument("-n", "--cluster_folder", type=str, help="Output folder", required=True)
    parser.add_argument("-o", "--output_directory", type=str, help="Output folder", required=True)
    parser.add_argument("-p", "--num_processors", type=int, help="Number of processors to use", required=True)
    parser.add_argument("-f", "--fdr", help="Perform false discovery rate")

    args = parser.parse_args()

    #Check for the output folder and also create the temporal folder

    if not os.path.exists(args.output_directory):
        os.makedirs(args.output_directory)

    #Output files

    output_file = open(args.output_directory + "/paml_results.txt", 'w')
    no_results_file = open(args.output_directory + "/no_results.txt", 'w')

    #Read the cluster file and group file
    clusters_to_analyze = [line.rsplit()[0] for line in open(args.cluster_list) if line.strip()]

    #Added on 3/29/14, to fix problems running processes in parallel

    pool = multiprocessing.Pool(args.num_processors)

    cluster_paml_results = []
    groups_no_data = []
    i = 1


    def store_results(result):
        output_file.write("\t".join(str(x) for x in result) + "\n")
        cluster_paml_results.append(result)

    run_results = []


    for cluster in clusters_to_analyze:
        temporal_folder = args.output_directory + "/temp_" + str(i)

        if not os.path.exists(temporal_folder):
            os.makedirs(temporal_folder)

        i += 1

        r = pool.apply_async(single_cluster_analysis, args=(cluster, args.cluster_folder,
                                                        args.output_directory, temporal_folder, no_results_file),
                                                        callback=store_results)

        run_results.append(r)

    pool.close()
    pool.join()

    #for result in cluster_paml_results:
    #    output_file.write("\t".join(str(x) for x in result) + "\n")

    #for entry in groups_no_data:
    #    no_results_file.write(entry + "\n")

    output_file.close()
    no_results_file.close()

    #Run False discovery rate analysis, if chosen. I need to correct both, the M1a vs M2 pvalues, and then
    #M7 vs M8
    if args.fdr:
        #Correct M1-M2
        corrected_pvalue_results = paml_stats.fdr(cluster_paml_results, 3)

        #Correct M7-M8
        corrected_pvalue_results = paml_stats.fdr(corrected_pvalue_results, 5)

        corrected_results_file = open(args.output_directory + "/paml_results_corrected.txt", "w")

        for result in corrected_pvalue_results:
            corrected_results_file.write("\t".join(str(x) for x in result) + "\n")
