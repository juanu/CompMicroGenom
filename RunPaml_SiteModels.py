def cluster_analysis(cluster_list, cluster_folder, output_folder, temporal_folder, results, not_found):
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
        paml_sites_results = paml_run.paml_sites(new_alignment_file, new_tree,
                                                      output_folder, temporal_folder)

        #Calculate pvalue

        pvalue_m1_m2 = paml_stats.lrt(paml_sites_results[1].get("lnL"), paml_sites_results[2].get("lnL"), 2)
        pvalue_m7_m8 = paml_stats.lrt(paml_sites_results[7].get("lnL"), paml_sites_results[8].get("lnL"), 2)

        #Store the omega and proportion of sites,based on the M8 model
        proportion_sites = float(paml_sites_results[8]["site_classes"][10]["proportion"])
        omega_value = float(paml_sites_results[8]["site_classes"][10]["omega"])

        #Store final results

        results.append([cluster, number_sequences, alignment_length, round(pvalue_m1_m2, 3), round(pvalue_m7_m8, 3),
                        proportion_sites, omega_value])


if __name__ == '__main__':
    import os
    import argparse
    import multiprocessing
    from multiprocessing import Manager
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

    #Read the cluster file and group file
    clusters_to_analyze = [line.rsplit()[0] for line in open(args.cluster_list) if line.strip()]

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

        p = multiprocessing.Process(target=cluster_analysis,
                                    args=(chunk, args.cluster_folder, args.output_directory,
                                          temporal_folder, cluster_paml_results, groups_no_data))

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
        no_results_file.write(entry + "\n")

    output_file.close()
    no_results_file.close()

    #Run False discovery rate analysis, if chosen. I need to correct both, the M1a vs M2 pvalues, and then
    #M7 vs M8
    if args.fdr:
        #Correct M1-M2
        corrected_pvalue_results = paml_stats.fdr(cluster_paml_results, 3)

        #Correct M7-M8
        corrected_pvalue_results = paml_stats.fdr(corrected_pvalue_results, 4)

        corrected_results_file = open(args.output_directory + "/paml_results_corrected.txt", "w")

        for result in corrected_pvalue_results:
            corrected_results_file.write("\t".join(str(x) for x in result) + "\n")


