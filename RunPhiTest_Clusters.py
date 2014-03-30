#!/usr/local/bin/python
#Created on 8/1/13

__author__ = 'Juan A. Ugalde'


#Over 0.05 to select. When is less than tat (p-value 0.05) we reject the null hypotehsis (no recombintaiton) and accept
#the hypothesis that there is recombination.
#In my case I'm interested in the recombination free groups, so the p-value has to be more than 0.05

def rename_phi_outfiles(folder, cluster):

    try:
        os.chdir(os.path.abspath(folder))

        phi_outfiles = ["Phi.poly.unambig.sites", "Phi.log", "Phi.inf.sites", "Phi.inf.list"]

        for filename in phi_outfiles:
            new_filename = cluster + "." + filename
            os.rename(filename, new_filename)

    except OSError:
        pass


def process_phi_log(cluster):
    import re

    pvalues = dict()

    log_filename = cluster + "." + "Phi.log"

    for line in open(log_filename, 'r'):
        if line.strip():
            line = line.rstrip()
            if line.startswith("NSS"):
                results = re.split("\s+", line)

                pvalues["NSS"] = results[1]

            elif line.startswith("Max"):
                results = re.split("\s+", line)
                pvalues["MaxChi2"] = results[2]

            elif line.startswith("PHI (Permutation)"):
                results = re.split("\s+", line)
                pvalues["PHI"] = results[2]

    return pvalues

if __name__ == '__main__':
    import argparse
    import os
    from collections import defaultdict

    program_description = "This script takes a folder with aligned clusters (either protein or DNA), " \
                          "and runs PhiTest"

    parser = argparse.ArgumentParser(description=program_description)

    parser.add_argument("-c", "--aligned_clusters_folder", type=str, help="Cluster file", required=True)
    #parser.add_argument("-t", "--sequence_type", type=str, help="Type of sequence. Nucleotide by default", required=True)
    parser.add_argument("-p", "--phitest_location", type=str, help="Location of phitest", required=True)
    parser.add_argument("-o", "--output_directory", type=str, help="Output folder", required=True)

    args = parser.parse_args()

    #Set the paths to the folders
    output_folder = os.path.abspath(args.output_directory)
    cluster_folder = os.path.abspath(args.aligned_clusters_folder)
    phitest_program = os.path.abspath(args.phitest_location) + "/Phi"

    #Create the output folders
    phitest_results_folder = output_folder + "/phitest_results"

    if not os.path.exists(output_folder):
        os.makedirs(args.output_directory)

    if not os.path.exists(phitest_results_folder):
        os.makedirs(phitest_results_folder)

    #Get the clusters
    cluster_list = os.listdir(cluster_folder)

    #Set the sequence type
    #I'm using nucleotide by default right now

    #Move to the phitest folder
    os.chdir(phitest_results_folder)

    #Run Phitest

    clusters_results = defaultdict(list)

    processed_clusters = 0

    for cluster_name in cluster_list:
        cluster_file = cluster_folder + "/" + cluster_name

        os.system("%s -f %s -p 1000 -w 50 -o" % (phitest_program, cluster_file))

        rename_phi_outfiles(phitest_results_folder, cluster_name[:-4])


        #Process the log file

        pvalues = process_phi_log(cluster_name[:-4])

        #I want that all the test validate the null hypothesis

        test_count = 0

        for test in pvalues:
            if float(pvalues[test]) >= 0.05:
                test_count += 1

        if test_count == 3:
            clusters_results["no_recombination"].append(cluster_name[:-4])
        else:
            clusters_results["recombination"].append(cluster_name[:-4])

        processed_clusters += 1

    #Print the results
    logfile = open(output_folder + "/logfile.txt", 'w')

    logfile.write("Processed clusters: %d \n" % processed_clusters)
    logfile.write("Clusters with evidence of recombination: %d\n" % len(clusters_results["recombination"]))
    logfile.write("Clusters with no recombination: %d\n" % len(clusters_results["no_recombination"]))

    #Print a list of each cluster
    clusters_no_recombination = open(output_folder + "/list_clusters_no_recombination.txt", 'w')
    clusters_recombination = open(output_folder + "/list_clusters_recombination.txt", 'w')

    for cluster in clusters_results["no_recombination"]:
        clusters_no_recombination.write(cluster + "\n")

    for cluster in clusters_results["recombination"]:
        clusters_recombination.write(cluster + "\n")

    logfile.close()
    clusters_no_recombination.close()
    clusters_recombination.close()

    #Go back to the output folder






