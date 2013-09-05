#!/usr/local/bin/python
#Created on 8/20/13

__author__ = 'Juan Ugalde'


def LRT_paml(H0, H1, df):
    """
    Perform the LRT test for the PAML results. I'm considering two calcuation, with 2 degrees of freedom (M1a vs M2a)
    and 4 degrees (M0 vs M3)
    """

    import math
    from scipy.stats import chi2

    LRT_value = 2 * math.fabs(H0 - H1)
    pval = 1 - chi2.cdf(LRT_value, df)

    return pval


def run_paml_ma_m1a(alignment, tree, output_dir, working_dir):
    """
    This is tu run PAML in each defined branch under models MA and M1a
    """
    from Bio.Phylo.PAML import codeml
    import os

    paml_results = dict()

    cml = codeml.Codeml()

    cml.alignment = alignment
    cml.tree = tree
    cml.out_file = output_dir + "/" + os.path.basename(tree)[:-4] + ".ma"
    cml.working_dir = working_dir

    cml.set_options(seqtype=1, CodonFreq=2, clock=0, model=2, NSsites=[2],  fix_kappa=0, kappa=2,
                    fix_omega=0, omega=5, verbose=1, fix_blength=1)

    print "Running codeml for model A in : %s" % os.path.basename(tree)

    results_ma = cml.run()

    #From here I need some of the results:

    ns_sites_ma = results_ma.get("NSsites")

    for site in ns_sites_ma:
        lnL = ns_sites_ma[site].get("lnL")
        parameters = ns_sites_ma[site].get("parameters")
        site_classes = parameters.get("site classes")

        model_results = {"lnL": lnL, "site_classes": site_classes}

        paml_results["Ma"] = model_results

    print "Running codeml for model 1A in : %s" % os.path.basename(tree)

    cml.out_file = output_dir + "/" + os.path.basename(tree)[:-4] + ".m1a"
    cml.set_options(seqtype=1, CodonFreq=2, clock=0, model=2, NSsites=[2],  fix_kappa=0, kappa=2,
                    fix_omega=1, omega=1, verbose=1, fix_blength=1)

    results_m1a = cml.run()

    ns_sites_m1a = results_m1a.get("NSsites")

    for site in ns_sites_m1a:
        lnL = ns_sites_m1a[site].get("lnL")
        parameters = ns_sites_m1a[site].get("parameters")
        site_classes = parameters.get("site classes")

        model_results = {"lnL": lnL, "site_classes": site_classes}

        paml_results["M1a"] = model_results

    return paml_results


def run_paml_per_group(groups, alignment, tree, output_dir, working_dir):
    """
    This is to take each defined group, modify the tree, and the run PAML
    """
    from Bio import Phylo

    cluster_tree = Phylo.read(tree, "newick")

    #Names have a pipe sign (|) with the organism|protein_id.
    # I need to keep track of everything to replace in the final tree
    #clades_in_tree = {str(clade).split("|")[0]: str(clade).split("|")[1] for clade in cluster_tree.get_terminals()}

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

            old_tree_informtion = open(tree).read()

            new_tree_information = multiple_replace(dict_new_clade_names, old_tree_informtion)

            group_tree = working_dir + "/" + group + ".tre"

            new_tree_file = open(group_tree, 'w')
            new_tree_file.write(new_tree_information)
            new_tree_file.close()

            #Run model for the new tree

            paml_results = run_paml_ma_m1a(alignment, group_tree, output_dir, working_dir)

            clade_results[group] = paml_results

        else:
            clade_results[group] = None

    return clade_results


def adjust_alignment(alignment_file, output_folder):
    """
    This adjust the alignment into an appropiate format for PAML. The format used is like this:
      NumberSeqs   Length Alignment
    ID
    SEQ
    ID
    SEQ
    """
    from Bio import AlignIO
    import os

    #Files are DNA, so always use extension fna
    alignment = AlignIO.read(open(alignment_file), "fasta")

    new_alignment_file = output_folder + "/" + os.path.basename(alignment_file[:-4]) + ".paml"

    output_alignment = open(new_alignment_file, 'w')

    output_alignment.write(" %d  %d\n" % (len(alignment), alignment.get_alignment_length()))

    for record in alignment:
        output_alignment.write("%s\n" % record.id)
        output_alignment.write("%s\n" % record.seq)

    return new_alignment_file, len(alignment), alignment.get_alignment_length()


def run_fasttree(alignment_file, output_folder):
    """
    Here we run FastTree with this options:
    -slow
    -gtr
    -nosupport

    """
    import os

    print "Making tree %s" % alignment_file

    tree_file = output_folder + "/" + os.path.basename(alignment_file)[:-4] + ".tre"

    os.system("FastTree -slow -nt -gtr -nosupport -quiet %s > %s" % (alignment_file, tree_file))

    return tree_file


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


if __name__ == '__main__':
    import os
    import argparse
    from collections import defaultdict
    from operator import itemgetter

    program_description = "Script that takes a list of clusters and runs PAML (codeml). The model used is a branch-site" \
                          "with relaxed test (MA vs M1a). "

    parser = argparse.ArgumentParser(description=program_description)

    parser.add_argument("-c", "--cluster_list", type=str, help="Cluster file", required=True)
    parser.add_argument("-n", "--cluster_folder", type=str, help="Output folder", required=True)
    parser.add_argument("-g", "--groups", type=str, help="Group constrains", required=True)
    parser.add_argument("-o", "--output_directory", type=str, help="Output folder", required=True)

    args = parser.parse_args()

    #Check for the output folder and also create the temporal folder
    #I'm using the PID to create the temporary folder, which should allow multiple instances of the script to run

    temporal_folder = args.output_directory + "/temp_" + str(os.getpid())

    if not os.path.exists(args.output_directory):
        os.makedirs(args.output_directory)

    if not os.path.exists(temporal_folder):
        os.makedirs(temporal_folder)

    #Read the cluster file and group file

    clusters_to_analyze = [line.rsplit()[0] for line in open(args.cluster_list) if line.strip()]

    group_constrains = defaultdict(list)

    for line in open(args.groups):
        if line.strip():
            line = line.rstrip()
            group_constrains[line.split("\t")[0]].append(line.split("\t")[1])

    #Result and output files

    cluster_paml_results = list()
    groups_no_data = list()

    for cluster in clusters_to_analyze:

        cluster_file = args.cluster_folder + "/" + cluster + ".fna"  # Add fna extension

        new_tree = run_fasttree(cluster_file, temporal_folder)  # make tree

        new_alignment_file, number_sequences, alignment_length = \
            adjust_alignment(cluster_file, temporal_folder)  # convert alignment to right format

        #m0_results = run_paml_m0(new_alignment_file, new_tree, args.output_directory, temporal_folder) # Run M0

        paml_site_branch_results = run_paml_per_group(group_constrains, new_alignment_file, new_tree,
                                                               args.output_directory, temporal_folder)

       # group_results = defaultdict(list)

        for group in paml_site_branch_results:

            if paml_site_branch_results[group] is None:
                groups_no_data.append([cluster, group])

            else:
                pvalue = LRT_paml(paml_site_branch_results[group]["Ma"].get("lnL"),
                                  paml_site_branch_results[group]["M1a"].get("lnL"), 1)

                qvalue = 0

                proportion_sites = float(paml_site_branch_results[group]["Ma"]["site_classes"][2]["proportion"]) + \
                                   float(paml_site_branch_results[group]["Ma"]["site_classes"][3]["proportion"])

                average_omega = (float(paml_site_branch_results[group]["Ma"]["site_classes"][2]["branch types"]["foreground"]) +
                                 float(paml_site_branch_results[group]["Ma"]["site_classes"][3]["branch types"]["foreground"])) / 2

                #Store the final results
                #Group, Nseqs, Length, p-value, q-value, P1 in Ma, Omega in W
                cluster_paml_results.append([cluster, group, number_sequences, alignment_length,
                                             round(pvalue, 3), qvalue, proportion_sites, average_omega ])


    #Perform FDR, based on Benjamini approach, at 5%. Add this to the qvalue in the dictionary
    #Based on this post:
    #http://stats.stackexchange.com/questions/870/multiple-hypothesis-testing-correction-with-benjamini-hochberg-p-values-or-q-va

    total_tests = len(cluster_paml_results)  # Total number of performed tests
    position = 1
    prev_adjusted_pvalue = 0

    print cluster_paml_results
    print total_tests

    sorted_results = sorted(cluster_paml_results, key=itemgetter(4))


    for entry in sorted_results:
        adjusted_pvalue = float(entry[4]) * (total_tests / position)

        print adjusted_pvalue

        #If the value is greater than 1, we set as one (0 < p < 1)
        adjusted_pvalue = min(adjusted_pvalue, 1)

        #Check that the value is not greater than the previous one
        adjusted_pvalue = max(adjusted_pvalue, prev_adjusted_pvalue)

        prev_adjusted_pvalue = adjusted_pvalue
        position += 1

        entry[5] = adjusted_pvalue



    print cluster_paml_results

    #print cluster_paml_results
    #Summary information






