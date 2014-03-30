__author__ = 'Juan A. Ugalde'


def ma_m1a(alignment, tree, output_dir, working_dir):
    """
    This is tu run PAML in each defined branch under models MA and M1a, with this options:
    model = 2
    NSsites = 2
    fix_omega = 0 (for Ma) and 1 (for M1a)
    fix_blength = 1 -> The supplied tree should have branch lengths, and PAML will use those as a starting point

    The output of this function is a dictionary containing the lnL value and the site_classes for each model (Ma and M1a)
    """
    from Bio.Phylo.PAML import codeml
    import os

    paml_results = dict()  # Store the results of the analysis

    cml = codeml.Codeml()  # Setup PAML

    #Parameters to PAML
    cml.alignment = alignment
    cml.tree = tree
    cml.out_file = output_dir + "/" + os.path.basename(alignment)[:-4] + os.path.basename(tree)[:-4] + ".ma"
    cml.working_dir = working_dir

    cml.set_options(seqtype=1, CodonFreq=2, clock=0, model=2, NSsites=[2],  fix_kappa=0, kappa=2,
                    fix_omega=0, omega=5, verbose=1, fix_blength=1)

    print "Running codeml for model A in : %s" % os.path.basename(tree)

    results_ma = cml.run()

    #Parse the results for the first run
    ns_sites_ma = results_ma.get("NSsites")

    for site in ns_sites_ma:
        lnL = ns_sites_ma[site].get("lnL")
        parameters = ns_sites_ma[site].get("parameters")
        site_classes = parameters.get("site classes")

        model_results = {"lnL": lnL, "site_classes": site_classes}

        paml_results["Ma"] = model_results

    #Run the second model
    print "Running codeml for model 1A in : %s" % os.path.basename(tree)

    #Parameters for the second model
    cml.out_file = output_dir + "/" + os.path.basename(alignment)[:-4] + os.path.basename(tree)[:-4] + ".m1a"
    cml.set_options(seqtype=1, CodonFreq=2, clock=0, model=2, NSsites=[2],  fix_kappa=0, kappa=2,
                    fix_omega=1, omega=1, verbose=1, fix_blength=1)

    results_m1a = cml.run()

    #Parse the results for the second run
    ns_sites_m1a = results_m1a.get("NSsites")

    for site in ns_sites_m1a:
        lnL = ns_sites_m1a[site].get("lnL")
        parameters = ns_sites_m1a[site].get("parameters")
        site_classes = parameters.get("site classes")

        model_results = {"lnL": lnL, "site_classes": site_classes}

        paml_results["M1a"] = model_results

    return paml_results

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
    from SelectionAnalysis import paml_stats
    from multiprocessing import Pool

    program_description = "Script that takes a list of clusters and runs PAML (codeml). The model used is a branch-site" \
                          "with relaxed test (MA vs MA with omega fixed at 1). "

    parser = argparse.ArgumentParser(description=program_description)

    parser.add_argument("-c", "--cluster_list", type=str, help="List of clusters to analyze. The format is"
                                                               "based on the other scripts of the pipeline", required=True)
    parser.add_argument("-n", "--cluster_folder", type=str, help="Output folder", required=True)
    parser.add_argument("-g", "--groups", type=str, help="Group constrains. Format: Group  Genome")
    parser.add_argument("-o", "--output_directory", type=str, help="Output folder", required=True)
    parser.add_argument("-p", "--num_processors", type=int, help="Number of processors to use (Default is 1)",
                        default=1)
    parser.add_argument("-f", "--fdr", help="Perform false discovery rate")

    args = parser.parse_args()

    #Check for the output folder and also creates the temporal folder

    if not os.path.exists(args.output_directory):
        os.makedirs(args.output_directory)

    #Read the cluster file
    cluster_to_analyze = [line.rstrip()[0] for line in open(args.cluster_list) if line.strip()]

    group_constrains = None

    if args.groups:
        group_constrains = defaultdict(list)
        for line in open(args.groups):
            if line.strip():
                line = line.rstrip()
                group_constrains[line.split("\t")[0]].append(line.split("\t")[1])



    pool = Pool(processes=2)
    #result = pool.apply_async(f, [10])















