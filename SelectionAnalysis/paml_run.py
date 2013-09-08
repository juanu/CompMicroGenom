__author__ = 'juan'


def paml_sites(alignment, tree, output_dir, working_dir):
    """
    This
    """
    from Bio.Phylo.PAML import codeml
    import os

    paml_results = dict()  # Store the results of the analysis

    cml = codeml.Codeml()  # Setup PAML

    #Parameters to PAML
    cml.alignment = alignment
    cml.tree = tree
    cml.out_file = output_dir + "/" + os.path.basename(alignment)[:-4] + os.path.basename(tree)[:-4] + ".sites.m1a"
    cml.working_dir = working_dir

    cml.set_options(seqtype=1, CodonFreq=2, clock=0, model=0, NSsites=[1, 2, 7, 8],  fix_kappa=0, kappa=2,
                    fix_omega=0, omega=5, verbose=1, fix_blength=1)

    print "Running codeml for models 1, 2, 7, 8 in : %s" % os.path.basename(tree)

    results_ma = cml.run()

    #Parse the results for the first run
    ns_sites_ma = results_ma.get("NSsites")

    for site in ns_sites_ma:
        lnL = ns_sites_ma[site].get("lnL")
        parameters = ns_sites_ma[site].get("parameters")
        site_classes = parameters.get("site classes")

        model_results = {"lnL": lnL, "site_classes": site_classes}

        paml_results[site] = model_results

    return paml_results


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