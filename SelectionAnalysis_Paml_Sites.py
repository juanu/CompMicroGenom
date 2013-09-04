
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