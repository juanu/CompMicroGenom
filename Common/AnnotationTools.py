def parse_annotation_folder(genome_jgi_list, annotation_folder):
    """
    Read the annotation for each genome. There are five possible options PFAM, Product, Feature_type
    COG, KO
    """

    from collections import defaultdict

    genome_annotation = defaultdict(lambda: defaultdict())
    function_defs = defaultdict()

    for genome in genome_jgi_list:
        #This are the genomes preprocessed, where the extension is txt
        genome_file = annotation_folder + "/" + genome + ".txt"

        for line in open(genome_file, 'r'):
            if line.strip():
                line = line.rstrip()

                info = line.split("\t")

                genome_annotation[info[0]][info[1]] = info[2]

                try:
                    function_defs[info[2]] = info[3]

                except IndexError:
                    continue

    return genome_annotation, function_defs