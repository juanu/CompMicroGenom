__author__ = 'juan'

def cog_definitions():
    import os
    from collections import defaultdict

    path = os.path.dirname(__file__)

    cog_one_letter = defaultdict(str)
    desc_cog_number = defaultdict(str)

    for line in open(path + "/cog_list.txt", 'r'):
        if line.strip():
            line = line.rstrip()
            letter_code, cog_number, cog_def = line.split("\t")
            letter_code = letter_code.replace("[", '')
            letter_code = letter_code.replace("]", '')

            cog_one_letter[cog_number] = letter_code
            desc_cog_number[cog_number] = cog_def

    desc_cog_letter = {
        "J": "Translation, ribosomal structure and biogenesis",
        "A": "RNA processing and modification",
        "K": "Transcription",
        "L": "Replication, recombination and repair",
        "B": "Chromatin structure and dynamics",
        "D": "Cell cycle control, cell division, chromosome partitioning",
        "Y": "Nuclear structure",
        "V": "Defense mechanisms",
        "T": "Signal transduction mechanisms",
        "M": "Cell wall/membrane/envelope biogenesis",
        "N": "Cell motility",
        "Z": "Cytoskeleton",
        "W": "Extracellular structures",
        "U": "Intracellular trafficking, secretion and vesicular transport",
        "O": "Posttranslational modification, protein turnover, chaperons",
        "C": "Energy production and conversion",
        "G": "Carbohydrate transport and metabolism",
        "E": "Amino acid transport and metabolism",
        "F": "Nucleotide transport and metabolism",
        "H": "coenzyme transport and metabolism",
        "I": "Lipid transport and metabolism",
        "P": "Inorganic ion transport and metabolism",
        "Q": "Secondary metabolites biosynthesis, transport and catabolism",
        "R": "General function prediction only",
        "S": "Function unknown"
    }

    return cog_one_letter, desc_cog_letter, desc_cog_number




