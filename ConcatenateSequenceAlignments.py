#!/Users/juan/anaconda/bin/python
#Created on 8/12/2013

__author__ = 'Juan A. Ugalde'

def trim_alignment_gblocks(alignment_file,sequence_type,temporal_folder, logfile):
    """Trim alignment with default parameters using Gblock"""
    import os
    import shutil

    new_file_location = temporal_folder + "/" + os.path.basename(alignment_file)
    shutil.copyfile(alignment_file, new_file_location)

    seq_type = None

    if sequence_type == "nucleotide":
        seq_type = "p"
    else:
        seq_type = "d"

    os.system("Gblocks %s -t=%s >> %s" % (new_file_location, seq_type, logfile))

    new_alignment_file = new_file_location + "-gb"

    return new_alignment_file

if __name__ == '__main__':
    #Take a folder with aligned clusters, and concatenate the alignemnt
    #Check that every cluster has the same number of sequences and the id of the genomes is the same

    import argparse
    from Bio import AlignIO
    import os
    import shutil

    program_description = "This script takes a list of aligned sequences and concatenate those sequences. Each alignments" \
                          "needs to have the same number of sequences in it"

    parser = argparse.ArgumentParser(description=program_description)

    parser.add_argument("-c", "--cluster_list", type=str,
                        help="List of clusters", required=True)
    parser.add_argument("-f", "--cluster_folder", type=str,
                        help="location of the cluster folder", required=True)
    parser.add_argument("-o", "--output_directory", type=str,
                        help="Output folder for the modified fasta files", required=True)
    parser.add_argument("-t", "--trimming", type=str, default="no", choices=['yes', 'no'], help="Performs"
                        "trimming of the alignment using Gblocks")

    args = parser.parse_args()

    if not os.path.exists(args.output_directory):
        os.makedirs(args.output_directory)


    alignment_list = [line.rsplit()[0] for line in open(args.cluster_list) if line.strip()]

    concatenated_alignment = None
    sorted_names = list()
    non_aligned_clusters = list()

    total_alignments = 0
    total_positions_non_trimmed = 0

    #Create temporal folder for gblocks trimming
    temporal_folder = "temp_scratch"
    gblocks_logfile = args.output_directory + "/gblocks_log.txt"
    if args.trimming == "yes":
        if not os.path.exists(temporal_folder):
            os.makedirs(temporal_folder)

    for aln in alignment_list:
        sequence_type = None
        fasta_info = args.cluster_folder + "/" + aln
        fasta_file = None

        if os.path.isfile(fasta_info + ".faa"):
            sequence_type = "aminoacid"
            fasta_file = fasta_info + ".faa"

        elif os.path.isfile(fasta_info + ".fna"):
            sequence_type = "nucleotide"
            fasta_file = fasta_info + ".fna"

        #Add the option of trimming the alignmnet using Gblocks
        if not args.trimming == "no":
            total_positions_non_trimmed += AlignIO.read(open(fasta_file), "fasta").get_alignment_length()
            fasta_file = trim_alignment_gblocks(fasta_file, sequence_type, temporal_folder, gblocks_logfile)

        alignment = AlignIO.read(open(fasta_file), "fasta")

        alignment.sort()

        if concatenated_alignment is None:
            concatenated_alignment = alignment
            for record in alignment:
                genome_name = record.id.split("|")[0]

                sorted_names.append(genome_name)

        else:
            try:
                concatenated_alignment = concatenated_alignment + alignment
            except ValueError:
                non_aligned_clusters.append(aln)
                continue

        total_alignments += 1

    #Print the output alignment
    output_alignment = open(args.output_directory + "/concatenated_alignment.fasta", "w")
    logfile = open(args.output_directory + "/logfile.txt", 'w')
    i = 0

    for sequence in concatenated_alignment:
        output_alignment.write(">%s\n" % sorted_names[i])
        output_alignment.write(str(sequence.seq) + "\n")
        i += 1

    #Logfile
    logfile.write("Total clusters analyzed: %d\n" % total_alignments)
    logfile.write("Clusters not added in the complete alignment: %d\n" % len(non_aligned_clusters))
    logfile.write("Concatenated alignment length: %d\n" % concatenated_alignment.get_alignment_length())

    if args.trimming=="yes":
        logfile.write("Positions without Gblocks trimming: %d\n" % total_positions_non_trimmed)
        shutil.rmtree(temporal_folder)
        #remove temporal folder


    #Non aligned clusters
    failed_clusters = open (args.output_directory + "/failed_clusters.txt", 'w')
    for cluster in non_aligned_clusters:
        failed_clusters.write(cluster + "\n")

    output_alignment.close()
    logfile.close()
    failed_clusters.close()

