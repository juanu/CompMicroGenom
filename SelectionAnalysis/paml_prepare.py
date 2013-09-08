__author__ = 'juan'


def adjust_alignment(alignment_file, output_folder):
    """
    This function is used to transform the fasta DNA alignment into the right format for
    PAML. The format used is this:

      NumberSeqs   Length Alignment
    ID
    SEQ
    ID
    SEQ
    """
    from Bio import AlignIO
    import os
    import sys

    #Files are DNA, so always use the extension fna. Check for this
    if not alignment_file.endswith(".fna"):
        sys.exit("Alignment does not have the fna extension")

    alignment = AlignIO.read(open(alignment_file), "fasta")  # Open the alignment file

    new_alignment_file = output_folder + "/" + os.path.basename(alignment_file[:-4]) + ".paml"  # New alignment file

    output_alignment = open(new_alignment_file, 'w')  # Output file for the new alignment

    #Modify the alignment and write into a new file
    output_alignment.write(" %d  %d\n" % (len(alignment), alignment.get_alignment_length()))

    for record in alignment:
        output_alignment.write("%s\n" % record.id)
        output_alignment.write("%s\n" % record.seq)

    return new_alignment_file, len(alignment), alignment.get_alignment_length()


def run_fasttree(alignment_file, output_folder):
    """
    This function runs FastTree and generates a new tree (based on DNA alignment) without
    support numbers (PAML does not like those numbers in the tree).
    The options used for FastTree are:

    -slow
    -gtr
    -nosupport

    """
    import os
    import subprocess
    import sys

    #Check that FastTree is in the system
    try:
        null = open("/dev/null", 'w')
        subprocess.Popen("FastTree", stdout=null, stderr=null)
        null.close()
    except OSError:
        sys.exit("Check that FastTree is installed and in the PATH")

    # Tell what is doing
    print "Making tree %s" % alignment_file

    #Tree file to be created
    tree_file = output_folder + "/" + os.path.basename(alignment_file)[:-4] + ".tre"

    #Run FastTree
    os.system("FastTree -slow -nt -gtr -nosupport -quiet %s > %s" % (alignment_file, tree_file))

    return tree_file
