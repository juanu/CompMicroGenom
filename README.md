Juan A. Ugalde
juanuu@gmail.com

Created on June 20th, 2013.

This is a collection of scripts used for the comparative analysis of microbial genomes. Starting from a GBK file,
either from Genbank, RAST or IMG, the scripts can be used to prepare data to run orthoMCL, analye the results, generate
alignments, trees and perform selection analysis.

There are also tools for importing dumps from IMG

The main requirements are Python 2.7+ (not Python 3), and libraries:
- PyCogent http://pycogent.org
- Biopython http://biopython.org
- Scipy

In addition, to run the alignment and tree script, Mafft, Fasttree and PAML (codeml) need to be installed and accessible globally.

#####Note (added on November 27th, 2013):
This code is still very preliminary. Our plans are to re-organize the code over the next few months, and have a first release on the first trimester of 2014. Keep an eye on the repository. The code currently here comes with no warranty :) But, any comments, ideas, etc, feel free to email me: juanuu@gmail.com



