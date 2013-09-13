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



