#Comparative genome analysis with OrthoMCL

Author: _Juan A. Ugalde_, juanuu@gmail.com

##Some considerations
- You need to have OrthoMCL installed on your machine, and the right permissions to create mysql databases


##Preparing the genome files

### Download the data from IMG

You need to download three files for each genome, and store them in separate folders, according to data type:
1. Annotation: the format is *GenomeNumber*.info.xls
2. Fasta proteins: the format is *GenomeNumber*.faa
3. Fasta genes (nucleotide): the format is *GenomeNumber*.fna

An important thing to keep in mind, is that sometimes when you are downloading a file from the IMG website, the file does not have the correct name and instead the name is main.cgi. Is that is the case, you need to appropiate rename the file.

Another detail, is that the extension of each file does not need to be faa or fna, but all the files in the folder must have the same extension.

### Create a genome list file

This will be used to replace the JGI names with appropiate (and more easy to remember) names. The format of this file is:

Taxon OID (which is the prefix of each file)  Full genome name  Prefix to use

where each file is separated by tabs

An example of this:
645058727       Acidithiobacillus caldus ATCC 51756     Acaldus51756

### Prepare the files
Use the __PrepareGenomes_OrthoMCL.py__ script. The required inputs for the script are:
- Genome list file
- Folder with either nucleotide or protein sequences
- Extension of the files in the folder
- Name of the output folder (it will be created if it doesn't exist)

### Run Blast
- Concatenate all of the genomes files into one single file
- Create blast database
- Blast all the sequences versus the blast database, the output must be in tabular format

Example command:

formatdb -i All_genomes.fasta -n AcidoProteins -p T

blastall -p blastp -d AcidoProteins -i All_genomes.fasta -F 'm S' -v 1000000 -b 1000000 -z 19055 -e 1e-5 -m 8 -a 6 -o blastp.Acidothio

### Create MySQL database and run MCL
1. Create MySQL database, and give privileges to the user running the analysis

2. Create the config file (example provided in the orthoMCL folder), and edit it with the information of the mysql database, user name and password.

3. Install orthoMCL schema:
orthomclInstallSchema __config_file__

4. Parse Blast Results:
orthomclBlastParser __blast_input__ __fasta_protein_folder__ >> similarSequences.txt

5. Load results in the database
orthomclLoadBlast __config_file__ similarSequences.txt

6. Run MCL pairs
orthomclPairs __config_file__ __log_file__ cleanup=no

7. Dump the pairs files
orthomclDumpPairsFiles __config_file__

8. MCL analysis
It is possible here to try with different inflation values (I)
mcl mclInput --abc -I 1.5 -o 1.5-mclOutput

9. Get group names for the MCL clusters
10. orthomclMclToGroups __group_name__ __number__ < __mcl_output_file__ > __group_output_file__

###Analyze the orthoMCL output with the custom python scripts
1. PrepareGenome_OrthoMCL.py
2. AnnotateOrthoMCL_Clusters.py
3. ClusterAlignmentTree.py





