The code in this repository analyses the data from the paper 
"Evolvability-enhancing mutations in the fitness landscape of an RNA and a protein"
by Andreas Wagner, Nature Communications (2023)

It is executable on a Windows 10 workstation running
spyder 5.3.3., python 3.8.10 (64 bit), and jupyter notebook 6.4.12. It requires
no special installation procedures or software packages except standard python 
modules that are loaded from within the scripts.

The expected runtime for each script on a standard desktop computer (in 2022)
does not exceed one hour.


This directory contains the following scripts

EE_funcs.py: miscellaneous functions to analyze the landscape data. Called
from the other scripts

prot_EE_mutations.py: 

Analyzes the protein (antitoxin) fitness landscape of Lite et al. (eLife 2020).
Requires as an input file the file protein_landscape_data.csv, which is identical to
file GSE153897_Variant_fitness.csv that is published as part of the supplement
to the Lite et al., paper.

Asks for all pairs of wild-type and mutants with two-sided one sample t tests 
whether the difference between the mean fitness of the 
1-neighbors of the mutant and the mean fitness of the 1-neighbors of the 
wild type differs significantly from zero.

Also computes for every wt-mut pair, the fraction of the
neighbors of the wt that are beneficial (have greater fitness than the wild-type),
and the fraction of neighbors of the mutant that have greater fitness than
the mutant, as well as the average benefits of these beneficial mutations. 

The output of this script is contained in file protein_EEmutations_out.txt


RNA_EE_mutations.py

Analyzes the yeast tRNA fitness landscape of Domingo et al. (Nature 2018).
Requires as an input file the file RNA_landscape_data.txt, which is identical to
Table S1 from this paper, but modified to 
eliminate the table caption, and to preserve only the following columns: 
seq, num_vars. pos_vars, ntd_vars, ntd_wt, fitness, SE.

Asks for all pairs of wild-type and mutants
with two-sided one sample t tests whether the difference 
between the mean fitness of the 
1-neighbors of the mutant and the mean fitness of the 1-neighbors of the wild type
differs significantly from zero.

Also computes for every wt-mut pair, the fraction of the
neighbors of the wt that are beneficial (have greater fitness than the wild-type),
and the fraction of neighbors of the mutant that have greater fitness than
the mutant, as well as the average benefits of these beneficial mutations.

The output of this script is contained in file RNA_EEmutations_out.txt

adaptive_walks_protein.ipynb (Jupyter notebook) 

builds an adaptive landscape graph from the (preprocessed) fitness data 
in file prot_EEmutations_out.txt, performs adaptive walk simulations 
on this graph, as well as some statistical analyses of the random walks.

The scripts writes the graph in gml format as protein_gnet.gml. 
Other output is displayed within the script. 

adaptive_walks_RNA.ipynb (Jupyter notebook) 

builds an adaptive landscape graph from the (preprocessed) fitness data 
in file RNA_EEmutations_out.txt, performs adaptive walk simulations 
on this graph, as well as some statistical analyses of the random walks.

The scripts writes the graph in gml format as RNA_gnet.gml. 
Other output is displayed within the script. 

