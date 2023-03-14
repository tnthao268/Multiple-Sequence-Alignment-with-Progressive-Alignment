# Report on the alignment example of micro RNA 126 (mir126) sequences 

Progressive alignment using this tool has been conducted with five mir126 sequences from different species from following FASTA files (saved in ```data``` folder):

```equus caballus mir126.fasta``` [Link](https://www.ncbi.nlm.nih.gov/nucleotide/NR_033047.1?report=genbank&log$=nucltop&blast_rank=73&RID=10V9NF6N016)

```homo sapien mir126.fasta``` [Link](https://www.ncbi.nlm.nih.gov/nucleotide/NR_029695.1?report=genbank&log$=nucltop&blast_rank=76&RID=10V9NF6N016)

```macaca mulatta mir126.fasta``` [Link](https://www.ncbi.nlm.nih.gov/nucleotide/NR_032398.1?report=genbank&log$=nucltop&blast_rank=74&RID=10V9NF6N016)

```pan troglodytes mir126.fasta``` [Link](https://www.ncbi.nlm.nih.gov/nucleotide/NR_035615.1?report=genbank&log$=nucltop&blast_rank=72&RID=10V9NF6N016)

```sus scrofa mir126.fasta``` [Link](https://www.ncbi.nlm.nih.gov/nucleotide/NR_038594.1?report=genbank&log$=nucltop&blast_rank=71&RID=10V9NF6N016)

and produced the alignment result ```example.txt```


## Discussion

Comparing the produced result to ClustalW alignment result ```data/mega_clustalw.fas``` (using MEGA version 11) [MEGA download link](https://www.megasoftware.net/]) it can be seen that ClustalW is more effective in reducing the open of gaps, which has more biological meaning. 

This can be due to the use of UPGMA algorithm in the project for creating a tree for clustering . UPGMA despites its simplicity in using, assumes a constant substitution rate, over time and phylogenetic lineages [[1]](#1), and has therefore not been very relevant in biological context.



## Source


<a id="1">[1]</a> 
Michael Weiß, Markus Göker,
Chapter 12 - Molecular Phylogenetic Reconstruction,
Editor(s): Cletus P. Kurtzman, Jack W. Fell, Teun Boekhout,
The Yeasts (Fifth Edition),
Elsevier,
2011,
Pages 159-174,
ISBN 9780444521491,
https://doi.org/10.1016/B978-0-444-52149-1.00012-4.
(https://www.sciencedirect.com/science/article/pii/B9780444521491000124)
