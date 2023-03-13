# Multiple-Sequence-Alignment with Progressive-Alignment
This projects contains multiple sequence alignments tool written in Julia language using Progressive Alignment 

It builds phylogenetic tree from sequences with UPGMA algorithm, then align these sequences according to the tree

Usage examples are in ```Project_example``` file, tests are in ```test``` folder

## Source code files and guideline to their usages

1. **Data Reader** 

To read data from Fasta files and save the information into Record datatype 

Call **readSequences_file** method to read many fasta files. **input_files_names** is a list of names from all fasta files: 

```
input_files_names = ["data/ChimpanzeeDND1.fasta","data/DogDND1.fasta","data/HomoSapienDND1.fasta","data/MouseDND1.fasta"] 

records_sequences = readSequences_file(input_files_names) 
 ```

**records_sequences** in this example is a list of Record objects. Record can be also created directly from description and DNA sequence string. 

```
fasta_record  = Record("test string1", "ACGT") 

```

Call **readSequence** method to read a fasta file (sequence.fasta): 


```
sequence = readSequences(dirname(@__FILE__()) * "/..data/sequence.fasta") 

```

2. **Distance Matrix and Dictionary** 

Use **records_sequences** as parameter to create a **Distance Matrix** and a **Dictionary**. 


```
distanceMatrix = createDistanceMatrix(records) 

dict = createDictionary(records) 

```

**dict** is a NamedTuple, which has **dict_records** and **leave_names**. **leave_names** is a list of letters in the alphabet. Each letter presents a **Record** from **records_sequences**. And they are described as keys and values of dictionary **dict_records**. 

```

dict_records = dict.dict_records 

leaf_names = dict.leaf_names 

```
3. **Guild Tree (UPGMA)**

Here **UPGMA** is selected as a tree-based algorithm to perform multiple alignment. With the help of this tree, alignment is performed based on cluster of sequences. Sequences which are firstly clustered are firstly aligned. 

Call **split_name_sequences** method to return list of clusters that need to be aligned: 

```
guildTree_instruction = split_name_sequences(cluster_list(change_cluster_name(tree))) 
```

**tree** is the created tree from UPGMA algorithm, using distance matrix and list of names. 

```
tree = upgma(distanceMatrix,leaf_names)
```
4. **Progressive Alignment**

**guildTree_instruction** shows pairs, which will be aligned, in order. A pair can be 2 sequences (pairwise alignment), or sequence-alignment or 2 alignments (multi sequences alignment).  

To run the Progressive Alignment call the method: 
```

result = progressiveAlignment(guildTree_instruction, dict_records) 
```

The result returns a dictionary of Records, which contain all aligned sequences. 

5. **DataWriter**

The aligned sequences are written in a text file format 

Call **writeSequences** method with output file name ("try_seq1.txt") and dictionary containing representation string of the records (dict_seq). *50* in the below example is the sequence's number of characters printed on each line. 

```
dict_seq = Dict("A" => record, "B" => record2) 
writeSequences("try_seq1.txt", dict_seq,50) 
```

## Library

BioAlignment [GitHub](https://github.com/BioJulia/BioAlignments.jl.git)

FastaIO [Link](https://docs.juliahub.com/FastaIO/i12XQ/1.0.0/)

MCPhyloTree [GitHub](https://github.com/erathorn/MCPhyloTree.jl.git)

## Language

Julia Version 1.8.2 [Link](https://julialang.org/downloads/)


