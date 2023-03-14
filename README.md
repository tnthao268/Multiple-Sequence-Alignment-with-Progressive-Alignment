# Multiple-Sequence-Alignment with Progressive-Alignment
This projects contains multiple sequence alignments tool written in Julia language using Progressive Alignment 

It builds phylogenetic tree from sequences with UPGMA algorithm, then align these sequences according to the tree

Usage examples are in ```Project_Example``` file, tests are in ```test``` folder

## Source code files and guideline to their usages

### 1. Data Reader

To read data from Fasta files and save the sequence's information (description, sequence)  into Record datatype 

Call ```readSequences_file``` method to read many fasta files. ```input_files_names``` is a list of names from all fasta files: 

```julia
input_files_names = ["data/ChimpanzeeDND1.fasta","data/DogDND1.fasta","data/HomoSapienDND1.fasta","data/MouseDND1.fasta"] 

records_sequences = readSequences_file(input_files_names) 
 ```

```records_sequences``` in this example is a list of Record objects. Record can be also created directly from description and DNA sequence string. 

```julia
fasta_record  = Record("test string1", "ACGT") 

```

Call ```readSequence``` method to read a fasta file (sequence.fasta): 


```julia
sequence = readSequences(dirname(@__FILE__()) * "/..data/sequence.fasta") 

```

### 2. Distance Matrix and Dictionary

Use ```records_sequences``` as a parameter to create a ```Distance Matrix``` and a ```Dictionary```. 


```julia
distanceMatrix = createDistanceMatrix(records_sequences) 

dict = createDictionary(records_sequences) 

```

```dict``` is a NamedTuple, which has ```dict_records``` and ```leaf_names```. ```leaf_names``` is a list of letters in the alphabet. Each letter presents a ```Record``` from ```records_sequences```. ```leaf_names``` and ```records_sequences``` are described as keys and values of the dictionary ```dict_records```. 

```julia
dict_records = dict.dict_records 
leaf_names = dict.leaf_names 
```

### 3. Guild Tree (UPGMA)

Here ```UPGMA``` is selected as a tree-based algorithm to perform multiple alignment. With the help of this tree, alignment is performed based on cluster of sequences. Sequences which are firstly clustered are firstly aligned. 

```tree``` is the created tree from UPGMA algorithm, using distance matrix and a list of names. 

```julia
tree = upgma(distanceMatrix,leaf_names)
```
Call ```guildTreeInstruction``` method to return list of clusters that need to be aligned: 

```julia
guildTree_instruction = guildTreeInstruction(tree)
```

### 4. Progressive Alignment

```guildTree_instruction``` shows pairs, which will be aligned, in order. A pair can be 2 sequences (pairwise alignment), or sequence-alignment or 2 alignments (multi sequences alignment).  

There are 2 ways to run the Progressive Alignment:

First way is using ```simpleGuildTree``` to simplify the ```guildTree_instruction```, then call ```progressiveAlignment```:
```julia
simpleinstruction = simpleGuildTree(guildTree_instruction)
result = progressiveAlignment(simpleinstruction, dict_records) 
```
Second one is using ```nestedInstruction``` to get the nested instruction from ```guildTree_instruction```. Then run ```progressiveAlignment2```:
```julia
nestedinstruction = nestedInstruction(guildTree_instruction)
result = progressiveAlignment2(nestedinstruction, dict_records)
```

Both ```result```s return a dictionary of Records, which contain all aligned sequences.

### 5. DataWriter

The aligned sequences are written in a text file format 

Call ```writeSequences``` method with output file name ("filename.txt") and dictionary containing representation string of the records (```result```). *50* in the below example is the sequence's number of characters printed on each line.

```
writeSequences("filename.txt", result,50) 
```

## Library

BioAlignment [GitHub](https://github.com/BioJulia/BioAlignments.jl.git)

FastaIO [Link](https://docs.juliahub.com/FastaIO/i12XQ/1.0.0/)

MCPhyloTree [GitHub](https://github.com/erathorn/MCPhyloTree.jl.git)

## Language

Julia Version 1.8.2 [Link](https://julialang.org/downloads/)


