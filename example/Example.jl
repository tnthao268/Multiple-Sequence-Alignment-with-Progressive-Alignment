include("../source/MultiSequencesAlignment.jl")

#read file and create list of Records
input_files_names = ["example/data/macaca mulatta mir126.fasta","example/data/pan troglodytes mir126.fasta","example/data/sus scrofa mir126.fasta","example/data/equus caballus mir126.fasta ", "example/data/homo sapien mir126.fasta"]

records_sequences = readSequences_file(input_files_names)

#createDistanceMatrix
distanceMatrix = createDistanceMatrix(records_sequences)

#create Dictionary and leaf_names
dict = createDictionary(records_sequences)
dict_records = dict.dict_records

leaf_names = dict.leaf_names

#create GuildTree
cluster = upgma(distanceMatrix, leaf_names)
print_ascii(cluster) # easy-to-view version of tree
guildtree = guildTreeInstruction(cluster)

#Progressive Alignment: There are two ways to run it and they return the same result
#First one: progressive alignment with a simplified instruction from the guild tree in list
simpleInstruction = simpleGuildTree(guildtree)
@time new_dict_records = progressiveAlignment(simpleInstruction,dict_records)

#Second one: progressive alignment with a nested instruction from the guild tree
nestedinstruction = nestedInstruction(guildtree)
@time new_dict_records2 = progressiveAlignment2(nestedinstruction, dict_records)

#Data Writer: Write the result of multi sequences alignment in a text file (.txt)

# first method of aligning
writeSequences("example/example.txt",new_dict_records,50)
writeSequences("example/example2.txt", new_dict_records2, 50)

