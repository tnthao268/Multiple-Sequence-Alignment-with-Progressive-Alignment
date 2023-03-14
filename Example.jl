include("source/MultiSequencesAlignment.jl")

#read file and create list of Records
input_files_names = ["data/macaca mulatta mir126.fasta","data/pan troglodytes mir126.fasta","data/sus scrofa mir126.fasta","data/equus caballus mir126.fasta "]

records_sequences = readSequences_file(input_files_names)

#createDistanceMatrix
distanceMatrix = createDistanceMatrix(records_sequences)

#create Dictionary and leaf_names
dict = createDictionary(records_sequences)
dict_records = dict.dict_records
leaf_names = dict.leaf_names

#create GuildTree
cluster = upgma(distanceMatrix, leaf_names)
nj = neighbor_joining(distanceMatrix,leaf_names) # try with tree built from neighbor joining algorithm
print_ascii(cluster) # easy-to-view version of tree

print_ascii(nj)
guildtree = guildTreeInstruction(cluster)
guildtree_nj = guildTreeInstruction(nj)

#there are 2 ways to align multi sequences
# first: progressive alignment with a simplified instruction from the guild tree in list
simpleInstruction = simpleGuildTree(guildtree)
@time new_dict_records = progressiveAlignment(simpleInstruction,dict_records) #lasts so long

simpleInstruction_nj = simpleGuildTree(guildtree_nj)
@time new_dict_records_nj = progressiveAlignment(simpleInstruction_nj,dict_records)

#second: progressive alignment with a nested instruction from the guild tree
nestedinstruction = nestedInstruction(guildtree)
@time new_dict_records = progressiveAlignment2(nestedinstruction, dict_records) #lasts maybe shorter, but still long

#Data Writer: Write the result of multi sequences alignment in a text file (.txt)
writeSequences("example.txt",new_dict_records,50)

writeSequences("example_nj.txt",new_dict_records,50)