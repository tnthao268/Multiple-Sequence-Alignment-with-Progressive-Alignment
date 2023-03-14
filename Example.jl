include("source/MultiSequencesAlignment.jl")

#read file and create list of Records
input_files_names = ["data/macaca mulatta mir126.fasta","data/pan troglodytes mir126.fasta","data/sus scrofa mir126.fasta","data/equus caballus mir126.fasta ", "data/homo sapien mir126.fasta"]

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


#there are 2 ways to align multi sequences
# first: progressive alignment with a simplified instruction from the guild tree in list
simpleInstruction = simpleGuildTree(guildtree)
@time new_dict_records = progressiveAlignment(simpleInstruction,dict_records) #lasts so long



#second: progressive alignment with a nested instruction from the guild tree
nestedinstruction = nestedInstruction(guildtree)
<<<<<<< HEAD
@time new_dict_records = progressiveAlignment2(nestedinstruction,dict_records) #lasts maybe shorter, but still long
=======
@time new_dict_records = progressiveAlignment2(nestedinstruction, dict_records) #lasts maybe shorter, but still long
>>>>>>> 8db4c90934dd11aa3164f121e3a35fb158288127

#Data Writer: Write the result of multi sequences alignment in a text file (.txt)

# first method of aligning
writeSequences("example.txt",new_dict_records,50)

# second method of aligning
writeSequences("example2.txt",new_dict_records2,50)
