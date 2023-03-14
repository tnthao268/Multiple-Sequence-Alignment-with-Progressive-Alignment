include("source/MultiSequencesAlignment.jl")

#read file and create list of Records
input_files_names = ["data/ChimpanzeeDND1.fasta","data/DogDND1.fasta","data/HomoSapienDND1.fasta","data/MouseDND1.fasta"]

records_sequences = readSequences_file(input_files_names)

#createDistanceMatrix
distanceMatrix = createDistanceMatrix(records_sequences)

#create Dictionary and leaf_names
dict = createDictionary(records_sequences)
dict_records = dict.dict_records
leaf_names = dict.leaf_names

#create GuildTree
cluster = upgma(distanceMatrix, leaf_names)
guildtree = guildTreeInstruction(cluster)

#there are 2 ways to align multi sequences
# first: progressive alignment with a simplified instruction from the guild tree in list
simpleInstruction = simpleGuildTree(guildtree)
@time new_dict_records = progressiveAlignment(simpleInstruction,dict_records) #lasts so long

#second: progressive alignment with a nested instruction from the guild tree
nestedInstruction = getNestedInstruction(guildtree)
@time new_dict_records = progressiveAlignment2(nestedInstruction) #lasts maybe shorter, but still long

#Data Writer: Write the result of multi sequences alignment in a text file (.txt)
writeSequences("example.txt",new_dict_records,100)