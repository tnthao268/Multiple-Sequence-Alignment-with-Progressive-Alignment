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
last_instruction = last(guildtree)

#progressive alignment
new_dict_records = align1(last_instruction,dict_records) #lasts so long