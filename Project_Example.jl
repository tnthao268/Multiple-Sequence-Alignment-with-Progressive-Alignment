# include("source/MultiSequencesAlignment.jl")
include("source/DataReader.jl")
# include("source/DataWriter.jl")
# include("source/UPGMA.jl")

input_files_names = ["data/ChimpanzeeDND1.fasta","data/DogDND1.fasta","data/HomoSapienDND1.fasta","data/MouseDND1.fasta"]

records_sequences = readSequences_file(input_files_names)

