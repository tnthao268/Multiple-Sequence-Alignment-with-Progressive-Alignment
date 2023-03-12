# include("source/MultiSequencesAlignment.jl")
include("source/DataReader.jl")
# include("source/DataWriter.jl")
# include("source/UPGMA.jl")

sequence1 = readSequences("data/ChimpanzeeDND1.fasta")
sequence2 = readSequences("data/DogDND1.fasta")
sequence3 = readSequences("data/HomoSapienDND1.fasta")
sequence4 = readSequences("data/MouseDND1.fasta")