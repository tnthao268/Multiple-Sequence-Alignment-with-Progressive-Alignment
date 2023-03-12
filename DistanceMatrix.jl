#Create Distance Matrix
include("ReadPairwiseAlignment.jl")
    
#=  create a distance matrix from scores of pairwise alignments of sequences
    return distance matrix
=#
function createDistanceMatrix(seqs::Vector{Record})
    len = length(seqs)
    
    dm = zeros(len,len) #create a matrix with number of sequences
    for i in 1:len
        for j in i+1:len
            pair = [seqs[i],seqs[j]]
            aln = readPairwiseAlignment(pair)
            dm[i,j] = dm[j,i] = aln.score
        end
    end
    dm
end

#=  create a Dictionary: Each key is a letter in alphabet and its value is a Record.
    The list of keys in Dictionary is used for Leaf Names in Guild Tree
    parameter list of Records
    return a Dictionary of Records and a list of Leaf Names
=#
function createDictionary(records::Vector{Record})
    #define a leaf names (as key of Dictionary)
    names = Array{Char}(undef,length(records))
    names[1] = 'A'
    for x in 2:length(names)
        names[x] = names[x-1] + 1
    end

    #convert Char[] to String[]
    leaf_names = [string(x) for x in name]

    return (leaf_names = leaf_names,dict_records = Dict(zip(leaf_names,records)))
end