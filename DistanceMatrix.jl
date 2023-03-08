#Create Distance Matrix
include("ReadAlignment.jl")
    
#=  create a distance matrix from scores of pairwise alignments of sequences
    return distance matrix
=#
function createDistanceMatrix(seqs::Vector{FastaRecord})
    len = length(seqs)
    
    dm = zeros(len,len) #create a matrix with number of sequences
    for i in 1:len
        for j in i+1:len
            pair = [seqs[i],seqs[j]]
            aln = readDNAAlignment(pair)
            dm[i,j] = dm[j,i] = aln.score
        end
    end
    dm
end