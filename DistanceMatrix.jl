include("ReadAlignment.jl")

#=  create a distance matrix from scores of alignments of sequences
    return distance matrix
=#
function createDistanceMatrix(dict_seqs::Dict{String,String})
    names_list = [x for x in keys(dict_seqs)]
    len = length(dict_seqs)
    
    dm = zeros(len,len) #create a matrix with number of sequences
    for i in 1:len
        for j in i+1:len
            dict_pair = filter(s -> s.first == names_list[i] || s.first == names_list[j], dict_seqs)
            aln = readDNAAlignment(dict_pair)
            dm[i,j] = aln.score
        end
    end
    dm
end
