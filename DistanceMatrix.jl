include("ReadAlignment.jl")

#=  create a distance matrix from scores of alignments of sequences
    return distance matrix
=#
function distanceMatrix(seqs::Dict{String,String})
    seqs_list= [x for x in values(seqs)] #list of only sequences
    len = length(seqs)
    
    dm = zeros(len,len) #create a matrix with number of sequences
    for i in 1:len
        for j in i+1:len
            aln = readDNAAlignment(seqs_list[i],seqs_list[j])
            dm[i,j] = aln.score
        end
    end
    dm
end
