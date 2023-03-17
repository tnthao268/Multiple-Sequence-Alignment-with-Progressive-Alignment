include("../ReadPairwiseAlignment.jl")
include("MSA.jl")

#--------------------------------------------
#=  convert each group in pair in String to simplify the guildTreeInstruction
    for example: complex pair [["A","B"],"C"] => simple pair ["AB","C"]
=#
function createGroup(vector)
    group = ""
    function readGroup(each)
        if typeof(each) == String
            group *= each
        else
            for x in each
                readGroup(x)
            end
        end
    end
    readGroup(vector)
    group
end

#-------------------------------------
# simplify the guild tree instruction
function simpleGuildTree(guildTree)
    simple = []
    for pair in guildTree
        simple_pair = [createGroup(each) for each in pair]
        push!(simple,simple_pair)
    end
    simple
end

#---------------------------
#=  Progressive Alignment
    parameter guildtree is a the order of alignment
              dict_records is a dictionary with records as values and their pseudo names as keys
    return result of progressive alignment
=#

function progressiveAlignment(simpleguildtree, dict_records::Dict{String,Record}, gap_open = -16, gap_extend = -1)
    dict = copy(dict_records)
    
    #check, if each pair has 2 indexes
    for guild in simpleguildtree
        @assert length(guild) == 2 "$guild is not a pair"
    end

    for guild in simpleguildtree

        #when the guild is pairwise alignment (for example guild = ["A","B"])
        if length(guild[1]) == length(guild[2]) == 1
            #pairwise alignment
            aln = readPairwiseAlignment([dict[x] for x in guild], gap_open, gap_extend)
            
            #dictionary of pair of sequences, which will be aligned
            dict_aln = Dict(zip(guild,aln.aln_pair))
            
            #update the dictionary dict_records new aligned sequences
            merge!(dict,dict_aln)
            #println(dict_records)
        
        #when the guild is multi sequences alignment (for example guild = ["AB","C"]) 
        else
            #Preapare 2 groups of keys of sequences for alignment (for example from guild[1] = "AB", we have g1 = ["A","B"])
            g1 = [string(x) for x in guild[1]]
            g2 = [string(y) for y in guild[2]]

            #two lists of Records: parameter for msa_globalAlignment methode
            aln_seqs1 = [dict[x] for x in g1]
            aln_seqs2 = [dict[y] for y in g2]
            
            #Multi Sequences Alignment: align more than 2 sequences 
            aln = msa_globalAlignment(aln_seqs1,aln_seqs2, gap_open, gap_extend)

            #2 dictionaries of sequences, which have been already aligned
            dict_aln1 = Dict(zip(g1,aln.aln_seqs1))
            dict_aln2 = Dict(zip(g2,aln.aln_seqs2))

            #update the dictionary dict_records new aligned sequences
            merge!(dict,dict_aln1)
            #println(dict_records)
            merge!(dict,dict_aln2)
            #println(dict_records)
        end
    end
    dict
end
