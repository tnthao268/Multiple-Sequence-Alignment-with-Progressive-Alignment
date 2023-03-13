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

#Examples
createGroup([["A","B"],"C"])
createGroup("A")

#progressiveAlignment
function progressiveAlignment1(guild, dict_records::Dict{String,Record})
    #check, if each pair has 2 indexes
    @assert length(guild) == 2 "$guild is not a pair"

    #when the guild is pairwise alignment (for example guild = ["A","B"])
    if length(guild[1]) == length(guild[2]) == 1
        #pairwise alignment
        aln = readPairwiseAlignment([dict_records[x] for x in guild])
        
        #dictionary of pair of sequences, which will be aligned
        dict_aln = Dict(zip(guild,aln.aln_pair))
        
        #update the dictionary dict_records new aligned sequences
        merge!(dict_records,dict_aln)
        #println(dict_records)
    
    #when the guild is multi sequences alignment (for example guild = ["AB","C"]) 
    else
        #Preapare 2 groups of keys of sequences for alignment (for example from guild[1] = "AB", we have g1 = ["A","B"])
        g1 = [string(x) for x in guild[1]]
        g2 = [string(y) for y in guild[2]]

        #two lists of Records: parameter for msa_globalAlignment methode
        aln_seqs1 = [dict_records[x] for x in g1]
        aln_seqs2 = [dict_records[y] for y in g2]
        
        #Multi Sequences Alignment: align more than 2 sequences 
        aln = msa_globalAlignment(aln_seqs1,aln_seqs2)

        #2 dictionaries of sequences, which have been already aligned
        dict_aln1 = Dict(zip(g1,aln.aln_seqs1))
        dict_aln2 = Dict(zip(g2,aln.aln_seqs2))

        #update the dictionary dict_records new aligned sequences
        merge!(dict_records,dict_aln1)
        #println(dict_records)
        merge!(dict_records,dict_aln2)
        #println(dict_records)
    end
    dict_records
end

# Multi Sequences Alignment
function align1(guildTree,dict_records)
    if typeof(guildTree) == Vector{String}
        return dict_records #return dict_records
    end

    for x in guildTree
        println(x)
        if length(x) > 1
            progressiveAlignment1([createGroup(y) for y in x],align1(x,dict_records))
        end
    end
    dict_records
end

#Using
align1(guildTree[4], dict_records)
    
