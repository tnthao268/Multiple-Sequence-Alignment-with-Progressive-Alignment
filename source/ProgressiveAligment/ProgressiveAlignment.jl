include("../ReadPairwiseAlignment.jl")
include("MSA.jl")

#=  Progressive Alignment
    parameter guildtree is a the order of alignment
              dict_records is a dictionary with records as values and their pseudo names as keys
    return result of progressive alignment
=#
function progressiveAlignment(guildtree::Vector{Vector{String}}, dict_records::Dict{String,Record})
    #check, if each pair has 2 indexes
    for guild in guildtree
        @assert length(guild) == 2 "$guild is not a pair"
    end

    for guild in guildtree

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
    end
    dict_records
end

#Example
s1 = Record("Seq1","TCAGGATGAAC")
s2 = Record("Seq2","ATCACGATGAACC")
s3 = Record("Seq3","ATCAGGAATGAATCC")
s4 = Record("Seq4","TCACGATTGAATCGC")
s5 = Record("Seq5","TCAGGAATGAATCGC")
records = [s1,s2,s3,s4,s5] #return from DataReader

include("../DistanceMatrix.jl")
dict_records = createDictionary(records).dict_records
dm = createDistanceMatrix(records)

guildTree = [["B","E"],["A","D"],["C","AD"],["CAD","BE"]]
guildTree2 = [["B","E"],["A","BE"],["C","D"],["CD","ABE"]]

progressiveAlignment(guildTree2,dict_records)

#=  Result with guildTree1
Alignment Score: 248
Dict{String, Record} with 5 entries:
  "B" => Record("Seq2", "ATCACGA-TGAA-C-C")
  "A" => Record("Seq1", "-TCAGGAT-GAA---C")
  "C" => Record("Seq3", "ATCAGGAATGAATC-C")
  "D" => Record("Seq4", "-TCACGATTGAATCGC")
  "E" => Record("Seq5", "-TCAGGAATGAATCGC")
=#

#= Result with guildTree2
Alignment Score: 266
Dict{String, Record} with 5 entries:
  "B" => Record("Seq2", "ATCACGA-TGAA-C-C")
  "A" => Record("Seq1", "-TCAGGA-TGAA---C")
  "C" => Record("Seq3", "ATCAGGAATGAATC-C")
  "D" => Record("Seq4", "-TCACGATTGAATCGC")
  "E" => Record("Seq5", "-TCAGGAATGAATCGC")
=#