#ReadAlignment
struct ReadAlignment
    score::Int64
    aln_seq1::Vector{Char}
    aln_seq2::Vector{Char}
end

#---------------------------
# check if sequence is DNA
function check_DNA(seq::String)
    nu = "ATGC"
    for i in seq
        if i ∉ nu; return false; end
    end
    return true
end

#--------------------------------------------------------
#=  aligns two sequences and save the aligned sequences
    returns a ReadAlignment(score,aln_seq1,aln_seq2)
=#
function readDNAAlignment(dict_pair::Dict{String,String})
    
    #check, if there are 2 sequences for pairwise alignment and if they are DNA seuqences
    @assert length(dict_pair) == 2 "there are more or less than two sequences to do a pairwise alignment"
    foreach(keys(dict_pair)) do s
        if !check_DNA(dict_pair[s])
            println("s is not a DNA sequence")
        end
    end  

    #pairwise alignment
    scoremodel = AffineGapScoreModel(EDNAFULL,gap_open = -10, gap_extend = -1)
    seq = collect(values(dict_pair))
    res = pairalign(GlobalAlignment(),seq[1],seq[2],scoremodel)
    aln = alignment(res)

    aln_seq1 = [x for (x,_) in aln]
    aln_seq2 = [y for (_,y) in aln]
    
    println("score: $(score(res))\n$aln")
    ReadAlignment(score(res),aln_seq1,aln_seq2)
end

#--------------
#Using example
include("BioAlignment_Versuch.jl")
dict_pair = Dict("A" => seq, "B" => ref)
aln1 = readDNAAlignment(dict_pair)

#Compare to the results in BioAlignment_Versuch
@assert aln_seq == aln1.aln_seq2
@assert aln_ref == aln1.aln_seq1

#Print the instance of an object ReadAlignment
aln1.aln_seq1
aln1.score