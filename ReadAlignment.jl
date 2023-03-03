export check_DNA; readDNAAlignment

#--------------
#ReadAlignment
struct ReadAlignment
    score::Int64
    dict_alnpair::Dict{String,String}
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
    name = reverse(collect(keys(dict_pair))) #list of names of sequences (use reverse because LIFO access)
    seq = reverse(collect(values(dict_pair))) #list of sequences
    res = pairalign(GlobalAlignment(),seq[1],seq[2],scoremodel)
    aln = alignment(res)

    #sequences in String and a dictionary of them (name, aligned sequence)
    aln_seq1 = reduce(*,[x for (x,_) in aln])
    aln_seq2 = reduce(*,[y for (_,y) in aln])
    dict_alnpair = Dict(zip(name,[aln_seq1,aln_seq2]))
    
    println("score: $(score(res))\n$aln")
    ReadAlignment(score(res),dict_alnpair)
end