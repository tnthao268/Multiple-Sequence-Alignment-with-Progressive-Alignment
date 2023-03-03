export check_DNA; readDNAAlignment
include("source/DataReader.jl")

#-------------------
#AlignmentRecord
struct AlignmentRecord
    description::String
    aln_seq::String
end

#--------------
#ReadAlignment
struct ReadAlignment
    score::Int64
    aln_pair::Vector{AlignmentRecord}
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
function readDNAAlignment(pair::Vector{FastaRecord})
    
    #check, if there are 2 sequences for pairwise alignment and if they are DNA seuqences
    @assert length(pair) == 2 "there are more or less than two sequences to do a pairwise alignment"
    foreach(pair) do s
         @assert check_DNA(s.sequence) "$(s.description) is not a DNA sequence"
    end

    #pairwise alignment
    scoremodel = AffineGapScoreModel(EDNAFULL,gap_open = -10, gap_extend = -1)
    res = pairalign(GlobalAlignment(),pair[1].sequence,pair[2].sequence,scoremodel)
    aln = alignment(res)

    #sequences in String and create list of AlignmentRecord(description, aligned sequence)
    aln_seq1 = reduce(*,[x for (x,_) in aln])
    aln_seq2 = reduce(*,[y for (_,y) in aln])
    aln_records = [AlignmentRecord(pair[1].description, aln_seq1),
                   AlignmentRecord(pair[2].description, aln_seq2)]
    
    println("score: $(score(res))\n$aln")
    ReadAlignment(score(res),aln_records)
end