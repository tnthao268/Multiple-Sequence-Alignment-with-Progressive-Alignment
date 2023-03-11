# Dictionary of EDNAFull

using Pkg
Pkg.add("BioAlignments")
using BioAlignments
Pkg.add("BioSymbols")
using BioSymbols

using Test

#-----
#Data
EDNAFULL
gap_open = -5
gap_extend = -1
#---------------------------
#=  get pair score
    paremeter n1, n2 can be nucleotides or gaps
    return pair score
=#
function getPairScore(n1::Char,n2::Char,gap_open=-5,gap_extend=-1)
    if '-' in [n1,n2]
        if n1 == n2; return gap_extend; end
        return gap_open
    end
    return EDNAFULL[n1,n2]
end

#Test 
@test EDNAFULL['A','A'] == getPairScore('A','A')
@test gap_open == getPairScore('A','-')
@test getPairScore('A','-',-10) == -10
@test getPairScore('-','-') == gap_extend

#------------
#=  get sum of pairs score (SoP) at a position in the progressive alignment
    parameter seqs1 is a list of first sequence(s)
            seqs2 is a list of second sequence(s)
            pos1, pos2 are positions on seqs1 and seqs2 for SoP calculation
=#

function getSoP_at_a_pos(seqs1::Vector{String},seqs2::Vector{String},pos1::Int64, pos2::Int64)
    chars1 = [x[pos1] for x in seqs1]
    chars2 = [y[pos2] for y in seqs2]

    local score = 0
    for x in chars1
        for y in chars2
            score += getPairScore(x,y)
        end
    end
    return score
end

#Test
seqs1 = ["ATAT";
         "ATT";
         "-TTA"]
seqs2 = ["CAA";
         "--AA"]
#(['A','A','-'],['C','-']) --> (AC + A-)* 2 + -C + -- 
sop = getSoP_at_a_pos(seqs1,seqs2,1,1)
test_sop = (getPairScore('A','C') + getPairScore('A','-')) * 2 + getPairScore('-','C') + getPairScore('-','-')
@test sop == test_sop

#---------------------------------
#=  set values in Traceback Matrix
    parameter score from diagonal,vertical and horizon direction in Score Matrix 
=#
function setTraceback(dia::Int64,ver::Int64,hor::Int64)
    m = maximum([dia, ver, hor])
    if m == dia; return "dia"
    elseif m == ver; return "ver"
    else return "hor"
    end
end

#Test 
@test setTraceback(5,4,3) == "dia"
@test setTraceback(4,5,3) == "ver"

#-----------------------------------------------------------------------
using .DataReader
#Multi Sequences Alignment mit  Global Alignment (Needleman Wunsch Alignment)
#=  get alignment from reading the Traceback Matrix
    parameter seqs1 and seqs2 are the sequences starts with '-' 
    return 2 lists of aligned sequences
=# 
function getAlignment(seqs1::Vector{String},seqs2::Vector{String},traceback::Matrix{String})
    new_aln_seqs1 = ["" for x in 1:length(seqs1)]
    new_aln_seqs2 = ["" for y in 1:length(seqs2)]

    #Prepare length of unaligned sequences
    j = length(seqs1[1])
    i = length(seqs2[1])

    #read the alignment
    while i > 1 && j > 1
        if traceback[i,j] == "dia"
            new_aln_seqs1 = [seqs1[x][j] * new_aln_seqs1[x] for x in 1:length(seqs1)]
            new_aln_seqs2 = [seqs2[y][i] * new_aln_seqs2[y] for y in 1:length(seqs2)]
            i -= 1; j -= 1
        
        elseif traceback[i,j] == "ver"
            new_aln_seqs1 = ['-' * x for x in new_aln_seqs1]
            new_aln_seqs2 = [seqs2[y][i] * new_aln_seqs2[y] for y in 1:length(seqs2)]
            i -= 1

        elseif traceback[i,j] == "hor"
            new_aln_seqs1 = [seqs1[x][j] * new_aln_seqs1[x] for x in 1:length(seqs1)]
            new_aln_seqs2 = ['-' * y for y in new_aln_seqs2]
            j -= 1
        end
    end

    #For the first gaps
    while i > 1
        new_aln_seqs1 = ['-' * x for x in new_aln_seqs1]
        new_aln_seqs2 = [seqs2[y][i] * new_aln_seqs2[y] for y in 1:length(seqs2)]
        i -= 1
    end
    while j > 1
        new_aln_seqs1 = [seqs1[x][j] * new_aln_seqs1[x] for x in 1:length(seqs1)]
        new_aln_seqs2 = ['-' * y for y in new_aln_seqs2]
        j -= 1
    end
    
    new_aln_seqs1,new_aln_seqs2
end

#=  Global Alignment of (aligned) sequences by Needlemann Wunsch methode:
    create a score matrix and then read the best alignment from the right bottom 
    to the left top and create a traceback matrix
    parameter two lists of Records. These lists contain groups of sequences, which are read from Guild Tree for Progressive Alignment
    return Traceback Matrix
=#
function globalAlignment(aln_seqs1::Vector{Record},aln_seqs2::Vector{Record}, gap_open = -5)
    #check, if 2 lists of Records has minimal 1 Records
    @assert (length(aln_seqs1) >= 1 && length(aln_seqs2) >= 2) ||
            (length(aln_seqs2) >= 1 && length(aln_seqs1) >= 2)

    #List of sequences in Vector{String}
    seqs1 = [x.sequence for x in aln_seqs1] 
    seqs2 = [y.sequence for y in aln_seqs2]
    #println(seqs1,"\n",seqs2)

    #check, if the align sequences in the same list has same length
    @assert allequal([length(x) for x in seqs1])
    @assert allequal([length(y) for y in seqs2])
    
    #add a gap (-) at the first position of every sequence
    seqs1 = ['-' * x for x in seqs1]
    seqs2 = ['-' * y for y in seqs2]
    #println(seqs1,"\n",seqs2)

    #Prepare length of seuqences for the size of matrixes
    len1 = length(seqs1[1])
    len2 = length(seqs2[1])
    
    #Prepare two leer matrixes: Traceback Matrix and Score Matrix
    traceback = Matrix{String}(undef,len2,len1)
    scorematrix = zeros(Int64,len2,len1)

    #Initial
    scorematrix[1,1] = 0
    for i in 2:len2
        scorematrix[i,1] = scorematrix[i-1,1] + getSoP_at_a_pos(seqs1,seqs2,1,i)
    end
    for j in 2:len1
        scorematrix[1,j] = scorematrix[1,j-1] + getSoP_at_a_pos(seqs1,seqs2,j,1)
    end
    #println(scorematrix)
    
    #add values in Score Matrix
    for i in 2:len2
        for j in 2:len1
            ver = scorematrix[i-1,j] + gap_open*length(seqs1)*length(seqs2) # = get SoP('-'*length(seqs2),length(seqs1))
            hor = scorematrix[i,j-1] + gap_open*length(seqs1)*length(seqs2) # = get SoP('-'*length(seqs1),length(seqs2))
            dia = scorematrix[i-1,j-1] + getSoP_at_a_pos(seqs1,seqs2,j,i)

            scorematrix[i,j] = max(dia,ver,hor)
            traceback[i,j] = setTraceback(dia,ver,hor)
        end
    end
    
    #update the lists of Records with the alignment
    new_aln_seqs = getAlignment(seqs1,seqs2,traceback)
    for x in 1:length(aln_seqs1)
        aln_seqs1[x].sequence = new_aln_seqs[1][x]
    end
    for y in 1:length(aln_seqs2)
        aln_seqs2[y].sequence = new_aln_seqs[2][y]
    end
    println("scorematrix:",typeof(scorematrix))
    for x in 1:len2; println(scorematrix[x,:]); end
    println("traceback:")
    for x in 2:len2; println(traceback[x,2:len1]); end
    aln_seqs1,aln_seqs2
end

#Examples
aln_seqs1 = Record[Record("seq1","AA-GC"), Record("seq2","AATGC")]
aln_seqs2 = Record[Record("seq3","ACTC")]
globalAlignment(aln_seqs1,aln_seqs2)
aln_seqs1
aln_seqs2

aln_seqs3 = Record[Record("s1","AATCGAG"), Record("s2","AA-CGAG")]
aln_seqs4 = Record[Record("s3","ACCGAG"), Record("s4","ACGGAG")]
new_aln_seqs34 = [Record[Record("s1","AATCGAG"), Record("s2","AA-CGAG")],
                Record[Record("s3","A-CCGAG"), Record("s4","A-CGGAG")]]
globalAlignment(aln_seqs3,aln_seqs4) 
aln_seqs3
aln_seqs4
