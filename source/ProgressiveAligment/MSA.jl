# Multi Sequences Alignment: Alignment of alignment and a sequence or of 2 alignments (not for Pairwise Alignment) 
using Pkg
Pkg.add("BioAlignments")
using BioAlignments
Pkg.add("BioSymbols")
using BioSymbols

#-----
#Data
EDNAFULL
gap_open = -5
gap_extend = -1
#---------------------------
#=  get pair score of second characters of 2 String parameters
    parameter each of n1, n2 is a string of 2 characters, which stand next to each other in the sequence.
              gap_open penalty and gap_extend penalty
    return pair score
=#
function getPairScore(n1::String,n2::String,gap_open=-5,gap_extend=-1)
    @assert length(n1) == length(n2) == 2
    c1 = n1[2]; c2 = n2[2]
    if '-' in [c1,c2]
        if c1 == c2
            return n1[1] == n2[1] == '-' ? gap_open : 0
        else
            return n1[1] != n2[1] && (n1[1] == c1 == '-' || n2[1] == c2 == '-') ?  gap_extend : gap_open
        end
    end
    return EDNAFULL[c1,c2]
end

#------------
#=  get sum of pairs score (SoP) at a position in the progressive alignment
    parameter seqs1 is a list of first sequence(s)
            seqs2 is a list of second sequence(s)
            pos1, pos2 are positions on seqs1 and seqs2 for SoP calculation
=#
function getSoP_at_a_pos(seqs1::Vector{String},seqs2::Vector{String},pos1::Int64, pos2::Int64)
    chars1 = pos1 > 1 ? [x[pos1-1:pos1] for x in seqs1] : ['A' * x[pos1] for x in seqs1]
    chars2 = pos2 > 1 ? [y[pos2-1:pos2] for y in seqs2] : ['A' * y[pos2] for y in seqs2]
    #println(chars1,"\n",chars2)
    local score = 0
    for x in chars1
        for y in chars2
            score += getPairScore(x,y)
        end
    end
    return score
end

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
    
    (new_aln_seqs1 = new_aln_seqs1, new_aln_seqs2 = new_aln_seqs2)
end

#=  Global Alignment of (aligned) sequences by Needlemann Wunsch methode:
    create a score matrix and then read the best alignment from the right bottom 
    to the left top and create a traceback matrix
    parameter two lists of Records. These lists contain groups of sequences, which are read from Guild Tree for Progressive Alignment
    return Traceback Matrix
=#
function msa_globalAlignment(aln_seqs1::Vector{Record},aln_seqs2::Vector{Record}, gap_open = -5)
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
            ver = scorematrix[i-1,j] + getSoP_at_a_pos(seqs1,[y[1:i-1] * '-' for y in seqs2],j,i) #add a gap in a copie of sequence, calculate SoP of opening gap.
            hor = scorematrix[i,j-1] + getSoP_at_a_pos([x[1:j-1] * '-' for x in seqs1],seqs2,j,i) #use copie to not change the sequence before the finish the scorematrix
            dia = scorematrix[i-1,j-1] + getSoP_at_a_pos(seqs1,seqs2,j,i)

            scorematrix[i,j] = max(dia,ver,hor)
            traceback[i,j] = setTraceback(dia,ver,hor)
        end
    end
    
    #update the lists of Records aln_seqs1, aln_seqs2 with the alignment
    new_aln_seqs = getAlignment(seqs1,seqs2,traceback)
    for x in 1:length(aln_seqs1)
        aln_seqs1[x].sequence = new_aln_seqs.new_aln_seqs1[x]
    end
    for y in 1:length(aln_seqs2)
        aln_seqs2[y].sequence = new_aln_seqs.new_aln_seqs2[y]
    end
    
    #print the score matrix, the traceback matrix and the alignment score
    println("scorematrix:")
    for x in 1:len2; println(scorematrix[x,:]); end
    println("traceback:")
    for x in 2:len2; println(traceback[x,2:len1]); end
    println("Alignment Score: $(scorematrix[len2,len1])")
    
    (aln_seqs1 = aln_seqs1, aln_seqs2 = aln_seqs2)
end

 #Example
aln_seqs1 = Record[Record("s1","AATCGAG"), Record("s2","AA-CGAG")]
aln_seqs2 = Record[Record("s3","ACCGAG"), Record("s4","ACGGAG")]
msa = msa_globalAlignment(aln_seqs1,aln_seqs2)
#=
scorematrix:
[0, -20, -40, -50, -70, -90, -110, -130]
[-20, 20, 0, -20, -40, -60, -70, -90]
[-40, 0, 4, -16, 0, -20, -40, -60]
[-60, -20, -16, -14, -14, 2, -18, -38]
[-80, -40, -36, -24, -30, 6, -14, 2]
[-100, -60, -20, -34, -40, -14, 26, 6]
[-120, -80, -40, -38, -50, -20, 6, 46]
traceback:
["dia", "dia", "hor", "hor", "hor", "dia", "hor"]
["ver", "dia", "hor", "dia", "hor", "hor", "hor"]
["ver", "dia", "dia", "dia", "dia", "hor", "dia"]
["ver", "dia", "ver", "dia", "dia", "dia", "dia"]
["dia", "dia", "ver", "dia", "ver", "dia", "hor"]
["ver", "ver", "dia", "dia", "dia", "ver", "dia"]
Alignment Score: 46
(Record[Record("s1", "AATCGAG"), Record("s2", "AA-CGAG")], Record[Record("s3", "AC-CGAG"), Record("s4", "AC-GGAG")])
=# 
msa.aln_seqs1 == aln_seqs1 # this means the variable is also updated after methode
#true

msa.aln_seqs2
#=
2-element Vector{Record}:
 Record("s3", "AC-CGAG")
 Record("s4", "AC-GGAG")
 =#
aln_seqs2
#=
2-element Vector{Record}:
 Record("s3", "AC-CGAG")
 Record("s4", "AC-GGAG")
 =#