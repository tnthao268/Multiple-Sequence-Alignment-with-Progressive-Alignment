using Pkg
Pkg.add("BioAlignments")
using BioAlignments

#macro to check, if the sequence is a dna
macro dna_str(dna)
    return quote
        nu = "ATGC" #indexes of a string is Char-Type
        
        #when there is a index, which isn't a DNA nucleotide
        for i ∈ $(esc(dna)) # \in
            if i ∉ nu # \notin
                return 0
            end
        end
        return $(esc(dna))
    end
end

#Using macro dna
dna"ATAV"
dna"ATAC"

@dna_str "ATTA"
@dna_str "AXT"

#Function to check DNA structure
function check_dna(dna::String)
    nu = "ATGC"
    for i in dna
        if i ∉ nu
            return "-"
        end
    end
    return dna
end

