import Pkg
Pkg.add("FastaIO")
using FastaIO

export readSequences, Record

# create an immutable data type to store sequence and its description 
mutable struct Record
    description :: AbstractString
    sequence :: AbstractString
end

# input_files_names = []

# function to return list of records based on name of the fasta filen(input_file_name)
function readSequences(input_file_name :: AbstractString)
    records = Record[] # records is type FastaRecord
    FastaReader(input_file_name) do fr
        for (desc, seq) in fr
            push!(records,Record("$desc","$seq"))
            
        end
    end
    records
end

# function to return a list of records based on a list of fasta files
function readSequences_file(input_files_names :: Vector{String}) # input_files_names is list of fasta files' names
    records = Record[] # records is type FastaRecord
    for name in input_files_names
        FastaReader(name) do fr
            for (desc, seq) in fr
                push!(records,Record("$desc","$seq"))
            
             end
        end
    end
    records
end



# print(sequences)