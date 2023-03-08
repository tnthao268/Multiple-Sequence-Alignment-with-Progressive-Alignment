module DataReader

import Pkg
Pkg.add("FastaIO")
using FastaIO

export readSequences, Record

# create an immutable data type to store sequence and its description 
struct Record
    description :: AbstractString
    sequence :: AbstractString
end

function readSequences(input_file_name :: AbstractString)
    records = Record[] # records is type FastaRecord
    FastaReader(input_file_name) do fr
        for (desc, seq) in fr
            push!(records,Record("$desc","$seq"))
            
        end
    end
    records
end


end