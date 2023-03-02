module DataReader

import Pkg
Pkg.add("FastaIO")
using FastaIO

export readSequences, FastaRecord

# create an immutable data type to store sequence and its description 
struct FastaRecord
    description :: AbstractString
    sequence :: AbstractString
end

function readSequences(input_file_name :: AbstractString)
    records = FastaRecord[] # records is type FastaRecord
    FastaReader(input_file_name) do fr
        for (desc, seq) in fr
            push!(records,FastaRecord("$desc","$seq"))
            
        end
    end
    records
end


end