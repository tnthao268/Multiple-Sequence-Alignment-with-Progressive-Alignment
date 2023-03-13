
include("DataReader.jl")

# using .DataReader

# import .DataReader.Record

export writeSequences



# function to find the key with maximum length in a dictionary
function key_max(dict::Dict{Any,Any})
    key_max = ""
    for key in keys(dict)
        if length(split(key," ")[1]) > length(key_max) # gets only "NC..." part of description
            key_max = split(key," ")[1]
        end
    end
    return key_max
end

# function to write with n line distances- with n sequences n+1 times \n 
# function line_distance()

# function to split sequences into chunks with char_each_line characters each
# result of this function is a list of substrings from a sequence
function split_sequence(str::String,char_each_line::Integer)
    substr = String[]
    for x in 1:char_each_line:length(str)
        if x+char_each_line-1 > length(str)
            push!(substr,str[x:length(str)])
        else
            push!(substr,str[x:x+char_each_line-1])
        end
    end
    return substr
end


# function to write aligned sequences into a text file, 
# with parameter as: output file's name and dictionary with key as alphabets, values as record
function writeSequences(output_file_name::String,dict_seq::Dict{String,Record},num_char_each_line::Integer)
    output_file = open(output_file_name, "w")
    sequences_dict = Dict() # new dict to store all aligned sequences, with key: description, value:sequence

    for val in values(dict_seq)
        sequence = val.sequence
        description = val.description
         # add sequence of each record to the sequences_dict dictionary
        merge!(sequences_dict,Dict(description => sequence))

    end
    
    x = 1 # x is the number of lines that each sequence will be printed out

    
    for sequence in values(sequences_dict)
        j = length(split_sequence(sequence,num_char_each_line)) 
        if j > x
            x = j # x: greatest number of lines that sequences need to be printed out
        end
    end

    
    for i in 1:x 

        for (description,sequence) in sequences_dict
       
                # write description right-justified in the beginning (write only "NM..." part of description)
                # and sequence with maximal length (num_char_each_line) next to it in each line
                write(output_file,rpad(split(description," ")[1],length(key_max(sequences_dict))) * " " * split_sequence(sequence,num_char_each_line)[i] ) 
                write(output_file,"\n")
            
        end
    write(output_file,"\n")


    end
    
    close(output_file)
    # return sequences_dict
    # println(length(values(sequences_dict)))


end







