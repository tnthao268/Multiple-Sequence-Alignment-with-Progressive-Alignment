using .DataReader

# check if 2 sequences are equal in length
function check(record1 :: FastaRecord, record2 :: FastaRecord)
    s = length(record1.sequence) == length(record2.sequence)
    return s
end

# test 
sequences1 = readSequences(dirname(@__FILE__()) * "/../data/input_test_sequences.faa")
check(sequences1[1],sequences1[2])
check(sequences1[1],sequences1[3])

