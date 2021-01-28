
def calculate_aminoacid_frequencies(fasta_filename,subsequences_filename,number_of_repetitions, output_filename):
    "Reads a FASTA file and a subsequences file and returns returns an output file with the proportion of proteins in the FASTA containing at least N-times each of the subsequences"
    
    # READING THE FASTA FILE

    fasta_sequences = open(fasta_filename,"r")
    seqs = []
    prots = []
    n = -1
    for line in fasta_sequences:
        line = line.rstrip()
        if line[0] == ">":
            my_seq = ""
            seqs.append(None)
            prots.append(line)
            n += 1
        if line[0] != ">":
            my_seq = my_seq + line
            seqs[n] = my_seq
    fasta_sequences.close 

    # READING THE SUBSEQUENCES FILE

    subsequences = open(subsequences_filename, "r")
    subseqs = []
    for sub in subsequences:
        sub = sub.rstrip()
        subseqs.append(sub) # All the subsequences on an array
    subsequences.close

    # WRITING THE RESULTS ON THE OUPUT FILE

    Outfile = open(output_filename, "w")
    print("{0: <25} {1: >10}".format("#Number of proteins: ", len(seqs)), file = Outfile)
    print("{0: <25} {1: >10}".format("#Number of subsequences: ",len(subseqs)), file = Outfile)
    print("#subsequence proportions: ", file = Outfile)

    for subsequence in subseqs:
        n_proteins = 0
        for sequence in seqs:
            #Number of reps
            repetitions = 0
            for x in range(len(sequence) - len(subsequence) + 1): # Use thin algorithm instead of count() to take into account the overlapping sequences
                if sequence[x:x+len(subsequence)] == subsequence:
                    repetitions += 1
            if (repetitions >= number_of_repetitions):
                n_proteins += 1
        relative_freq = n_proteins/len(seqs)
        print("{0: <5} {1: >5}".format(subsequence, n_proteins),"\t", '%.4f' % relative_freq, file = Outfile)
    Outfile.close()



