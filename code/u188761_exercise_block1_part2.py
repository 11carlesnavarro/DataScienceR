
def get_proteins_ratio_by_residue_threshold(filename,residue,relative_threshold=0.03, absolute_threshold=10):
    "Reads a FASTA file and returns the ratio of proteins with a relative threshold for a given residue higher than the one in the input and the same for the absolute"
    # Reading the file, put the names and the sequences in two different arrays.
    fasta = open(filename,"r")
    seqs = []
    prots = []
    n = -1
    for line in fasta:
        line = line.rstrip()
        if line[0] == ">":
            length = 0
            my_seq = ""
            seqs.append(None)
            prots.append(line)
            n += 1
        if line[0] != ">":
            # Join the lines of the sequence and put the full sequence in the position n of the arrray
            my_seq = my_seq + line 
            seqs[n] = my_seq
    fasta.close 
    
    n_higher = 0 # Number of proteins with abs and relative frequencies higher than the given ones

    for element in seqs:
        abs_freq = element.count(residue)
        length = len(element)
        rel_freq = abs_freq/length
        if (rel_freq >= relative_threshold and abs_freq >= absolute_threshold):
            n_higher += 1 # Number of proteins that it's over the thresholds
    return n_higher/len(seqs)


def print_sequence_summary(filename, output_filename, first_n=10, last_m=10):
    "Reads a FASTA file and returns a file with the protein identifier, the first N-aminoacids, the last M-aminoacids and the absolute frequency in the protein of all te aminoacids found in the protein" 
    
    # READING THE FILE, put the names and the sequences in two different arrays.

    fasta = open(filename,"r")
    seqs = []
    prots = []
    n = -1
    for line in fasta:
        line = line.rstrip()
        if line[0] == ">":
            length = 0
            my_seq = ""
            seqs.append(None)
            prots.append(line)
            n += 1
        if line[0] != ">":
            my_seq = my_seq + line
            seqs[n] = my_seq
    fasta.close 

    # WRITING THE OUTPUT FILE

    Outfile= open(output_filename,'w')
    m = 0
    for element in seqs:
        freqs = []
        first = element[0:first_n] # First N aminoacids
        last = element[-last_m:] # Last M aminoacids
        for amino in element:
            freq = amino + ':' + str(element.count(amino))
            freqs.append(freq) 
        freqs = set(freqs) # Total n of amino in the sequence. set() to eliminate repetitions
        print(prots[m][1:], '\t', first, '\t', last," \t ", end = " ", file=Outfile) 
        for f in freqs:
            print( f, end = " ", file = Outfile) # Print the frequencies
        print("", file = Outfile)
        m += 1
    Outfile.close()


if __name__ == "__main__":

    get_proteins_ratio_by_residue_threshold("example_fasta_file.fa.txt","Q")