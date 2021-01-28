
def FASTA_iterator(fasta_filename):

    fasta_sequences = open(fasta_filename,"r")
    identifier = None
    sequence = []
    for line in fasta_sequences:
        line = line.rstrip()
        if line[0] == ">":
            if identifier: yield (identifier, ''.join(sequence))
            identifier, sequence = line, []
        else:
            sequence.append(line)
    if identifier: yield (identifier, ''.join(sequence))

    fasta_sequences.close

### EXERCISE 2 REPEATED ###

def get_proteins_ratio_by_residue_threshold(filename,residue,relative_threshold=0.03, absolute_threshold=10):
    
    """ 
    Reads a FASTA file and returns the ratio of proteins with a relative threshold for a 
    given residue higher than the one in the input and the same for the absolute

    """
    
    n_higher = 0 # Number of proteins with abs and relative frequencies higher than the given ones
    n_proteins = 0
    for identifier, sequence in FASTA_iterator(filename):
        abs_freq = sequence.count(residue)
        length = len(sequence)
        n_proteins += 1
        rel_freq = abs_freq/length
        if (rel_freq >= relative_threshold and abs_freq >= absolute_threshold):
            n_higher += 1 # Number of proteins that it's over the thresholds
    return n_higher/n_proteins


def print_sequence_summary(filename, output_filename, first_n=10, last_m=10):
    
    """
    Reads a FASTA file and returns a file with the protein identifier, the first N-aminoacids, 
    the last M-aminoacids and the absolute frequency in the protein of all te aminoacids found in the protein
    """ 

    # WRITING THE OUTPUT FILE

    Outfile= open(output_filename,'w')
    for identifier, sequence in FASTA_iterator(filename):
        freqs = []
        first = sequence[0:first_n] # First N aminoacids
        last = sequence[-last_m:] # Last M aminoacids
        for amino in sequence:
            freq = amino + ':' + str(sequence.count(amino))
            freqs.append(freq) 
        freqs = set(freqs) # Total n of amino in the sequence. set() to eliminate repetitions
        print(identifier[1:], '\t', first, '\t', last," \t ", end = " ", file=Outfile) 
        for f in freqs:
            print( f, end = " ", file = Outfile) # Print the frequencies
        print("", file = Outfile)
    Outfile.close()



### MAXIMUM LENGTH ###

def get_max_sequence_length_from_FASTA_file (fasta_filename):

    
    """
    Given a multiline fasta file, returns the length of the sequence with the maximum length
    
    """
    max_length = 0
    for identifier, sequence in FASTA_iterator(fasta_filename):
        length = len(sequence)
        if length > max_length:
            max_length = length
    return max_length


### MINIMUM LENGTH ###

def get_min_sequence_length_from_FASTA_file (fasta_filename):
    
    """
    Given a multiline fasta file, returns the length of the sequence with the maximum length
    """
    min_length = 0
    for identifier, sequence in FASTA_iterator(fasta_filename):
        length = len(sequence)
        if length < min_length or min_length ==0 :
            min_length = length
    return min_length



### LONGEST SEQUENCE ###

def get_longest_sequences_from_FASTA_file (fasta_filename):

    """
    Given a multiline fasta file, returns a list of tuples corresponding to the sequences with maximum length
    """
    max_length = 0
    my_longest_sequences = []
    for identifier, sequence in FASTA_iterator(fasta_filename):
        length = len(sequence)
        my_sequence = (identifier,sequence)

        if length > max_length:
            max_length = length
            my_longest_sequences = [] # If a sequence is longer than the max, a new empty list is created
            my_longest_sequences.append(my_sequence)
        elif length == max_length:
            my_longest_sequences.append(my_sequence)
    
    my_longest_sequences = sorted(my_longest_sequences, key=lambda tup: tup[0].upper())
    
    return my_longest_sequences



### SHORTEST SEQUENCE ###

def get_shortest_sequences_from_FASTA_file (fasta_filename):

    """
    Given a multiline fasta file, returns a list of tuples corresponding to the sequences with minimum length
    """

    min_length = 0
    my_shortest_sequences = []
    for identifier, sequence in FASTA_iterator(fasta_filename):
        length = len(sequence)
        my_sequence = (identifier, sequence)

        if length < min_length or min_length == 0:
            min_length = length
            my_shortest_sequences = []
            my_shortest_sequences.append(my_sequence)
        elif length == min_length:
            my_shortest_sequences.append(my_sequence)
    
    my_shortest_sequences = sorted(my_shortest_sequences, key=lambda tup: tup[0].upper())

    return my_shortest_sequences



### MOLECULAR WEIGHTS ###
    
def get_molecular_weights(fasta_filename):

    """
    Given a protein fasta file, returns a dictionary with the molecular weights of all the proteins in the file.
    """

    aminoacid_mw = {'A': 89.09, 'C': 121.16, 'E': 147.13, 'D': 133.1, 'G': 75.07, 'F': 165.19, 'I': 131.18, 'H': 155.16, 'K': 146.19, 'M': 149.21, 'L': 131.18, \
         'N': 132.12, 'Q': 146.15, 'P': 115.13, 'S': 105.09, 'R': 174.2, 'T': 119.12, 'W': 204.23, 'V': 117.15, 'Y': 181.19}
    
    protein_and_weight = {}
    for identifier, sequence in FASTA_iterator(fasta_filename):
        unique_aminos = set(sequence) # With this set we get the different aminos in the protein. Without repeats.
        molecular_weight = 0

        for amino in unique_aminos:
            if amino in aminoacid_mw:
                molecular_weight += sequence.count(amino)*aminoacid_mw[amino] # Each amino multiplied by its mol weight
        protein_and_weight[identifier] = molecular_weight

    return protein_and_weight



### LIGHTER MOLECULE ###

def get_sequence_with_min_molecular_weight (fasta_filename):
    
    """
    Given a protein fasta file, returns a tuple with the (identifier, sequence) of the protein with the lowest molecular weight.
    """

    aminoacid_mw = {'A': 89.09, 'C': 121.16, 'E': 147.13, 'D': 133.1, 'G': 75.07, 'F': 165.19, 'I': 131.18, 'H': 155.16, 'K': 146.19, 'M': 149.21, 'L': 131.18, \
        'N': 132.12, 'Q': 146.15, 'P': 115.13, 'S': 105.09, 'R': 174.2, 'T': 119.12, 'W': 204.23, 'V': 117.15, 'Y': 181.19}
    
    min_molecular_weight = 0

    for identifier, sequence in FASTA_iterator(fasta_filename):
        unique_aminos = set(sequence)
        molecular_weight = 0

        for amino in unique_aminos:
            if amino in aminoacid_mw:
                molecular_weight += sequence.count(amino)*aminoacid_mw[amino] 

        if molecular_weight < min_molecular_weight or min_molecular_weight == 0:
            min_molecular_weight = molecular_weight
            lighter_protein = (identifier,sequence)

    return lighter_protein


### MEAN MOLECULAR WEIGHT ###

def get_mean_molecular_weight(fasta_filename):

    """
    Given a protein fasta file, returns the mean of the molecular weights of all the proteins
    """
    aminoacid_mw = {'A': 89.09, 'C': 121.16, 'E': 147.13, 'D': 133.1, 'G': 75.07, 'F': 165.19, 'I': 131.18, 'H': 155.16, 'K': 146.19, 'M': 149.21, 'L': 131.18, \
        'N': 132.12, 'Q': 146.15, 'P': 115.13, 'S': 105.09, 'R': 174.2, 'T': 119.12, 'W': 204.23, 'V': 117.15, 'Y': 181.19}
    
    total_molecular_weight = 0    
    n = 0

    for identifier, sequence in FASTA_iterator(fasta_filename):
        single_aminos = set(sequence)
        for amino in single_aminos:
            if amino in aminoacid_mw:
                total_molecular_weight += sequence.count(amino)*aminoacid_mw[amino] 
        n += 1
    mean_molecular_weight = total_molecular_weight/n
    return mean_molecular_weight


