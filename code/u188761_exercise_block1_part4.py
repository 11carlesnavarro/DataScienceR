

def FASTA_iterator(fasta_filename):
    "Reads a fasta file and return a tuple (identifier, sequence)"
    fasta_sequences = open(fasta_filename,"r")
    identifier = None
    sequence = []
    for line in fasta_sequences:
        line = line.rstrip()
        if line[0] == ">":
            if identifier: yield (identifier, ''.join(sequence)) # Return the the previous tuple each time it founds a new protein
            identifier, sequence = line, []
        else:
            sequence.append(line)
    if identifier: yield (identifier, ''.join(sequence)) # For the last protein

    fasta_sequences.close


def compare_fasta_file_identifiers (fasta_filenames_list):
    """
    Given a list of fasta files, returns a dictionary that contains:
    intersection: a set with the common identifiers found in all the files
    union: a set with all the identifiers (unique) found in all the files
    frequency: a dictionary with all the identifiers as keys and the number of files in which appears as values
    specific: a dictionary with the name of the input files as keys and a set with the specific identifiers as values
        
    """
    #### LOOKING FOR THE UNION ####
    
    identifiers = []
    union = set()
    for fasta in fasta_filenames_list:
        filenames = []
        for identifier, sequence in FASTA_iterator(fasta):
            # Create a list of lists. Each sublist is a file
            filenames.append(identifier.upper()) # Case-intensive
            union.add(identifier.upper())
        identifiers.append(filenames)
    
    #### LOOKING FOR THE INTERSECTION ####

    intersection = set()
    for name in identifiers[0]:
        if all(name in fasta for fasta in identifiers):
            intersection.add(name)
    
    #### LOOKING FOR THE FREQUENCIES ####

    frequency = {}
    for fasta in identifiers:
        fasta = set(fasta)
        for name in fasta:
            if name not in frequency: 
                frequency[name] = 1
            elif name in frequency:
                frequency[name] += 1
    
    #### LOOKING FOR SPECIFIC ####

    specific = {}
    for fasta in fasta_filenames_list:
        specific_id = set()
        for identifier, sequence in FASTA_iterator(fasta):
            # If the value in frequencies is 1 it is specific of that file
            if frequency[identifier.upper()] == 1:
                specific_id.add(identifier.upper())

        if specific_id: specific[fasta] = specific_id
        else: specific[fasta] = None
                
    return (intersection, union, frequency, specific)



