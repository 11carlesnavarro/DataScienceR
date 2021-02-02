
class Protein():
    
    def __init__(self, identifier, sequence):
        
        self.identifier = identifier
        self.sequence = sequence
    
    def get_identifier(self): 
        """
        Returns the identifier of the protein
        """
        return self.identifier

    def get_sequence(self):
        """
        Returns the sequence of the protein
        """
        return self.sequence

    def get_mw(self):
        """
        Returns the molecular weight of the protein
        """
        aminoacid_mw = {'A': 89.09, 'C': 121.16, 'E': 147.13, 'D': 133.1, 'G': 75.07, 'F': 165.19, 'I': 131.18, 'H': 155.16, 'K': 146.19, 'M': 149.21, 'L': 131.18, \
        'N': 132.12, 'Q': 146.15, 'P': 115.13, 'S': 105.09, 'R': 174.2, 'T': 119.12, 'W': 204.23, 'V': 117.15, 'Y': 181.19}
        
        unique_aminos = set(self.sequence)
        molecular_weight = 0
        
        for amino in unique_aminos:
            if amino in aminoacid_mw:
                molecular_weight += self.sequence.count(amino)*aminoacid_mw[amino]
        return molecular_weight
    
    def has_subsequence(self, subsequence):
        """
        Returns true if the parsed sequence is a subsequence of the protein
        """ 
        if subsequence in self.sequence:
            return True
        else:
            return False
    
    def get_length(self):
        """
        Returns the lenght of the sequence
        """
        return len(self.sequence)


def FASTA_iterator(fasta_filename):
    """"
    Reads a fasta file and return an object of the class Protein()
    """

    fasta_sequences = open(fasta_filename,"r")
    identifier = None
    sequence = []
    for line in fasta_sequences:
        line = line.rstrip()
        if line[0] == ">":
            if identifier: yield Protein(identifier, ''.join(sequence))
            identifier, sequence = line, []
        else:
            sequence.append(line)
    if identifier: yield Protein(identifier, ''.join(sequence))

    fasta_sequences.close


   