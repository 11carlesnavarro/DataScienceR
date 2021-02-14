

class Sequence (object):


    protein_letters = 'ACDEFGHIKLMNPQRSTVWY'
    rna_letters = 'GAUC'
    dna_letters = 'GATC'
    molecular_weight_dict = {}

    def __init__(self, identifier, sequence):
        self.identifier = identifier
        self.sequence = sequence

    def get_identifier(self): 
        """
        Returns the identifier of the sequence
        """
        return self.identifier

    def get_sequence(self):
        """
        Returns the sequence of the protein or nucleic acid
        """
        return self.sequence

    def get_mw(self):
        """
        Returns the molecular weight of the sequence
        """

        unique_aminos = set(self.sequence)
        molecular_weight = 0
        
        for amino in unique_aminos:
            if amino in self.molecular_weight_dict:
                molecular_weight += self.sequence.count(amino)*self.molecular_weight_dict[amino]
        return molecular_weight
    
    def has_subsequence(self, subsequence):

        """
        Returns true if the parsed sequence is a subsequence of the protein
        """ 
        if subsequence in self.sequence:
            return True
        else:
            return False

    def __eq__(self, other):
        
        """
        Returns True if the sequences are equal, else return False
        """
        return self.get_sequence() == other.get_sequence()

    def __len__(self):
        
        """
        Returns the length of the sequence
        """
        return len(self.get_sequence())

    def __contains__(self, subsequence):
        
        """
        Returns a boolean deppending if the sequence is subsequence of the sequence. 
        """
        if isinstance(subsequence, Sequence):
            return subsequence.get_sequence() in self.get_sequence()
        else:
            raise ValueError("This is not a sequence instance")
    
    def __getitem__(self, key):
        """
        Returns the desired index of the sequence
        """
        return self.get_sequence()[key]
       
    def __add__(self, other):
        """
        If apply the + operator to two sequence objects it creates a new instance with the two joined sequences and an identifier with the format: id1+id2 
        """
        if type(self) == type(other):
            new_sequence = Sequence(identifier= f'{self.get_identifier()}+{other.get_identifier()}', sequence = self.get_sequence() + other.get_sequence())
            return new_sequence
        else:
            raise ValueError(f"{self.get_identifier()} and {other.get_identifier()} do not belong to the same Class")

    def __lt__(self, other):
        """
        Defines how to compare sequence instances. By molecular weight
        """
        return self.get_mw() > other.get_mw()

    def __hash__(self):
        """
        Function to allow the sequence instances to be keys of a hash or values of a set
        """
        return self.identifier.__hash__()

class NucleotideSequence (Sequence):

    # Define empty attributes that will be filled in the child classes

    translate_dict = {}
    stop_codons = []
    start_codons = []

    def __init__(self, identifier, sequence):

        super().__init__(identifier = identifier, sequence = sequence)

    def translate(self):
        
        ProteinSequence = "" 
        i = 0
        stop = True
        while i < len(self.sequence):
            (x,y,i) = (i, i+3,i+3) # Indexes to retrieve codons
            if self.sequence[x:y] in self.start_codons or stop == False: # Check if it is a start codon or if we can continue transcribing
                # Check for stop codons
                if self.sequence[x:y] in self.stop_codons:
                    yield ProteinSequence # I use a generator because a sequence can code more than one protein
                    ProteinSequence = ""
                    stop = True
                else:
                    ProteinSequence += self.translate_dict[self.sequence[x:y]]
                    stop = False

class DNASequence (NucleotideSequence):

    # Defining all the class atributes

    dna_table = {'CTT': 'L', 'ATG': 'M', 'ACA': 'T', 'ACG': 'T', 'ATC': 'I', 'ATA': 'I', 'AGG': 'R', 'CCT': 'P', 'AGC': 'S', 'AGA': 'R', 'ATT': 'I', 'CTG': 'L', 'CTA': 'L', 'ACT': 'T', 'CCG': 'P', 'AGT': 'S', 'CCA': 'P', 'CCC': 'P', 'TAT': 'Y', 'GGT': 'G', 'CGA': 'R', 'CGC': 'R', 'CGG': 'R', 'GGG': 'G', 'GGA': 'G', 'GGC': 'G', 'TAC': 'Y', 'CGT': 'R', 'GTA': 'V', 'GTC': 'V', 'GTG': 'V', 'GAG': 'E', 'GTT': 'V', 'GAC': 'D', 'GAA': 'E', 'AAG': 'K', 'AAA': 'K', 'AAC': 'N', 'CTC': 'L', 'CAT': 'H', 'AAT': 'N', 'CAC': 'H', 'CAA': 'Q', 'CAG': 'Q', 'TGT': 'C', 'TCT': 'S', 'GAT': 'D', 'TTT': 'F', 'TGC': 'C', 'TGG': 'W', 'TTC': 'F', 'TCG': 'S', 'TTA': 'L', 'TTG': 'L', 'TCC': 'S', 'ACC': 'T', 'TCA': 'S', 'GCA': 'A', 'GCC': 'A', 'GCG': 'A', 'GCT': 'A'}
    translate_dict = dna_table

    dna_weights = {'A': 347.0, 'C': 323.0, 'T': 322.0, 'G': 363.0}
    molecular_weight_dict = dna_weights

    dna_start_codons = ['TTG', 'CTG', 'ATG']
    start_codons = dna_start_codons

    dna_stop_codons = ['TAA', 'TAG', 'TGA']
    stop_codons = dna_stop_codons

    def __init__(self, identifier, sequence):

        super().__init__(identifier = identifier, sequence = sequence)

        for letter in self.sequence:
            if letter not in self.dna_letters:
                raise ValueError("Impossible to create instance: %s not possible" %(letter))
            
    
    def transcribe(self):

        RNASequence = ""

        for coding_nucleotide in self.sequence:
            if coding_nucleotide == 'T':
                RNASequence += 'U'
            else:
                RNASequence += coding_nucleotide
        
        
        return RNASequence

class RNASequence (NucleotideSequence):

    # Defining all the class atributes

    rna_table = {'GUC': 'V', 'ACC': 'T', 'ACA': 'T', 'ACG': 'T', 'GUU': 'V', 'AAC': 'N', 'AGG': 'R', 'UGG': 'W', 'AGC': 'S', 'AUC': 'I', 'AGA': 'R', 'AAU': 'N', 'ACU': 'T', 'CAC': 'H', 'GUG': 'V', 'CCG': 'P', 'CCA': 'P', 'AGU': 'S', 'CCC': 'P', 'GGU': 'G', 'UCU': 'S', 'GCG': 'A', 'CGA': 'R', 'CAG': 'Q', 'CGC': 'R', 'UAU': 'Y', 'CGG': 'R', 'UCG': 'S', 'CCU': 'P', 'GGG': 'G', 'GGA': 'G', 'GGC': 'G', 'GAG': 'E', 'UCC': 'S', 'UAC': 'Y', 'CGU': 'R', 'GAA': 'E', 'AUA': 'I', 'GCA': 'A', 'CUU': 'L', 'UCA': 'S', 'AUG': 'M', 'CUG': 'L', 'AUU': 'I', 'CAU': 'H', 'CUA': 'L', 'GCC': 'A', 'AAA': 'K', 'AAG': 'K', 'CAA': 'Q', 'UUU': 'F', 'GAC': 'D', 'GUA': 'V', 'UGC': 'C', 'GCU': 'A', 'UGU': 'C', 'CUC': 'L', 'UUG': 'L', 'UUA': 'L', 'GAU': 'D', 'UUC': 'F'}
    translate_dict = rna_table

    rna_weights = {'A': 363.0, 'C': 339.0, 'U': 340.0, 'G': 379.0}
    molecular_weight_dict = rna_weights

    rna_start_codons = ['UUG', 'CUG', 'AUG']
    start_codons = rna_start_codons

    rna_stop_codons = ['UAA', 'UAG', 'UGA']
    stop_codons = rna_stop_codons


    def __init__(self, identifier, sequence):

        super().__init__(identifier = identifier, sequence = sequence)

        for letter in self.sequence:
            if letter not in self.rna_letters:
                raise ValueError("Impossible to create instance: %s not possible" %(letter))
    
    def reverse_transcribe(self):

        DNASequence = ""

        for nucleotide in self.sequence:
            if nucleotide == 'U':
                DNASequence += 'T'
            else:
                DNASequence += nucleotide
        
        return DNASequence

class ProteinSequence (Sequence):
    
    # Defining all the class atributes

    protein_weights = {'A': 89.09, 'C': 121.16, 'E': 147.13, 'D': 133.1, 'G': 75.07, 'F': 165.19, 'I': 131.18, 'H': 155.16, 'K': 146.19, 'M': 149.21, 'L': 131.18, 'N': 132.12, 'Q': 146.15, 'P': 115.13, 'S': 105.09, 'R': 174.2, 'T': 119.12, 'W': 204.23, 'V': 117.15, 'Y': 181.19}
    molecular_weight_dict = protein_weights

    def __init__(self, identifier, sequence):

        super().__init__(identifier = identifier, sequence = sequence)

        for letter in self.sequence:
            if letter not in self.protein_letters:
                raise ValueError("Impossible to create instance: %s not possible" %(letter))

if __name__ == "__main__":

    sequence = "TTGAGATAA"
    identifier = ">B6I411"

    sequence2 = "CUGGUCACCACAUAAUAAUUGGUCUAAUAGUUG"

    my_rna= RNASequence(">BAXE", sequence2 )
    protein = DNASequence(identifier, sequence)
    







