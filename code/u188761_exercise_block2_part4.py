class IncorrectSequenceLetter(ValueError):
    """ Exception when a letter that not belongs to the sequence alphabet is found"""
    def __init__(self, letter, class_name):
        self.letter = letter
        self.class_name = class_name.__class__.__name__
    
    def __str__(self):
        return "The sequence item %s is not found in the alphabet of class %s" %(self.letter, self.class_name)
