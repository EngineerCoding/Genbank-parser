from re import match, split, IGNORECASE

from .location_parser import RemoteLocation


def parse_origin(gbp):
    """ The main method which parses the ORIGIN to a Sequence object.

    Parameters:
        gbp - GenbankParser object
            The parser which holds the file pointer of the genbank
            file.
    Returns:
        A Sequence object.
    """
    # Check if the header is there
    gbp.handle_keyword('ORIGIN', do_split=False, remove_keyword=False)
    sequence = ''
    # Keep reading lines when they start with a number
    line = gbp.read_valid_line().strip()
    while match('^\d+.*', line):
        # Remove the whitespace
        splitted = split('\s+', line)
        del splitted[0]  # Removes the number which is irrelevant to us
        sequence += ''.join(splitted).upper()
        line = gbp.read_valid_line().strip()
    return Sequence(sequence)


class Sequence(object):
    """ A Sequence object can be any sequence a string can represent,
    however the most likely sequences will be a DNA, RNA or protein
    sequence. Along with this sequence an accession can be set, just
    for the sake of relating a RemoteLocation to a sequence.

    A Location object can be used to retrieve the sequence which
    corresponds with the actual locations of the Location object.
    However, this object only does this based on the 'get_range'
    method of a Location object. It is recommended to use the
    method on a Location object to turn a Location into a sequence.

    With DNA and RNA sequences also the complement of the sequence
    can be retrieved.
    """

    def __init__(self, sequence):
        """ Creates a Sequence with the given string sequence and
        an accession which can be used by 'get_accession'/
        'set_accesion'.

        Parameters:
            sequence - string
                A string representing this Sequence object
        """
        self.sequence = sequence
        self.accession = None

    def get_sequence(self):
        """ Retrieves the string sequence of this object """
        return self.sequence

    def set_accession(self, accession):
        """ Sets the accession for this object.

        Parameters:
            accession - string
                A string representing the accession
        """
        self.accession = accession

    def get_accession(self):
        """ Retrieves the accession string for this object

        Returns:
            When the accession has been set, this will return the
            accession string. Otherwise, this will return None.
        """
        return self.accession

    def length(self):
        """ Retrieves the length of this sequence

        Returns:
            An int representing the length of this sequence
        """
        return len(self.sequence)

    def get_sequence_from_location(self, location, sequence=None):
        """ Gets the sequence for the location. For this the method
        'get_range' is used for to determine what parts to use.
        TODO: Currently RemoteLocation is a special case, but this
        should be moved to the actual RemoteLocation object.

        Parameters:
            location - Location object
                A location object which specifies the range of the
                sequence
            sequence - Sequence object. Default: None.
                Alternative sequence to use for a RemoteLocation, will
                be removed in the future.
        Returns:
            A string sequence representing the sequence for the
            location.
        """
        # Special case for the RemoteLocation
        if isinstance(location, RemoteLocation):
            if location.accession == self.accession:
                first, last = location.get_range()
                return self.sequence[first - 1:last]
            elif sequence is not None:
                return sequence.get_sequence_from_location(location)
        # Get the range and show the sequence according to that
        first, last = location.get_range()
        return self.sequence[first - 1:last]

    def get_complement_sequence(self):
        """ Creates a new Sequence object which represents the
        complement of this sequence. This is only possible for DNA and
        RNA.

        Returns:
            A Sequence object representing the complement of this
            Sequence object.
        Raises:
            ValueError when this sequence represents something else
            than DNA or RNA.
        """
        # Determinant for dna or rna
        dna_valid = match('^[ATCG]+$', self.sequence, IGNORECASE)
        # Check whether this sequence is RNA or DNA
        if not (dna_valid or match('^[AUCG]+$', self.sequence, IGNORECASE)):
            raise ValueError('Sequence is no DNA or RNA')
        # Use Uracil or Thymine dependent on dna_valid
        t = 'T' if dna_valid else 'U'
        reverse = ''
        # Create the complement in reverse
        for char in self.sequence.upper():
            if char == t:
                reverse = 'A' + reverse
            elif char == 'A':
                reverse = t + reverse
            elif char == 'C':
                reverse = 'G' + reverse
            elif char == 'G':
                reverse = 'C' + reverse
        return Sequence(reverse)
