from enum import Enum


class LocationType(Enum):
    """ All the location types which can occur in a Genbank file """
    single_base = 'single_base'
    adjoining = 'adjoining'
    single_base_range = 'single_base_range'  # Unclear what this is
    range = 'range'
    remote = 'remote'


class AdjoiningLocationType(Enum):
    """ All the types of Adjoining"""
    endonucleolytic = 'endonucleolytic'
    circulair = 'circulair'


class Location:
    """ The base type of location which implement the most basic methods
    for sub classes.
    """

    type = None  # Type of the location

    def __init__(self, location_string):
        """ In the base class implementation this does nothing, but in sub
        classes this method should parse the given location_string.

        Parameters:
            location_string - string
                The string that is going to be parsed
        """
        self.first = -1
        self.second = -1

    def get_range(self):
        """ Retrieves the range of the current location.

        Returns:
            A tuple, containing:
             1. The first location
             2. The second location
        """
        return self.first, self.second

    def _check_not_position(self, position):
        """ This method is used to check if the position is not contained in
        this location.

        Parameters:
            position - Location object
        Returns:
            boolean suc
        """
        return position not in self and isinstance(position, Location)

    def get_diff(self, position):
        """ Calculates the amount of residues the gap is between this location
        and the given location

        Parameters:
            position - Location object
        Returns:
            The difference between this location and the other location. When
            the parameter is not a Location, a 0 is returned.
        """
        if self._check_not_position(position):
            if position.first < self.first:
                return self.first - position.first
            else:  # position.first > self.second
                return position.first - self.second
        return 0

    def is_left(self, position):
        """ Checks if the position is left to this location

        Parameters:
            position - Location object
        Returns:
            Boolean whether it is True or not.
        """
        if self._check_not_position(position):
            return position.first < self.first
        return False

    def is_right(self, position):
        """ Checks if the position is right to this location

        Parameters:
            position - Location object
        Returns:
            Boolean whether it is True or not.
        """
        if self._check_not_position(position):
            return position.first > self.second

    def to_sequence(self, sequence, alt_sequence=None):
        """ Converts this location to the correct position on the given
        sequence object.

        Parameters:
            sequence - Sequence object
                Used to get the string sequence from
            alt_sequence - Sequence object. Default: None
                Should be used when this location represents a RemoteLocation
                object and the first sequence does not have that accession.
        Returns:
            A string sequence representing this location
        """
        return sequence.get_sequence_from_location(self, alt_sequence)

    def __contains__(self, item):
        if isinstance(item, Location):
            return self.first < item.first < self.second
        return False

    def __len__(self):
        return self.second - self.first + 1


class SingleBaseLocation(Location):
    """ Represents a single base, which points to a single number.
    The representation is:
        n
        Where n is a number in a sequence
    """

    type = LocationType.single_base

    def __init__(self, location_string):
        super().__init__(location_string)
        self.first = self.second = _convert(int, location_string)

    def __str__(self):
        return str(self.first)


class DelimitedLocation(Location):
    """ A DelimitedLocation is a sub class which can be used to easily
    parse locations like:
        42..56
        4^6
        etc.
    Only the delimiter on the class has to be set on a sub class of
    this class
    """

    delimiter = None  # The delimiter to use

    def __init__(self, location_string):
        """ Parses the location string using the set delimiter. When
        the class does not define one, a ValueError will be risen.
        Otherwise, it will create a list with all the values using
        the delimiter as separator.

        This class defines 3 methods which should be defined in sub
        classes:
            1. _parse_left: parses the most left item of the list
            2. _parse_in_between: parses everything but the most left
               and most right element in the list.
            3. _parse_right: parses the most right item of the list
        """
        super().__init__(location_string)
        # check the delimiter
        if not self.delimiter:
            raise ValueError('Delimiter has to be set!')
        # Create the list
        split = location_string.split(self.delimiter)
        # Call the correct parse methods with the correct elements
        self._parse_left(split[0])
        if len(split) == 1:
            raise ValueError('Expected values on both sides of ' +
                             self.delimiter)
        for i in range(1, len(split) - 1):
            self._parse_in_between(split[i])
        self._parse_right(split[len(split) - 1])

    def _parse_left(self, string):
        raise NotImplementedError

    def _parse_right(self, string):
        raise NotImplementedError

    def _parse_in_between(self, string):
        raise ValueError('Too many ' + self.delimiter + ' have been used!')

    def __str__(self):
        return str(self.first) + self.delimiter + str(self.second)


class AdjoiningLocation(DelimitedLocation):
    """ Parses an Adjoining location which can look like:
            n^n+1 -> endonucleolytic cleavage
            n^1 -> circulair molecule
    """

    type = LocationType.adjoining
    delimiter = '^'

    subtype = ''

    def _parse_left(self, string):
        self.first = _convert(int, string)

    def _parse_right(self, string):
        self.second = _convert(int, string)
        # Set the correct type
        if self.second == 1:
            self.subtype = AdjoiningLocationType.circulair
        elif self.second - 1 == self.first:
            self.subtype = AdjoiningLocationType.endonucleolytic
        else:
            raise ValueError('Invalid adjoining location: {}^{}'
                             .format(self.first, self.second))

    def __len__(self):
        return 2


class RangeLocation(DelimitedLocation):
    """ Parses a range location which looks like:
         x..y
            Where x and y or both integers and y is greater than x
         <x..y
            Same as above, but x can be smaller than the actual
            defined value
         x..>y
            Same as above, but y can be greater than the actual
            defined value
    This represents a range of residues including both x and y.
    """

    type = LocationType.range
    delimiter = '..'

    def __init__(self, string):
        super().__init__(string)
        self.can_be_lesser = False
        self.can_be_greater = False

    def _parse_left(self, string):
        self.can_be_lesser = string[0] == '<'
        if self.can_be_lesser:
            string = string[1:]
        self.first = _convert(int, string)

    def _parse_right(self, string):
        self.can_be_greater = string[0] == '>'
        if self.can_be_greater:
            string = string[1:]
        self.second = _convert(int, string)

    def __str__(self):
        return '{}..{}'.format(self.first, self.second)


class RemoteLocation(DelimitedLocation):
    """ Parses a RemoteLocation which looks like this:
         accession:location
            Where accession is a regular accesion number and location is
            an actual type of Location. This means to look in the given
            accession sequence at that given location
    """

    type = LocationType.remote
    delimiter = ':'

    accession = ''
    location = None

    def _parse_left(self, string):
        self.accession = string.strip()
        if len(string) == 0:
            raise ValueError('Expected accession string!')

    def _parse_right(self, string):
        self.location = parse_location(string)

    def get_range(self):
        return self.location.get_range()

    def get_accession(self):
        return self.accession

    def to_sequence(self, sequence, alt_sequence=None):
        return self.location.to_sequence(sequence, alt_sequence)

    def __contains__(self, item):
        if isinstance(item, RemoteLocation) and \
                        item.accession == self.accession:
            return item in self.location
        return False

    def __len__(self):
        return len(self.location)

    def __str__(self):
        return self.accession + ':' + str(self.location)


class JoinedLocation(Location):
    """ Representation of a joined operator in the location object.
    Looks like:
        joined(range, range)
         Where range represents a RangeLocation
    """

    def __init__(self, *locations):
        super().__init__(None)
        self.locations = locations

    def calculate_inversed_locations(self, genome_length):
        """ Calculates the locations that are not represented by this
        JoinedLocation object. Those regions can also be referred to,
        when residues are DNA or RNA, introns.

        Parameters:
            genome_length - int
                The length of the sequence to compare to
        Returns:
            A list of RangeLocation objects representing the inversed
            locations of this object.
        """
        # Length must be positive
        if genome_length < 1:
            raise ValueError('Genome length must be set!')
        last_seq_index = 1
        locations = []

        for first, last in self.get_ranges():
            # When the last border is the same as the first of this
            # location, there is no intron in front of here
            if first == last_seq_index:
                continue
            locations.append(RangeLocation('{}..{}'.format(last_seq_index,
                                                           first - 1)))
            last_seq_index = last + 1
        # Check for remaining intron to the right of the last exon
        if genome_length > last_seq_index:
            locations.append(RangeLocation('{}..{}'.format(last_seq_index,
                                                           genome_length)))
        return locations

    def get_range(self):
        raise NotImplementedError

    def get_ranges(self):
        """ This is a generator object which returns all the ranges
        of the locations which this JoinedLocation object contains.

        Returns:
            A generator which yields the range of a location
        """
        for location in self.locations:
            yield location.get_range()

    def to_sequence(self, sequence, alt_sequence=None):
        generated_sequence = ''
        for location in self.locations:
            generated_sequence += location.to_sequence(sequence, alt_sequence)
        return generated_sequence


class ComplementLocation(JoinedLocation):
    def __init__(self, location):
        super().__init__(location)

    def get_translated_joined(self, genome_length):
        new_locations = []
        for location in self.locations:
            first = genome_length - location.second + 1
            second = location.second - location.first + first
            new_locations.append(RangeLocation('{}..{}'.format(first, second)))
        return JoinedLocation(*new_locations)

    def get_complement_sequences(self, sequence):
        sequences = []
        locations = self.get_translated_joined(sequence.lenght()).locations
        complement_sequence = sequence.get_complement_sequence()
        for location in locations:
           sequences.append(complement_sequence
                            .get_sequence_from_location(location))
        return sequences


def _convert(var_type, string):
    try:
        return var_type(string)
    except ValueError:
        raise ValueError('Expected {}: {}'.format(var_type.__class__.__name__,
                                                  string))


def parse_location(location_string):
    """ The main parser function which parses a string to a location
    object.

    Parameters:
        location_string - string
            The string to parse to a Location object
    """
    location = __parse_string_arguments(location_string.strip())
    if len(location) > 1:
        raise ValueError('Cannot parse {} to a Location!'
                         .format(location_string))
    # Check for a function
    result = __execute_function(location[0])
    # If not a function, a regular location string
    return result or __parse_location_string(location[0])


def __execute_function(argument, top_level=True):
    """ Executes an argument string which represents a function.

    Parameters:
        argument - string
            The function to parse and execute
        top_level - boolean
            A boolean to determine whether functions are allowed at
            this level. Operators cannot nest deep in locations.
    Returns:
        The return value of the executed function (which is most
        likely a JoinedLocation or ComplementLocation)
    """
    for func_map in function_mapping:
        # check whether it is this function
        if argument.startswith(func_map['name']):
            # Parse the arguments of the functions
            remaining = argument[len(func_map['name']):].lstrip()
            if remaining[0] != '(':
                raise ValueError('Expected \'()\'')
            remaining = remaining[1:]
            arguments = __parse_string_arguments(remaining)
            # Check if the function can handle multiple arguments
            if func_map['one_arg'] and len(arguments) > 1:
                raise ValueError('1 argument allowed for ' + func_map['name'])
            # Check if we can execute another functions
            if top_level:
                for i in range(len(arguments)):
                    result = __execute_function(arguments[i], False)
                    if result is not None:
                        arguments[i] = result
            # Parse remaining string values in the arguments to
            # Location objects
            for i in range(len(arguments)):
                if isinstance(arguments[i], str):
                    arguments[i] = __parse_location_string(arguments[i])
            return func_map['function'](*arguments)


def __parse_string_arguments(location_string):
    """ Parses the argument for a function. Note that the first
    parenthesis should be removed to avoid going up a parsing
    level. This parsing level defines when to stop parsing.

    Parameters:
        location_string - string
            The string to parse from
    Returns:
        A list of arguments (string)
    """
    argument = ''
    arguments = []
    parentheses_level = 0
    for char in location_string:
        if char == '(':
            parentheses_level += 1
            argument += char
        elif char == ')' and parentheses_level != 0:
            parentheses_level -= 1
            argument += char
        elif char == ')' and parentheses_level == 0:
            break
        elif char == ',' and parentheses_level == 0:
            arguments.append(argument.strip())
            argument = ''
        else:
            argument += char
    arguments.append(argument.strip())
    return arguments


def __parse_location_string(location_string):
    """ Selects the correct parses class """
    for (char, cls) in location_mapping:
        if char in location_string:
            return cls(location_string)
    return SingleBaseLocation(location_string)


location_mapping = ((':', RemoteLocation), ('..', RangeLocation),
                    ('^', AdjoiningLocation))


def __create_func_mapping(name, func, one_arg=False, sub_func=False,
                          contain_sub_func=False):
    return dict(name=name, function=func, one_arg=one_arg, sub_func=sub_func,
                contain_sub_func=contain_sub_func)


# Definition of the operators/functions
function_mapping = [
    __create_func_mapping('complement', ComplementLocation, True, False,
                          True),
    __create_func_mapping('join', JoinedLocation, False, True, False),
    # __create_func_mapping('order', __parse_order, False, True, False)
]
