from re import split

from .location_parser import parse_location

FEATURE_START_SPACE = ' ' * 5


def parse_features(gbp):
    """ The main method which parses the FEATURES to a list with
    Feature objects.
    Parameters:
        gbp - GenbankParser object
            The parser which holds the file pointer of the genbank
            file.
    Returns:
        A list of Feature objects
    """
    # The FEATURES line must be the next line
    gbp.handle_keyword('FEATURES', do_split=False, remove_keyword=False)
    features = []  # The list to fill
    # The old position before the line in case there is no feature
    # available anymore
    old_position = gbp.filehandle.tell()
    line = gbp.read_valid_line()  # The next line
    len_spaces = len(FEATURE_START_SPACE)
    # Keep checking whether the file has the required spaces to be a
    # Feature
    while ((line.startswith(FEATURE_START_SPACE) and
            not line[len_spaces:len_spaces + 1].isspace())):
        # Get the name and location string
        name, location = split('\s+', line.strip())
        # Start parsing the attributes of this Feature
        # Create a Feature object and append it to the list
        features.append(Feature(name, location, __parse_attributes(gbp)))
        # Read the next line
        old_position = gbp.filehandle.tell()
        line = gbp.read_valid_line()
    # A line without feature has been hit, set the file pointer to
    # before that happening
    gbp.filehandle.seek(old_position)
    return features


def __parse_attributes(gbp):
    """ This method will parse the attributes of a Feature.

    Parameters:
        gbp - GenbankParser object
            The parser which holds the file pointer of the genbank
            file.
    Returns:
        A dictionary of attributes where the values all are strings
    """
    attributes = {}
    # Check if this is an attribute
    attribute = gbp.handle_keyword('/', do_split=False, raise_error=False)
    while attribute is not None:
        # Parse the key
        key = ''
        for char in attribute:
            if char == '=':
                break
            key += char
        # Parse the value
        value = attribute[len(key) + 1:]
        # When the value is a string, parse it as a string (which can be
        # multiline)
        if value[0:1] == '"':
            value = __parse_string(gbp, value)
        attributes[key] = value
        # Try for a next attribute
        attribute = gbp.handle_keyword('/', do_split=False, raise_error=False)
    return attributes


def __parse_string(gbp, value):
    """ Parses a string over multiple lines """
    remaining = value[1:]
    # Keep reading until a " has been hit
    # TODO: a quote can be escaped in Genbank files
    while remaining[-1] != '"':
        remaining += ' ' + gbp.read_valid_line().strip()
    return remaining[:-1]


class Feature(object):
    """ A Feature object has three things:
    1. A name of the Feature
    2. A location string
    3. A dictionary with attributes
    """

    def __init__(self, name, location, attributes):
        """ This constructor will parse the location to a Location
        object, which is more useful then the simple location string.

        Parameters:
            name - string
                The name of the Feature
            location - string
                The string which represents a location
            attributes - dict
                A dictionary full of attributes wich are related to this
                Feature.
        """
        self.name = name
        self.location = parse_location(location)
        self.attributes = attributes

    def has_attribute(self, attribute):
        """ Checks whether an attribute aexists or not

        Parameters:
            attribute - string
                The name of attribute which will be checked for
        Returns:
            A boolean representing whether the attribute is in the
            dictionary.
        """
        return attribute in self.attributes

    def get_attribute(self, attribute):
        """ Retrieves an attribute from the attributes dictionary

        Parameters:
            attribute - string
                The name of attribute which will be checked for
        Returns:
            A string representing the attribute as a literal from the
            genbank file. An exception for this is a string in the
            genbank file, as those are parsed off.
        Raises:
            KeyError when the attribute does not exist in the attributes
            dictionary.
        """
        return self.attributes[attribute]
