from os.path import exists
from re import split

from .features_parser import parse_features as parse_actual_features
from .metadata_parser import parse_metadata as parse_actual_metadata
from .origin_parser import parse_origin as parse_actual_origin

CONTINUE_LINE_SPACING = ' ' * 12


class GenbankParser(object):
    """ This class is used to parse a Genbank file.

    The parsing of a file contains 3 stages:
      1. Parsing of the metadata (Data of the sequences such as length,
         name, what organism and the publications.)
      2. Parsing of the features. These are the main data which most
         likely will be used for scientific research. This contains all
         fields (for instance source or CDS) with all of their
         attributes.
      3. Parsing of the origin, or better known as the sequence where
         the whole Genbank file is based on.

    These parsing stages should be parsed in this sequence. The option
    exists for the stage to skip the data (and thus not load it into
    memory) but have the file pointer set at the correct position for
    the next stage.

    Note that this object can be considered as file: it needs to be
    closed which is supported with a 'with' statement are a regular
    'close' method.
    """

    def __init__(self, filename):
        """ Creates a new parser from the given file.

        Parameters:
            filename - string
                The name of the file pointing to the file which needs
                to be parsed.
        Raises:
            ValueError when the file does not exist on the filesystem.
        """
        if not exists(filename):
            raise ValueError('File {} does not exist.'.format(filename))
        self.filehandle = open(filename, 'r')

    def parse_metadata(self, return_meta=True):
        """ Parses the metadata as described in the docstring of this
        class.

        Parameters:
            return_meta - boolean. Default: True
                This is a boolean which determines whether to store the
                parsed data or not. True for storing data, False for
                not storing the data.
        Return:
            When return_meta is set to True, this will return a Metadata
            object. If set to False, this will simply return True.
        """
        if return_meta:
            return parse_actual_metadata(self)
        # Read until the FEATURES are hit
        self.read_until('FEATURES')
        return True

    def parse_features(self, return_features=True):
        """ Parses the features as described in the docstring of this
        class.

        Parameters:
            return_features - boolean. Default: True
                This is a boolean which determines whether to store the
                parsed data or not. True for storing data, False for
                not storing the data.
        Return:
            When return_features is set to True, this will return a
            list of Feature objects. If set to False, this will
            simply return True.
        """
        if return_features:
            return parse_actual_features(self)
        # Read until the ORIGIN is hit
        self.read_until('ORIGIN')
        return True

    def parse_origin(self, return_origin=True):
        """ Parses the origin as described in the docstring of this
        class.

        Parameters:
            return_origin - boolean. Default: True
                This is a boolean which determines whether to store the
                parsed data or not. True for storing data, False for
                not storing the data.
        Return:
            When return_origin is set to True, this will return a
            Sequence object. If set to False, this will simply
            return True.
        """
        if return_origin:
            return parse_actual_origin(self)
        # Read until the end of the file has been reached
        # Might be over eager, as this is actually pointless
        self.read_until('//')
        return True

    def read_until(self, keyword):
        """ Reads until a keyword has been hit. When this keyword is
        hit, it will set the file pointer back to right before the
        keyword, so the parser can double check whether the keyword
        is even there.

        Parameters:
            keyword - string
                This is the string to look for at the begin of a line.
                Whitespace will be stripped away, so it is not
                necessary to account for that.

        """
        # A variable to store the last line in
        line = ' '
        # A varibale to store the position in before the last read line
        last_position = 0
        # Keep checking until we have line which starts with the
        # keyword
        while not line.startswith(keyword) and line:
            # Still not found, so read a line and store the position
            # before this line
            last_position = self.filehandle.tell()
            line = self.read_valid_line().strip()
        # If the line evaluates to false, the end of the file has
        # been reached
        if not line:
            raise ValueError('Invalid GenBank file')
        # Set the next line to be read to the line with the keyword
        self.filehandle.seek(last_position)

    def read_valid_line(self):
        """ Keeps reading a line until a line is not an empty line.
         This method will return the last line which is not equal to
         '\n' since this is an empty line as defined by python.
         Note that this method can return an empty string (without
         '\n') which means that the end of file has been reached.
        """
        content = '\n'
        while content == '\n':
            content = self.filehandle.readline()
            # Check if the line is full of whitespace
            if len(content) > 1 and len(content.strip()) == 0:
                content = '\n'
        return content

    def get_continuing_line(self):
        """ In genbank a line can continue on the next line, however
        the line will then be preceded with 12 spaces. This method
        will try to find that first occuring line.

        Returns:
            A string which represents the line when it is a so-called
            continued line. Returns None when no line has been found.
        """
        # Save the position in case we don't find a valid line
        old_position = self.filehandle.tell()
        line = self.read_valid_line()
        if line and line.startswith(CONTINUE_LINE_SPACING):
            return line
        # Not found, thus return to the reading of said line
        self.filehandle.seek(old_position)

    def handle_keyword(self, keyword, do_split=True, remove_keyword=True,
                       raise_error=True):
        """ Handles a keyword of the Genbank specification which also
        can automatically do some basic processing.

        Parameters:
            keyword - string
                The string which a line should start with to be
                applicable for that keyword.
            do_split - boolean. Default: True
                This determines whether the contents should be split by
                whitespace. When set to True, this will be done.
            remove_keyword - boolean. Default: True
                This determines whether to remove the keyword from the
                string, so the data of that line can be accessed right
                away. When set to True, this will be done.
            raise_error - boolean. Default: True
                This determines whether to error if the keyword has not
                been found. When set to True, this will happen.
        Returns:
            When do_split has been set to True this will return a list
            of items which represent the data of the line. No
            whitespace is included in this setting.
            When do_split is not set, the method will return the read
            string.
        """
        # Store the last position in case we need to return to before
        # this line.
        last_position = self.filehandle.tell()
        # Read the line
        line = self.read_valid_line().strip()
        if not (line.startswith(keyword) and line):
            # Raise an error dependent of the raise_error parameter
            if raise_error:
                raise ValueError('Did not find {}'.format(keyword))
            # Set the file pointer back to before reading this line
            self.filehandle.seek(last_position)
            return
        # Remove the keyword dependent on the remove_keyword parameter
        if remove_keyword:
            line = line[len(keyword):].strip()
        # Split on whitespace dependent on the do_split parameter
        if do_split:
            return split('\s+', line)
        return line

    def handle_continuing_lines(self, base, delimiter=' '):
        """ Reads multiple continuing lines automatically.

        Parameters:
            base - string
                The base string to append the line to
            delimiter - string. Default: ' '
                A delimiter which is appended in between lines.
        Returns:
            The base string with the available lines and delimiters
            appended to it.
        """
        line = self.get_continuing_line()
        while line is not None:
            base += delimiter + line.strip()
            line = self.get_continuing_line()
        return base

    def handle_multiline_keyword(self, keyword, delimiter=' ', **kwargs):
        """ Handles a keyword which may need multiple strings.
        Basically the same as the 'handle_keyword' method but then with
        multiple lines automatically adjoined to them.

        Parameters:
            keyword - string
                The string which a line should start with to be
                applicable for that keyword.
            delimiter - string
                A delimiter which is appended in between lines.
            remove_keyword - boolean. Default: None
                This determines whether to remove the keyword from the
                string, so the data of that line can be accessed right
                away. When set to True, this will be done.
            raise_error - boolean. Default: None
                This determines whether to error if the keyword has not
                been found. When set to True, this will happen.

        Returns:
            When do_split has been set to True this will return a list
            of items which represent the data of the line. No
            whitespace is included in this setting.
            When do_split is not set, the method will return the read
            string.

            The line with the keyword with the available lines and
            delimiters appended to it.
        """
        base = self.handle_keyword(keyword, **kwargs)
        if base:
            return self.handle_continuing_lines(base, delimiter)

    def close(self):
        """ Closes the file handle """
        self.filehandle.close()

    def __enter__(self):
        return self

    def __exit__(self, *args):
        self.close()
