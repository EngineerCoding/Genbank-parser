from enum import Enum


class LocationType(Enum):
    single_base = 'single_base'
    adjoining = 'adjoining'
    single_base_range = 'single_base_range'  # Unclear what this is
    range = 'range'
    remote = 'remote'


class AdjoiningLocationType(Enum):
    endonucleolytic = 'endonucleolytic'
    circulair = 'circulair'


class Location:
    type = None

    def __init__(self, location_string):
        self.first = -1
        self.second = -1

    def get_range(self):
        return self.first, self.second

    def _check_position(self, position):
        return position not in self and isinstance(position, Location)

    def get_diff(self, position):
        if self._check_position(position):
            if position.first < self.first:
                return self.first - position.first
            else:  # position.first > self.second
                return position.first - self.second
        return 0

    def is_left(self, position):
        if self._check_position(position):
            return position.first < self.first
        return False

    def is_right(self, position):
        if self._check_position(position):
            return position.first > self.second

    def to_sequence(self, sequence, alt_sequence=None):
        return sequence.get_sequence_from_location(self, alt_sequence)

    def __contains__(self, item):
        if isinstance(item, Location):
            return self.first < item.first < self.second
        return False

    def __len__(self):
        return self.second - self.first + 1


class SingleBaseLocation(Location):
    type = LocationType.single_base

    def __init__(self, location_string):
        super().__init__(location_string)
        self.first = self.second = _convert(int, location_string)

    def __str__(self):
        return str(self.first)


class DelimitedLocation(Location):
    delimiter = None

    def __init__(self, location_string):
        super().__init__(location_string)
        if not self.delimiter:
            raise ValueError('Delimiter has to be set!')
        splitted = location_string.split(self.delimiter)
        self._parse_left(splitted[0])
        if len(splitted) == 1:
            raise ValueError('Expected values on both sides of ' +
                             self.delimiter)
        for i in range(1, len(splitted) - 1):
            self._parse_in_between(splitted[i])
        self._parse_right(splitted[len(splitted) - 1])

    def _parse_left(self, string):
        raise NotImplementedError

    def _parse_right(self, string):
        raise NotImplementedError

    def _parse_in_between(self, string):
        raise ValueError('Too many ' + self.delimiter + ' have been used!')

    def __str__(self):
        return str(self.first) + self.delimiter + str(self.second)


class AdjoiningLocation(DelimitedLocation):
    type = LocationType.adjoining
    delimiter = '^'

    subtype = ''

    def _parse_left(self, string):
        self.first = _convert(int, string)

    def _parse_right(self, string):
        self.second = _convert(int, string)
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
    # Not supported <x..y and x..>y yet

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
    def __init__(self, *locations):
        super().__init__(None)
        self.locations = locations

    def calculate_inversed_locations(self, genome_length):
        if genome_length < 1:
            raise ValueError('Genome length must be set!')
        last_seq_index = 1
        locations = []
        for first, last in self.get_ranges():
            if first == last_seq_index:
                continue
            locations.append(RangeLocation('{}..{}'.format(last_seq_index,
                                                           first - 1)))
            last_seq_index = last + 1
        if genome_length > last_seq_index:
            locations.append(RangeLocation('{}..{}'.format(last_seq_index,
                                                           genome_length)))
        return locations

    def get_range(self):
        raise NotImplementedError

    def get_ranges(self):
        for location in self.locations:
            yield location.get_range()

    def to_sequence(self, sequence, alt_sequence=None):
        generated_sequence = ''
        for location in self.locations:
            generated_sequence += location.to_sequence(sequence, alt_sequence)
        return generated_sequence


class ComplementLocation(RangeLocation):
    def __init__(self, range_location):
        super().__init__(str(range_location))
        self._5_3_sequence = range_location

    def calculate_3_5_positions(self, genome_length):
        self.first = genome_length - self._5_3_sequence.second + 1
        self.second = (self._5_3_sequence.second -
                       self._5_3_sequence.first + self.first)
        return RangeLocation('{}..{}'.format(self.first, self.second))

    def get_complement(self, sequence):
        complement_location = self.calculate_3_5_positions(sequence.length())
        return (sequence.get_complement_sequence()
                .get_sequence_from_location(complement_location))


def _convert(var_type, string):
    try:
        return var_type(string)
    except ValueError:
        raise ValueError('Expected {}: {}'.format(var_type.__class__.__name__,
                                                  string))


def parse_location(location_string):
    location = __parse_string_arguments(location_string.strip())
    if len(location) > 1:
        raise ValueError('Cannot parse {} to a Location!'
                         .format(location_string))
    # Check for a function
    result = __execute_function(location[0])
    return result or __parse_location_string(location[0])


def __execute_function(argument, top_level=True):
    for func_map in function_mapping:
        if argument.startswith(func_map['name']):
            remaining = argument[len(func_map['name']):].lstrip()
            if remaining[0] != '(':
                raise ValueError('Expected \'()\'')
            remaining = remaining[1:]
            arguments = __parse_string_arguments(remaining)
            if func_map['one_arg'] and len(arguments) > 1:
                raise ValueError('1 argument allowed for ' + func_map['name'])
            if top_level:
                for i in range(len(arguments)):
                    result = __execute_function(arguments[i], False)
                    if result is not None:
                        arguments[i] = result
            for i in range(len(arguments)):
                if isinstance(arguments[i], str):
                    arguments[i] = __parse_location_string(arguments[i])
            return func_map['function'](*arguments)


def __parse_string_arguments(location_string):
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


function_mapping = [
    __create_func_mapping('complement', ComplementLocation, True, False,
                          True),
    __create_func_mapping('join', JoinedLocation, False, True, False),
    # __create_func_mapping('order', __parse_order, False, True, False)
]
