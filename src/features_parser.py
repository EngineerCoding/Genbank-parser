from re import split

FEATURE_START_SPACE = ' ' * 5


def parse_features(gbp):
	# Check if the header is there
	gbp.handle_keyword('FEATURES', do_split=False, remove_keyword=False)
	features = []
	old_position = gbp.filehandle.tell()
	line = gbp.read_valid_line()
	len_spaces = len(FEATURE_START_SPACE)
	while (line.startswith(FEATURE_START_SPACE) and
			   not line[len_spaces:len_spaces + 1].isspace()):
		name, location = split('\s+', line.strip())
		features.append(Feature(name, location, __parse_attributes(gbp)))
		old_position = gbp.filehandle.tell()
		line = gbp.read_valid_line()
	gbp.filehandle.seek(old_position)
	return features


def __parse_attributes(gbp):
	attributes = {}
	attribute = gbp.handle_keyword('/', do_split=False, raise_error=False)
	while attribute is not None:
		# Parse the key
		key = ''
		for char in attribute:
			if char == '=':
				break
			key += char
		value = attribute[len(key) + 1:]
		if value[0:1] == '"':
			value = __parse_string(gbp, value)
		attributes[key] = value
		attribute = gbp.handle_keyword('/', do_split=False, raise_error=False)
	return attributes


def __parse_string(gbp, current_remaining):
	str = current_remaining[1:]
	while str[-1] != '"':
		str += ' ' + gbp.read_valid_line().strip()
	return str[:-1]


class Feature(object):
	def __init__(self, name, location, attributes):
		self.name = name
		self.location = location
		self.attributes = attributes

	def has_attribute(self, attribute):
		return attribute in self.attributes

	def get_attribute(self, attribute):
		return self.attributes[attribute]
