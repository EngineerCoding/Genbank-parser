from os.path import exists
from re import split
from metadata_parser import parse_metadata as parse_actual_metadata
from features_parser import parse_features as parse_actual_features
from origin_parser import parse_origin as parse_actual_origin

CONTINUE_LINE_SPACING = ' ' * 12


class GenbankParser(object):
	def __init__(self, filename):
		if not exists(filename):
			raise ValueError("File {} does not exist.".format(filename))
		self.filehandle = open(filename, 'r')

	def parse_metadata(self, return_meta=True):
		if return_meta:
			return parse_actual_metadata(self)
		self.read_until('FEATURES')
		return True

	def parse_features(self, return_features=True):
		if return_features:
			return parse_actual_features(self)
		self.read_until('ORIGIN')
		return True

	def parse_origin(self, return_origin=True):
		if return_origin:
			return parse_actual_origin(self)
		self.read_until('//')
		return True

	def read_until(self, keyword):
		last_position = 0
		line = ''
		while not line.startswith(keyword) and line:
			last_position = self.filehandle.tell()
			line = self.read_valid_line()
		# If the line evaluates to false, we reached the end of the file
		if not line:
			raise ValueError("Invalid GenBank file")
		# Set the next line to be read to the FEATURES line
		self.filehandle.seek(last_position)

	def read_valid_line(self):
		content = '\n'
		while content == '\n':
			content = self.filehandle.readline()
		return content

	def get_continuing_line(self):
		old_position = self.filehandle.tell()
		line = self.read_valid_line()
		if line and line.startswith(CONTINUE_LINE_SPACING):
			return line
		# Not found, thus return to the reading of said line
		self.filehandle.seek(old_position)

	def handle_keyword(self, keyword, do_split=True, remove_keyword=True,
					   raise_error=True):
		last_position = self.filehandle.tell()
		line = self.read_valid_line().strip()
		if not (line.startswith(keyword) and line):
			if raise_error:
				raise ValueError('Did not find {}'.format(keyword))
			self.filehandle.seek(last_position)
			return
		# Remove the keyword
		if remove_keyword:
			line = line[len(keyword):].strip()
		# Split if necessary
		if do_split:
			return split('\s+', line)
		return line

	def handle_continuing_lines(self, base, delimiter=' '):
		line = self.get_continuing_line()
		while line is not None:
			base += delimiter + line.strip()
			line = self.get_continuing_line()
		return base

	def handle_multiline_keyword(self, keyword, delimiter=' ', **kwargs):
		base = self.handle_keyword(keyword, **kwargs)
		if base:
			return self.handle_continuing_lines(base, delimiter)

	def close(self):
		self.filehandle.close()

	def __exit__(self, *args):
		self.close()
