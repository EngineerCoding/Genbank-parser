from re import match, split


def parse_origin(gbp):
	# Check if the header is there
	gbp.handle_keyword('ORIGIN', do_split=False, remove_keyword=False)
	sequence = ''
	line = gbp.read_valid_line().strip()
	while match('^\d+.*', line):
		splitted = split('\s+', line)
		del splitted[0]
		sequence += ''.join(splitted).upper()
		line = gbp.read_valid_line().strip()
	return Sequence(sequence)


class Sequence(object):
	def __init__(self, sequence):
		self.sequence = sequence
