from re import match, split, IGNORECASE

from location_parser import RemoteLocation


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
		self.accession = None

	def get_sequence(self):
		return self.sequence

	def set_accession(self, accession):
		self.accession = accession

	def get_accession(self):
		return self.accession

	def length(self):
		return len(self.sequence)

	def get_sequence_from_location(self, location, sequence=None):
		if isinstance(location, RemoteLocation):
			if location.accession == self.accession:
				first, last = location.get_range()
				return self.sequence[first - 1:last]
			elif sequence is not None:
				return sequence.get_sequence_from_location(location)
		first, last = location.get_range()
		return self.sequence[first - 1:last]

	def get_complement_sequence(self):
		dna_valid = match('^[ATCG]+$', self.sequence, IGNORECASE)
		if not (dna_valid or match('^[AUCG]+$', self.sequence, IGNORECASE)):
			raise ValueError("Sequence is no DNA or RNA")
		t = 'T' if dna_valid else 'U'
		reverse = ''
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
