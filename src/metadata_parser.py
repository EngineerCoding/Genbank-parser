def parse_metadata(gbp):
	part_argument_gens = [__parse_locus, __parse_definition,
						  __parse_accession_version_keywords,
						  __parse_source, __parse_publications,
						  __parse_comment]
	args = []
	for argument_gen in part_argument_gens:
		args += argument_gen(gbp)
	return Metadata(*args)


def __parse_locus(gbp):
	parts = gbp.handle_keyword('LOCUS')
	# Delete unnecessary parts
	del parts[2]  # Bp
	# Check if we have a type
	if len(parts) == 5:
		parts.insert(3, '')
	return parts


def __parse_definition(gbp):
	return [gbp.handle_multiline_keyword('DEFINITION', do_split=False)]


def __parse_accession_version_keywords(gbp):
	accession = gbp.handle_keyword('ACCESSION')
	versions = tuple(gbp.handle_keyword('VERSION'))
	# only eat the dblink when it is available
	gbp.handle_keyword('DBLINK', do_split=False, remove_keyword=False,
					   raise_error=False)
	keywords = gbp.handle_multiline_keyword('KEYWORDS', do_split=False)
	return accession + [versions] + [keywords]


def __parse_source(gbp):
	source = gbp.handle_keyword('SOURCE', do_split=False)
	organism = gbp.handle_multiline_keyword('ORGANISM', delimiter='\n',
											do_split=False)
	return [source, organism]


def __parse_publications(gbp):
	publications = []
	reference = gbp.handle_multiline_keyword('REFERENCE', do_split=False,
											 raise_error=False)
	while reference is not None:
		authors = gbp.handle_multiline_keyword('AUTHORS', do_split=False)
		title = gbp.handle_multiline_keyword('TITLE', do_split=False)
		journal = gbp.handle_multiline_keyword('JOURNAL', do_split=False)
		pubmed = gbp.handle_multiline_keyword('PUBMED', do_split=False,
											  raise_error=False)
		publications.append(Publication(reference, authors, title, journal,
										pubmed))
		reference = gbp.handle_multiline_keyword('REFERENCE', do_split=False,
												 raise_error=False)
	return [publications]


def __parse_comment(gbp):
	# Only eat the comment when it is available
	gbp.handle_multiline_keyword('COMMENT', do_split=False,
								 remove_keyword=False, raise_error=False)
	return []


class Metadata(object):
	def __init__(self, locus_name, seq_length, molecule_type, type,
				 gb_division, modification_date, description, accession,
				 version, keywords, source, organism,
				 publications):
		self.locus_name = locus_name
		self.seq_length = int(seq_length)
		self.molecule_type = molecule_type
		self.molecule_formation = type
		self.division = gb_division
		self.modification_date_str = modification_date
		self.description = description
		self.accession = accession
		self.version = version
		self.keywords = keywords
		self.source = source
		self.organism = organism
		self.publications = publications


class Publication(object):
	def __init__(self, reference, authors, title, journal, pubmed):
		self.reference = reference
		self.authors = authors
		self.title = title
		self.journal = journal
		self.pubmed = pubmed
