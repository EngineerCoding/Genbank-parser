def parse_metadata(gbp):
    """ The main method which parses the metadata to a Metadata object.

    Parameters:
        gbp - GenbankParser object
            The parser which holds the file pointer of the genbank
            file.
    Returns:
        A Metadata object
    """
    # The order to parse
    part_argument_gens = [__parse_locus, __parse_definition,
                          __parse_accession_version_keywords,
                          __parse_source, __parse_publications,
                          __parse_comment]
    args = []
    # Execute each parsing part and create a Metadata object with
    # the collected return values.
    for argument_gen in part_argument_gens:
        args += argument_gen(gbp)
    return Metadata(*args)


def __parse_locus(gbp):
    """ Parses the LOCUS tag

    Parameters:
        gbp - GenbankParser object
            The parser which holds the file pointer of the genbank
            file.
    Returns:
        A list with the following values:
         1. Locus name - string
         2. Sequence length - string
         3. The unit which defines the length (usually bp) - string
         4. molecule type - string
         5. Genbank division - string
         6. Modification date - string
    """
    parts = gbp.handle_keyword('LOCUS')
    # Delete unnecessary parts
    del parts[2]  # Bp
    # Check if we have a type
    if len(parts) == 5:
        parts.insert(3, '')
    return parts


def __parse_definition(gbp):
    """ Parses a definition which is more like a description.

    Parameters:
        gbp - GenbankParser object
            The parser which holds the file pointer of the genbank
            file.
    Returns:
        A list of length 1 with the description
    """
    return [gbp.handle_multiline_keyword('DEFINITION', do_split=False)]


def __parse_accession_version_keywords(gbp):
    """ Parses the ACCESSION, VERSION and KEYWORDS part.

    Parameters:
        gbp - GenbankParser object
            The parser which holds the file pointer of the genbank
            file.
    Returns:
        A list with the accessions (which can be multiple in a string),
        a tuple with the version information and a string of keywords
    """
    accession = gbp.handle_keyword('ACCESSION', do_split=False)
    versions = tuple(gbp.handle_keyword('VERSION'))
    # only eat the dblink when it is available, currently not supported
    gbp.handle_keyword('DBLINK', do_split=False, remove_keyword=False,
                       raise_error=False)
    keywords = gbp.handle_multiline_keyword('KEYWORDS', do_split=False)
    return [accession, versions, keywords]


def __parse_source(gbp):
    """ Parses the SOURCE and ORGANISM keywords

    Parameters:
        gbp - GenbankParser object
            The parser which holds the file pointer of the genbank
            file.
    Returns:
        A list with the organism string and complete string of the ORGANISM
        keyword
    """
    source = gbp.handle_keyword('SOURCE', do_split=False)
    organism = gbp.handle_multiline_keyword('ORGANISM', delimiter='\n',
                                            do_split=False)
    return [source, organism]


def __parse_publications(gbp):
    """ Parses the publications of this GenBank file.

    Parameters:
        gbp - GenbankParser object
            The parser which holds the file pointer of the genbank
            file.
    Returns:
        A list in a list which contains all the publications
    """
    publications = []
    # Initialise the parsing of references; publications
    reference = gbp.handle_multiline_keyword('REFERENCE', do_split=False,
                                             raise_error=False)
    # Keep parsing all references
    while reference is not None:
        # Parse the data of a publication
        authors = gbp.handle_multiline_keyword('AUTHORS', do_split=False)
        title = gbp.handle_multiline_keyword('TITLE', do_split=False)
        journal = gbp.handle_multiline_keyword('JOURNAL', do_split=False)
        # PUBMED is optional
        pubmed = gbp.handle_multiline_keyword('PUBMED', do_split=False,
                                              raise_error=False)
        publications.append(Publication(reference, authors, title, journal,
                                        pubmed))
        reference = gbp.handle_multiline_keyword('REFERENCE', do_split=False,
                                                 raise_error=False)
    return [publications]


def __parse_comment(gbp):
    """ Parses the comment and PRIMARY keyword of this GenBank file.

        Parameters:
            gbp - GenbankParser object
                The parser which holds the file pointer of the genbank
                file.
        Returns:
            An empty list, we are currently not interested in the comments
        """
    # Only eat the comment when it is available
    gbp.handle_multiline_keyword('COMMENT', do_split=False,
                                 remove_keyword=False, raise_error=False)
    gbp.handle_multiline_keyword('PRIMARY', do_split=False,
                                 remove_keyword=False, raise_error=False)
    return []


class Metadata(object):
    def __init__(self, locus_name, seq_length, molecule_type, formation,
                 gb_division, modification_date, description, accession,
                 version, keywords, source, organism,
                 publications):
        self.locus_name = locus_name
        self.seq_length = int(seq_length)
        self.molecule_type = molecule_type
        self.division = gb_division
        self.molecule_formation = formation
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
