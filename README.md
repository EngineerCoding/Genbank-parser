# Genbank-parser

This is a very minimal Genbank file parser, which is not really optimised for various Genbank files as this is kind of hardcoded. Plans are being made to upgrade this parser to do dynamic parsing, however currently this will do for a project. 

This parser integrates a location parser which can be given to a Sequence object, which then retrieves the correct sequence. It is recommend though, to use the direct methodes of the location because a `JoinLocation` is for instance a wrapper for other `Location` objects.

Current usage:
```
from genbank_parser import GenbankParser
with GenbankParser('myfile.gb') as parser:
  # The following methods have to be called to properly parse a Genbank file
  parser.parse_metadata() # Optional to actually store this data
  parser.parse_features() # Optional to actually store this data
  parser.parse_origin() # Optional to actually store this data
```

Currently there is about no documentation, so code has to be read to understand what it does. Documentation is the current priority though.
