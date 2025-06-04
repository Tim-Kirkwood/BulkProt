# Motivation
UniProt is a great database, but there limited options for bulk searching.  Users can download a lot of entries if they have accessions for specific entries already in hand, and they can perform lots of flexible searches if they can code, but non-coders have no tools for bulk searching.  

BulkProt fills this gap by allowing users to specify a list of search terms in a CSV input file, and presents the output of multiple searches as a single consolidated CSV for easy analysis.  It also allows users to systematically search for entries that are likely to be equivalent to their search term, but are hard to detect via traditional searches due to inconsistent naming conventions.

# Approach



# Tutorial
```
usage: BulkProt [-h] -csv [-e] [-q] [-sd] [-f] [-o]

options:
  -h, --help            show this help message and exit
  -csv, --csv_path  (REQUIRED) Full filepath to the CSV-format file containing your queries. (default: None)
  -e, --excel_compatible
                        If set, and you include a sequence field, then sequence will be truncated to 32000 chars (the
                        maximum number of chars in an excel cell). (default: False)
  -q, --quick           If set, will run much faster, but a single query within your CSV of queries can only return up
                        to 500 results. (default: False)
  -sd, --seed_only      If set, will only perform a seed search. (default: False)
  -f, --fields      (OPTIONAL) Result fields you would like to store as columns in the output CSV. Full list of
                        available entries is at https://www.uniprot.org/help/return_fields ("Returned Field"). Input
                        fields should be seperated by a comma with no spaces (e.g. "field1,field2,..."). For example,
                        if you only want Entry and Entry Name for each output entry, put "accession,id". (default: acc
                        ession,id,protein_name,gene_names,organism_name,length,sequence,go_p,go_c,go,go_f,ft_topo_dom,
                        ft_transmem,cc_subcellular_location,ft_intramem)
  -o, --organism_id
                        (OPTIONAL) Default ID codes for homo sapiens. All seed and main search results will be
                        specific to this organism. Including a specific organism as part of your csv queries may cause
                        an error if it conflicts with the organism specified here. See here for all available IDs:
                        https://www.uniprot.org/taxonomy/?query=* (default: 9606)
```

# Examples
See application note ("Citation") and associated data - https://github.com/Tim-Kirkwood/BulkProt_application_note. 

# Citation
TBC 
