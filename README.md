# Motivation
UniProt is a great database, but there limited options for bulk searching.  Users can download a lot of entries if they have accessions for specific entries already in hand, and they can perform lots of flexible searches if they can code, but non-coders have no tools for bulk searching.  

BulkProt fills this gap by allowing users to specify a list of search terms in a CSV input file, and presents the output of multiple searches as a single consolidated CSV for easy analysis.  It also allows users to systematically search for entries that are likely to be equivalent to their search term, but are hard to detect via traditional searches due to inconsistent naming conventions.

# Approach

<img src="https://github.com/Tim-Kirkwood/BulkProt/blob/main/Supplementary_Figure_1.png" width="600">

**Supplementary Figure 1: The BulkProt pipeline.**  A CSV-formatted file of UniProt queries is used as input to the BulkProt pipeline.  Each search is queried against the UniProt database (via the API) to generate  a seed table per search term, containing the entries associated with that term.  For each seed table, the protein and gene names are extracted and used to construct a second query, which is again queried against UniProt via the API.  The resulting “main search” table will contain all of the proteins and genes associated with the initial seed table search term, but will also contain irrelevant entries as a result of the unsupervised query construction.  To remove these, entries in the main search table are dropped if they do not have a gene name that is present in the initial seed search table.  Dropped and filtered entries are written to separate tables.  The tables for all queries are concatenated to consolidate results.  The final output of the BulkProt pipeline includes four CSV files:  the seed table, the main table, the filtered table, and the dropped table.

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
