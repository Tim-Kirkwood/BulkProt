# Motivation
UniProt is a great database, but there are limited options for bulk searching.  Users can download a lot of entries, if they have the associated accessions already in hand, and they can perform lots of flexible searches if they can code.  However, non-coders have no tools for bulk searching.  

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
*(These examples are copied from the "Usage" section of the BulkProt Application Note ("Citation") and the associated Supplementary Data is provided at https://github.com/Tim-Kirkwood/BulkProt_application_note.)*

Two example searches were performed based on previous work in the lab. In the first example, we aimed to annotate our internal proteomics database with ADME categories taken from supplementary Tables S3 and S4 in (Achour et al., 2021). To do this, we had to search for all Uniprot entries linked to genes listed in the supplementary tables. These Uniprot entries could then be cross-referenced with our internal database, allowing us to annotate the relevant database entries with the associated ADME category. This was achieved using the command:

**bulkprot -csv “path/to/query.csv” –seed only**

The ”seed only” parameter indicates that a main search was not needed (as this was a simple search task, rather than an attempt to link entries with different names). The placeholder ”path/to/query.csv” refers to the file path for either Supplementary Data 1 or Supplementary Data 2, which contain the BulkProt queries derived from Supplementary Table 3 (173 searches) and Supplementary Table 4 (362 searches) in (Achour et al., 2021), respectively. The respective seed table outputs generated by BulkProt are given in Supplementary Data 3 and Supplementary Data BulkProt 3 4. Performing 535 separate UniProt searches via the web portal would take approximately 4 to 5 hours, assuming roughly 30 seconds per search (this includes time taken to perform a single search, export that search’s results to excel, and to collate all search excel files into a single, interpretable excel file). Whilst this will be dependent on internet speed, generating Supplementary Data 3 and Supplementary Data 4 from Supplementary Data 1 and Supplementary Data 2 took
under 9 minutes in total.  

Next, we will give an example of target searching using two separate search terms, (protein name:”High mobility group protein B1”) and (protein name:”Heterogeneous nuclear ribonucleoprotein Q”). This was performed using the same command as above, but without the “-seed only” parameter – see Supplementary Data 5 for the query file. A comparison of the seed table (Supplementary Data 6) with the filtered table (Supplementary Data 7) shows that the filtered BulkProt results for Heterogeneous nuclear ribonucleoprotein Q identified many cognate entries that were deposited with a different name (e.g. Synaptotagmin binding cytoplasmic RNA interacting protein), which would not have been retrieved using the original search term alone. No single protein name was present across all entries, making it impossible to identify all entries by searching UniProt for a single protein name. Similarly, some entries lack “Gene Names”, and so no gene name searches would identify all entries in the table. In contrast, BulkProt was able to find all entries associated with this protein - the full output is available in Supplementary Data 7. This highlights the utility of BulkProt’s automated expansion step for capturing entries annotated under diverse nomenclature – see Figure 1 for a coloured example that highlights different entry names. Missing entries are also observed for High mobility group protein B1 (Supplementary Data 7, accessions A0A286R9D9, A0A286R9F1 and Q9NYD7). However, this search term also highlights the need for manual curation of the BulkProt results. Also present in the filtered table are accessions B2RPK0 (derived from the seed search) and B2RDE8 (derived from the main search), which derive from pseudogenes that have similar names to the query, but which are unlikely to be considered equivalent to the protein in question.


# Citation
TBC 
