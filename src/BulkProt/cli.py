# -*- coding: utf-8 -*-
"""
Created on Mon Oct 14 21:51:00 2024

@author: u03132tk
"""
from argparse import ArgumentParser
import os
import pandas as pd

ALLOWED_FIELDS = ["accession", "id", "gene_names", "gene_primary", "gene_synonym",
                  "gene_oln", "gene_orf", "organism_name", "organism_id", "protein_name",
                  "xref_proteomes", "lineage", "lineage_ids", "virus_hosts",
                  "cc_alternative_products", "ft_var_seq", "error_gmodel_pred", 
                  "fragment", "organelle", "length", "mass", "cc_mass_spectrometry",
                  "ft_variant", "ft_non_cons", "ft_non_std", "ft_non_ter",
                  "cc_polymorphism", "cc_rna_editing", "sequence", "cc_sequence_caution",
                  "ft_conflict", "ft_unsure", "sequence_version", "absorption",
                  "ft_act_site", "cc_activity_regulation", "ft_binding", 
                  "cc_catalytic_activity", "cc_cofactor", "ft_dna_bind", "ec", 
                  "cc_function", "kinetics", "cc_pathway", "ph_dependence", 
                  "redox_potential", "rhea", "ft_site", "temp_dependence", "annotation_score",
                  "cc_caution", "comment_count", "feature_count", "keywordid",
                  "keyword", "cc_miscellaneous", "protein_existence", "reviewed",
                  "tools", "uniparc_id", "cc_interaction", "cc_subunit",
                  "cc_developmental_stage", "cc_induction", "cc_tissue_specificity", 
                  "go_p", "go_c", "go", "go_f", "go_id", "cc_allergen", "cc_biotechnology",
                  "cc_disruption_phenotype", "cc_disease", "ft_mutagen", "cc_pharmaceutical",
                  "cc_toxic_dose", "ft_intramem", "cc_subcellular_location", "ft_topo_dom",
                  "ft_transmem", "ft_chain", "ft_crosslnk", "ft_disulfid",
                  "ft_carbohyd", "ft_init_met", "ft_lipid", "ft_mod_res",
                  "ft_peptide", "cc_ptm", "ft_propep", "ft_signal", "ft_transit", 
                  "structure_3d", "ft_strand", "ft_helix", "ft_turn", "lit_pubmed_id",
                  "date_created", "date_modified", "date_sequence_modified", "version",
                  "ft_coiled", "ft_compbias", "cc_domain", "ft_domain", "ft_motif", 
                  "protein_families", "ft_region", "ft_repeat", "ft_zn_fing"]

def read_args(arg_list: list[str] | None = None):
    parser = ArgumentParser()
    parser.add_argument("-csv", 
                        "--csv_path", 
                        type = str,
                        required=True,
                        metavar='\b',
                        help = '''(REQUIRED) Full filepath to the CSV-format file containing 
                        your queries.''')
                        
    parser.add_argument("-e", 
                        "--excel_compatible", 
                        action = 'store_true',
                        help = '''If set, and you include a sequence field, then 
                        sequence will be truncated to 32000 chars (the maximum
                        number of chars in an excel cell).''')
    parser.add_argument("-q", 
                        "--quick", 
                        action = 'store_true',
                        help = '''If set, will run much faster, but a single query 
                        within your CSV of queries can only return up to 500 results.''')
    parser.add_argument("-sd", 
                        "--seed_only", 
                        action = 'store_true',
                        help = '''If set, will only perform a seed search.''')
    parser.add_argument("-f", 
                        "--fields", 
                        type = str,
                        required=False,
                        default = 'accession,id,protein_name,gene_names,organism_name,'\
                                  'length,sequence,go_p,go_c,go,go_f,ft_topo_dom,'\
                                  'ft_transmem,cc_subcellular_location,ft_intramem',
                        metavar='\b',
                        help = f'''(OPTIONAL - default is accession,id,protein_name,
                        gene_names,organism_name,length,sequence,go_p,go_c,go,go_f,
                        ft_topo_dom,ft_transmem,cc_subcellular_location,ft_intramem) 
                        Result fields you would like to store as columns 
                        in the output CSV (see here: https://www.uniprot.org/help/return_fields).  
                        Input fields should be seperated by a 
                        comma with no spaces (e.g. "field1,field2,...").  
                        
                        POSSIBLE fields are: 
                        {','.join(ALLOWED_FIELDS)}''')
    parser.add_argument("-o", 
                        "--organism_id", 
                        type = int,
                        required=False,
                        default = 9606, #human
                        metavar='\b',
                        help = '''(OPTIONAL - default is 9606, which codes for 
                        homo sapiens) All main search results will be in this organism.  
                        NOTE: This has no effect on the seed search unless you specify 
                        "AND(organism_id:{organism_id})" in the query taken from your 
                        CSV.  See here for all available IDs: 
                        https://www.uniprot.org/taxonomy/?query=*''')
    args = parser.parse_args(arg_list)
    return args, parser

def check_args(args, parser):
    #check csv path - is it a file, does it have one column, does the result dir
    #path derived from the csv filepath already exist
    if not os.path.isfile(args.csv_path):
        parser.exit(f'csv_path {args.csv_path} does not exist')
    table = pd.read_csv(args.csv_path, header = None)
    if len(table.columns) != 1:
        parser.exit(f'{args.csv_path} has {len(table.columns)} columns - it should have 1')
    dir_path = args.csv_path[0:args.csv_path.rindex('.')]
    if os.path.isdir(dir_path):
        parser.exit(f'This folder will be made to hold your results:  "{dir_path}".  '\
                    'It already exists (perhaps from a previous run).  Please move '\
                    '(or change the name of) this folder to prevent it being '\
                    'overwritten.')
    
    #check the organism id - is it an int
    try: 
        int(args.organism_id)
    except ValueError:
        parser.exit(f'organism_id {args.organism_id} is not a recognisable integer')
    
    #check user fields are all recognised by uniprot and have gene/protein name
    #information required to build queries
    user_fields = args.fields.split(',')
    if not all ([i in user_fields for i in ['protein_name', 
                                            'gene_names']
                 ]):
        parser.exit("'protein_name' and 'gene_names' are required fields - "\
                    f"supplied fields are {user_fields}")
                 
    for f in user_fields:
        if f not in ALLOWED_FIELDS:
            parser.exit(f'fields - {f} for recognised')