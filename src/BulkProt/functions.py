# -*- coding: utf-8 -*-
"""
Spyder Editor

This is a temporary script file.
"""
import pandas as pd
import requests
import io
import numpy as np
import os
import time
import copy
from pandas.api.types import is_string_dtype
'''
future plans:
    at some point would be good to add format options etc rather than hard coding
    add in pagination to deal with > 500 hits (https://www.uniprot.org/help/api_queries)
    add in filter step to get rid of eg '{query}-associated' hits
    #TODO filter by {query}_associated, {query}-associated, {query} associated,
    #{query}_interacting, {query}-interacting, {query} interacting,
    #anything that has at least one gene but does not have query gene (will 
    #be sub/superset)
    #TODO add survivers to master filtered df - dont forget to include column 
    #for user search term and expanded search term
    #TODO add culled to master culled df - dont forget to include column 
    #for user search term and expanded search term
    #TODO def sanitise input - check all things are recognised. I think this is done 
    #already via the API fail codes but could be worth checking for white space etc 
'''

#TODO add pattern options etc to the cli so user has control

def split_string(s):
    #it is assumed that s is formatted as 'n1 (n2) (n3)...'
    #regex doesnt work well as some of the names can have matching patterns
    #') (' is a fairly stringent pattern, but ' (' can have hits in the protein 
    #names - for example, '3-alpha-(17-beta)-hydroxysteroid dehydrogenase (NAD(+))'.
    #This function looks at ') (' as a split pattern, and then looks for ' (' in
    #the first name only, as this is the only place where you get ' (' instead of ') ('.
    #NOTE this will fail when the first name matches this pattern.  For example, 
    #if '3-alpha-(17-beta)-hydroxysteroid dehydrogenase (NAD(+))' is the first name, 
    #then (NAD(+)) looks like n2 in the expected format.  I can't see a way to 
    #distinguish this as it is a fundamentally inconsistent pattern.  However, 
    #if '3-alpha-(17-beta)-hydroxysteroid dehydrogenase (NAD(+))', or another 
    #name with the ' (' pattern, is not the first name, then no issues.
    
    #TODO add a check that if you have brackets then the next set follow expected pattern.
    #e.g. this one doesnt for Q59HF2_HUMAN:
    #'D site of albumin promoter (Albumin D-box) binding protein variant'
    #'HNF1 homeobox A (Transcription factor 1, hepatic LF-B1, hepatic nuclear factor (HNF1), albumin proximal factor, isoform CRA_b)'
    #-> ['HNF1 homeobox A',
    #    'Transcription factor 1, hepatic LF-B1, hepatic nuclear factor',
    #    'HNF1), albumin proximal factor, isoform CRA_b']
    #i think you need a bracket counting strategy as in find_delete_indicies() - perhaps look at making it possible to define delims
    
    if s[-1] == ')': #do not include this as it is not a split pattern
        s = s[0:-1]
    partially_split_string = s.split(') (')
    split_first_obj = partially_split_string[0].split(' (')
    full_split_string = split_first_obj + partially_split_string[1:]
    return full_split_string        

def remove_square_bracket_info(s, patterns = ['[Includes:', '[Cleaved into:']):
    def find_delete_indicies(s, pattern):
        open_bracket = False
        internal_brackets = 0
        delete_indicies = []
        for i, c in enumerate(s):
            if s[i:i+len(pattern)] == pattern:
                open_bracket = True
            if open_bracket:
                if c == '[':
                    internal_brackets += 1 
                elif c == ']':
                    internal_brackets -= 1
                    if internal_brackets == 0:
                        open_bracket = False
                        delete_indicies +=[i]#add final close bracket
            #CODE SMELL - delete 'if open bracket' i think, then can also delete 
            #the extra 'delete indicies' add as it will be executred in the loop
            if open_bracket:
                delete_indicies += [i]
        assert not open_bracket
        assert internal_brackets == 0
        return delete_indicies
                
    delete_indicies = [] 
    for pattern in patterns:
        delete_indicies += find_delete_indicies(s, pattern)
    return ''.join([c for i, c in enumerate(s) if i not in delete_indicies]).strip()
        
        

def is_ec(s):
    if s[0:3] == 'EC ':
        numbers = s.split(' ', 1)[1].split('.')
        ok_format = len(numbers) == 4
        #EC numbers can have dashes and letters eg Aflatoxin B1 aldehyde reductase member 2 (EC 1.1.1.n11) (AFB1 aldehyde reductase 1) (AFB1-AR 1) (Aldoketoreductase 7) (Succinic semialdehyde reductase) (SSA reductase)
        #try:
        #    list(map(int, numbers))
        #except ValueError:
        #    ok_format = False 
    else:
        ok_format = False
    return ok_format
        
def clean_proteins(df, exclude_ec = True, exclude_specific = []):#['NAD(+)']):
    clean_proteins = []
    for hit_proteins in df['Protein names']:
        if not hit_proteins is np.nan:
            simple_name = remove_square_bracket_info(hit_proteins)
            split_proteins = split_string(simple_name)
            clean_proteins += [i for i in split_proteins if i not in exclude_specific]
    #dont include ec numbers - these refer to reactions not proteins, so potential for false positives
    if exclude_ec:
        clean_proteins = [p for p in clean_proteins if not is_ec(p)]
    
    return set(clean_proteins)

def clean_genes(df, exclude_specific = []):
    genes = []
    for hit_genes in df['Gene Names']:
        if not pd.isna(hit_genes):
            genes += [i for i in hit_genes.split(' ') if i not in exclude_specific] 
    return set(genes) 

def df_to_query(df):#, gene_only = True):
    exclude_specific = ['putative', 'NAD(+)', 'variant','protein']
    cl_genes = clean_genes(df, exclude_specific)
    queries = [f'(gene:"{gene.strip()}")OR' for gene in cl_genes]
    #if not gene_only:
    cl_proteins = clean_proteins(df, exclude_specific)
        #'#' was breaking my url 
    queries += [f'(protein_name:"{protein.strip().replace("#", "%23")}")OR' for protein in cl_proteins]
    query = ''.join(queries)[0:-2] #remove final OR
    return query

def queries_to_table(base, query, organism_id):
    rest_url = base + f'query=({query})AND(organism_id:{organism_id})'
    print (f'REST URL:  {rest_url}')
    for attempt in range(5):
        try:
            response = requests.get(rest_url)
        except requests.exceptions.ConnectionError:
            continue
        if response.status_code == 200:
            return pd.read_csv(io.StringIO(response.text), 
                               sep = '\t')
        else:
            raise ValueError(f'The uniprot API returned a status code of {response.status_code}.  '\
                             'This was not 200 as expected, which may reflect an issue '\
                             f'with your query:  {query}.\n\nSee here for more '\
                             'information: https://www.uniprot.org/help/rest-api-headers.  '\
                             f'Full url: {rest_url}')


def build_dir(filepath):
    assert os.path.isfile(filepath)
    dir_path = filepath[0:filepath.rindex('.')]
    os.makedirs(dir_path, exist_ok=False)
    return dir_path

def get_drop_indexes(df_seed, df_main):
    seed_genes = clean_genes(df_seed)
    drop_indexes = []
    for index, row in enumerate(df_main['Gene Names']):
        if pd.isna(row):
            #keep for manual checking
            continue
        genes = set(row.split(' '))
        num_shared_genes = len(genes.intersection(seed_genes))
        if num_shared_genes == 0:
            drop_indexes += [index]
    return drop_indexes

def update_cols_inplace(df, col_map):
    num_row = len(df.index)
    num_col = len(df.columns)
    for col, val in col_map.items():
        df.insert(0, col,'')
        #you know col index is 0 as you just added it
        try:
            df.iloc[0, 0] = val
        except IndexError:
            assert df.empty
    #print(df)
    assert num_row == len(df.index)
    assert num_col == len(df.columns) - len(col_map)

def build_url_base(quick, fields):
    if quick:
        api_func = 'search'
    else:
        api_func = 'stream'
    url_base = f'https://rest.uniprot.org/uniprotkb/{api_func}?fields={fields}&format=tsv&'  
    return url_base

def check_input(table):
    one_col = len(table.columns) == 1 
    if not one_col:
        msg = "input table has multiple columns (you should check for whitespace chars)"
        raise ValueError(msg)    
    all_values = list(table.iloc[:,0])
    unique_values = len(all_values) == len(set(all_values))
    if not unique_values:
        redundant_values = [i for i in set(all_values) if all_values.count(i) > 1]
        msg = f"input table has non-unique queries - {redundant_values})"
        raise ValueError(msg)

def process_query(input_csv_row, url_base, organism_id, seed_only):
    #if there are no hits following seed or main search
    error_message = 'No UniprotKB hits'
    
    #Initialise variables - sometimes they will be updated, othertimes they 
    #will not (e.g. when their are no hits to a given seed query, seed will 
    #be 'No UniprotKB hits', main/filtered will be None, main_query will be 
    #'Not performed').
    df_main = None
    df_filtered = None
    df_dropped = None
    seed_query = input_csv_row[0]
    #less redundant code if you just define at the start and update
    main_query = 'Not performed'  
    
    #send user query to uniprot and convert response to SEED dataframe 
    df_seed = queries_to_table(url_base, 
                               input_csv_row[0], 
                               organism_id)
    
    #check you have hits - if not, (i) update dataframe to reflect issues 
    #and (ii) do not make MAIN and FILTERED dataframes
    if len(df_seed.index) == 0:
        df_seed = pd.DataFrame([[error_message]*len(df_seed.columns)],
                                 columns = df_seed.columns)
        if not seed_only:
            print ('SEED QUERY: Error - no hits\n'\
                   'MAIN QUERY: Not performed\n'\
                   'FILTERED: Not performed')
                
    #If the SEED dataframe is ok 
    else:
        #and the user wants to perform a MAIN search
        if not seed_only:
            
            #convert the SEED data to a new query
            main_query = df_to_query(df_seed)
            print (f'MAIN QUERY: {main_query}')
            
            #send query to uniprot and convert response to MAIN dataframe 
            df_main = queries_to_table(url_base, main_query, organism_id)
    
            #check if the MAIN dataframe is acceptable - if not, (i) update 
            #dataframe in place to reflect issues and (ii) do not make FILTERED 
            #dataframe.  
            if len(df_main.index) == 0:
                df_main = pd.DataFrame([[error_message]*len(df_main.columns)],
                                         columns = df_main.columns)
                print ('MAIN QUERY: Error - no hits \n'\
                       'FILTERED: Not performed')
                    
            #Otherwise, filter MAIN dataframe by removing rows that
            #do not have a SEED gene or protein.  Surviving rows are written 
            #to FILTERED dataframe, dropped rows are written to DROPPED
            #dataframe
            else:
                drop_indexes = get_drop_indexes(df_seed, df_main)
                df_filtered = df_main.drop(drop_indexes)
                df_dropped = df_main.loc[drop_indexes]
                print (f'FILTERED: Dropped {len(drop_indexes)} out of '\
                       f'{len(df_main.index)} entries.  New df has '\
                       f'{len(df_filtered.index)} entries.')
    return {'df_seed' : df_seed,
            'df_main' : df_main,
            'df_filtered' : df_filtered,
            'df_dropped' : df_dropped,
            'seed_query' : seed_query,
            'main_query' : main_query}


def BulkProt(filepath : str, fields, organism_id, seed_only, excel_compatible, 
             quick):
    
    #initialise lists and constants
    seed_all = []
    main_all = []
    filtered_all = []
    dropped_all = []  
    url_base = build_url_base(quick, fields)  
    
    #read in data
    new_dir = build_dir(filepath)
    table = pd.read_csv(filepath, header = None)
    check_input(table)
            
    #process data
    for row_index, input_csv_row in table.iterrows():
        
        
        
        print (f'\nUSER QUERY: {input_csv_row[0]}')
        
        #get seed, main, filtered and dropped tables for this query
        query_results = process_query(input_csv_row, url_base, organism_id, seed_only)
        
        #assign data to vars
        seed_query = query_results['seed_query']
        main_query = query_results['main_query']
        df_main = query_results['df_main']
        df_filtered = query_results['df_filtered']
        df_dropped = query_results['df_dropped']
        df_seed = query_results['df_seed']
        
        
        #update dataframes with SEED and MAIN query details.  Add one column each,
        #and write query to first row in the column.  Leave other columns blank 
        #to aid with manual curation of the results.  Add formatted dfs to their
        #respective master lists
        seed_col_map = {'Seed query' : seed_query}
        general_col_map = {'Main query' : main_query,
                           'Seed query' : seed_query}
        for df, col_map, parental_list in [(df_dropped, general_col_map, dropped_all), 
                                           (df_filtered, general_col_map, filtered_all),
                                           (df_main, general_col_map, main_all),
                                           (df_seed, seed_col_map, seed_all)]:
            if df is not None:
                update_cols_inplace(df, col_map)
                parental_list += [df]

    #Concatenate all dataframes and write to CSV format in a input-specific 
    #results dir 
    
    for name, df_list, path in [('seed', seed_all, f'{new_dir}/seed.csv'),
                                ('main', main_all, f'{new_dir}/main.csv'),
                                ('filtered', filtered_all, f'{new_dir}/filtered.csv'),
                                ('dropped', dropped_all, f'{new_dir}/dropped.csv')]:
        if len(df_list) > 0:
            df = pd.concat(df_list)
            print (f'{name.upper()} - Writing {len(df.index)} rows to {path}')
            if excel_compatible:
                print ('Trimming all values with <32000 chars')
                #https://stackoverflow.com/a/45270483/11357695 - note applymap replaced by map
                df = df.map(lambda x: x[0:32000] if isinstance(x, str) else x)
            df.to_csv(path)

    return {'seed_all' : seed_all, 
            'main_all' : main_all, 
            'filtered_all' : filtered_all, 
            'dropped_all' : dropped_all}
            
    