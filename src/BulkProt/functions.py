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
                               #include space as protein names can have brackets eg Delta(4)-3-oxosteroid 5-beta-reductase
            #print (clean_proteins)
    #dont include ec numbers - these refer to reactions not proteins, so potential for false positives
    if exclude_ec:
        clean_proteins = [p for p in clean_proteins if not is_ec(p)]
    
    #print (clean_proteins)
    return set(clean_proteins)

def clean_genes(df, exclude_specific = []):
    genes = []
    for hit_genes in df['Gene Names']:
        if not hit_genes is np.nan:
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
    rest_url = base + f'query=(({query})AND(organism_id:{organism_id}))'
    response = requests.get(rest_url)
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
        if row is np.nan:
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
    #print(df)
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
# def check_df(df, size):
#     return {'no_hits' : len(df.index) == 0,
#             'too_many_hits' : len(df.index) == size}

# def update_bad_df(df, status):
#     #not(0 and 0) == 1, not 0 and 0 == 0
#     assert not (status['no_hits'] and status['too_many_hits'])
#     if status['no_hits']:
#         message = 'No UniprotKB hits'
#     elif status['too_many_hits']:
#         message = 'Too many UniprotKB hits - make this query more stringent'
#     return pd.DataFrame([[message]*len(df.columns)],
#                         columns = df.columns)

def BulkProt(filepath : str, fields, organism_id, seed_only, excel_compatible):
    
    #initialise super lists and constants
    error_message = 'No UniprotKB hits'
    seed_all = []
    main_all = []
    filtered_all = []
    dropped_all = []  
    url_base = f'https://rest.uniprot.org/uniprotkb/stream?fields={fields}&format=tsv&'    
    #read in data
    table = pd.read_csv(filepath, header = None)
    num_rows_pre = {'dropped' : 0,
                'filtered' : 0,
                'main' : 0,
                'seed' : 0}
    num_rows_post = {'dropped' : 0,
                'filtered' : 0,
                'main' : 0,
                'seed' : 0}
    #process data
    start = time.time()
    for row_index, input_csv_row in table.iterrows():
        if row_index % 500 == 0:
            print (f'Done {row_index} queries ({int(time.time() - start)} seconds)')
        #set these as None in the loop - if they are not reset then they carry 
        #over the loops (e.g. df_main)
        df_main = None
        df_filtered = None
        df_dropped = None
        main_query = 'Not performed' #needs to be initialised for general_col_map 
        print (f'USER QUERY: {input_csv_row[0]}')
        
        #send user query to uniprot and convert response to SEED dataframe 
        df_seed = queries_to_table(url_base, 
                                   input_csv_row[0], 
                                   organism_id)
        #print(df_seed.index)
        #check you have hits - if not, (i) update dataframe to reflect issues 
        #and (ii) do not make MAIN and FILTERED dataframes
        if len(df_seed.index) == 0:#seed_status['no_hits'] or seed_status['too_many_hits']:
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
                #print(df_main.index)
                #check if the MAIN dataframe is acceptable - if not, (i) update 
                #dataframe in place to reflect issues and (ii) do not make FILTERED 
                #dataframe.  Otherwise, filter MAIN dataframe by removing rows that
                #do not have a SEED gene or protein.  Surviving rows are written 
                #to FILTERED dataframe, dropped rows are written to DROPPED
                #dataframe
                if len(df_main.index) == 0:
                    df_main = pd.DataFrame([[error_message]*len(df_main.columns)],
                                             columns = df_main.columns)
                    print ('MAIN QUERY: Error - no hits \n'\
                           'FILTERED: Not performed')
                else:
                    drop_indexes = get_drop_indexes(df_seed, df_main)
                    df_filtered = df_main.drop(drop_indexes)
                    #print (df_filtered.index)
                    df_dropped = df_main.loc[drop_indexes]
                    print (f'FILTERED: Dropped {len(drop_indexes)} out of {len(df_main.index)} entries.  New df has {len(df_filtered.index)} entries.')
        
        #update dataframes with SEED and MAIN query details.  Add one column each,
        #and write query to first row in the column.  Leave other columns blank 
        #to aid with manual curation of the results.  Add formatted dfs to their
        #respective master lists
        seed_col_map = {'Seed query' : input_csv_row[0]
                        }
        general_col_map = {'Main query' : main_query,
                           'Seed query' : input_csv_row[0]
                           }
        for name, df, col_map, parental_list in [('dropped',df_dropped, general_col_map, dropped_all), 
                                                 ('filtered',df_filtered, general_col_map, filtered_all),
                                                 ('main',df_main, general_col_map, main_all),
                                                 ('seed', df_seed, seed_col_map, seed_all)]:
            if df is not None:
                #print(name)
                #pre_df = len(df.index)
                #print (f'PRE {name}: {len(df.columns)} cols, {len(df.index)} rows')
                #pre_df = copy.deepcopy(df)
                #pre_df.to_csv('D:/BulkProt/src/BulkProt/FAIL_pre_df.csv')
                num_rows_pre[name] += len(df.index)
                #this is adding a row in certain situations - not sure why
                update_cols_inplace(df, col_map)
                #print (f'POST {name}: {len(df.columns)} cols, {len(df.index)} rows')
                #post = len(df.index)
                #post_df = copy.deepcopy(df)
                #post_df.to_csv('D:/BulkProt/src/BulkProt/FAIL_post_df.csv')
                num_rows_post[name] += len(df.index)
                #if name != 'dropped':
                    #try:
                    #    assert len(pre_df.index) == len(post_df.index)
                    #except AssertionError as e:
                    #    print (f'PRE:\n\n{pre_df}\n\nPOST:\n\n{post_df}\n\nmain_query \n\n{main_query} \n\ncol_map\n\n{col_map}')
                    #    
                    #    raise e
                parental_list += [df]
        print ()
    print (num_rows_pre) 
    print (num_rows_post)
    #Concatenate all dataframes and write to CSV format in a input-specific 
    #results dir 
    new_dir = build_dir(filepath)
    for name, df_list, path in [('seed', seed_all, f'{new_dir}/seed.csv'),
                                ('main', main_all, f'{new_dir}/main.csv'),
                                ('filtered', filtered_all, f'{new_dir}/filtered.csv'),
                                ('dropped', dropped_all, f'{new_dir}/dropped.csv')]:
        if len(df_list) > 0:
            print (f'Writing results to {path}')
            #TODO 
            print (name.upper())
            print (f'number of rows before concat {sum([len(d.index) for d in df_list])}')
            df = pd.concat(df_list)
            print (f'number of rows after concat {len(df.index)}')
            if excel_compatible:
                #https://stackoverflow.com/a/49451329/11357695
                if 'Sequence' in df.columns:
                    #https://stackoverflow.com/a/79095162/11357695
                    df['Sequence'] = df['Sequence'].str.slice(0,32000)
            print (f'number of rows after excel formatting {len(df.index)}')
            df.to_csv(path)
            print (f'{name} {len(df.index)}')
    return {'seed_all' : seed_all, 
            'main_all' : main_all, 
            'filtered_all' : filtered_all, 
            'dropped_all' : dropped_all}
            
    