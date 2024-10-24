# -*- coding: utf-8 -*-
"""
Created on Wed Oct 16 17:44:30 2024

@author: u03132tk
"""


#BUG 1 - seems to be an excel issue.  see discusiion about maximum cell length here https://stackoverflow.com/a/79095162/11357695
import requests 
import pandas as pd
import io 

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

size = 500
fields = 'accession,id,protein_name,gene_names,organism_name,'\
          'length,sequence,go_p,go_c,go,go_f,ft_topo_dom,'\
          'ft_transmem,cc_subcellular_location,ft_intramem'
url_base = f'https://rest.uniprot.org/uniprotkb/search?'\
           f'fields={fields}&format=tsv&'
query = '(gene:"GPX8")OR(gene:"UNQ847/PRO1785")OR(protein_name:"Probable glutathione peroxidase 8")OR(protein_name:"GSHPx-8")OR(protein_name:"GPx-8")OR(protein_name:"Glutathione peroxidase 8")OR(protein_name:"Glutathione peroxidase")OR(protein_name:"putative")'
organism_id = 9606
df = queries_to_table(url_base, query, organism_id)
#-> df looks fine - one row and 15 columns
dfout = pd.concat([df])#
###this fixes it https://stackoverflow.com/a/49451329/11357695
dfout['Sequence'] = df['Sequence'].str.slice(0,32000)
##############
dfout.to_csv('test2_error.csv')
#-> opening in excel this is broken - it splits df['Sequence'] into two rows at 
#the junction between 'RLLANAECQEGQSVCFEIRVSGIPPPTLKWEKDG' and 
#'PLSLGPNIEIIHEGLDYYALHIRDTLPEDTGYY'. In df['Sequence'], this sequence is joined 
#by a 'q':
#tdstlrpmfkRLLANAECQEGQSVCFEIRVSGIPPPTLKWEKDGqPLSLGPNIEIIHEGLDYYALHIRDTLPEDTGYYrvtatntags

dfin = pd.read_csv('test2_error.csv').fillna(0).drop(['Unnamed: 0'], axis = 1)
res = dfin == dfout.fillna(0) #all true

#CONCLUSION - excel is mangling it 