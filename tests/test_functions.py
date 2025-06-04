# -*- coding: utf-8 -*-
"""
Created on Fri Jan  3 14:14:16 2025

@author: u03132tk
"""

import pytest 
import pandas as pd
import functions as f
import numpy as np

@pytest.fixture
def uniprot_results_df():
    return pd.DataFrame([['na(m)e (name#) (EC 1.1.1.n11) (EC1.1.1.n11)', 'gene1 gene2 gene3'],
                         ['name (name) (EC 1.1.1.n11) (EC1.1.1.n11)', 'gene2 gene3 gene6'],
                         [np.nan, np.nan],
                         ['name (cofactor) (name) (EC 1.1.1.n11) (EC1.1.1.n11)', 'gene1 gene2 gene3'],
                         ['name(cofactor) (EC 1.1.1.n11) (EC1.1.1.n11)', ''],
                         [np.nan, np.nan],
                         ['name (EC 1.1.1.n11) (EC1.1.1.n11)', 'gene2 gene3 gene6'],
                         ['3-alpha-(17-beta)-hydroxysteroid dehydrogenase (NAD(+)) (EC 1.1.1.n11) (EC1.1.1.n11)', np.nan],
                         ['another protein (3-alpha-(17-beta)-hydroxysteroid dehydrogenase (NAD(+))) (EC 1.1.1.n11) (EC1.1.1.n11)', np.nan],
                         ['another protein (another protein) (3-alpha-(17-beta)-hydroxysteroid dehydrogenase (NAD(+))) (EC 1.1.1.n11) (EC1.1.1.n11)', np.nan]
                         ['name [Includes: x y z] [not_pattern a] (EC 1.1.1.n11) (EC1.1.1.n11)', ''],
                         [np.nan, np.nan],
                         ['name [not_pattern a] [Cleaved into: x y z] (EC 1.1.1.n11) (EC1.1.1.n11)', np.nan],
                         ['name [not_pattern a] (EC 1.1.1.n11) (EC1.1.1.n11)', '']
                         ],
                        columns = ['Protein names', 'Gene Names']
                        )

#protein name parsing is very fragile and relies a lot on uniprot format consistency 
#i am not going to optimise or fully test this as there are obvious edge tests 
#that will fail without the aforsaid optimisation, and less than ideal outputs 
#(e.g. trailing spaces for square bracket removal).  It is much easier and more
#robust to get the names from a different file format (JSON) which is intended 
#for name extraction and does not need explit string parsing at the user end.  
#I will do this in the next version
@pytest.mark.parametrize('protein_names, expected_split_protein_names', 
                         [('na(m)e (name) (name)', ['na(m)e', 'name', 'name'])
                          ('name (name) (name)', ['name', 'name', 'name'])
                          ('name (cofactor) (name)', ['name', 'cofactor', 'name'])
                          ('name(cofactor)', ['name(cofactor'])
                          ('name', ['name']),
                          ('3-alpha-(17-beta)-hydroxysteroid dehydrogenase (NAD(+))', 
                           ['3-alpha-(17-beta)-hydroxysteroid dehydrogenase', 
                            'NAD(+)'],
                           )
                          ('another protein (another protein) (3-alpha-(17-beta)-hydroxysteroid dehydrogenase (NAD(+)))', 
                           ['another protein', 
                            'another protein', 
                            '3-alpha-(17-beta)-hydroxysteroid dehydrogenase (NAD(+))']
                           ),
                          ('another protein (3-alpha-(17-beta)-hydroxysteroid dehydrogenase (NAD(+)))', 
                           ['another protein', 
                            '3-alpha-(17-beta)-hydroxysteroid dehydrogenase', 
                            'NAD(+)']
                           )
                           ]
                          )
def test_split_string(protein_names, expected_split_protein_names):
    assert f.split_string(protein_names) == expected_split_protein_names

@pytest.mark.parametrize('name, clean_name', 
                         [('name [pattern: x y z] [not_pattern a]', 
                           'name  [not_pattern a]'),
                          ('name [not_pattern a] [pattern: x y z]', 
                           'name [not_pattern a] ')
                          ('name [not_pattern a]'), 'name [not_pattern a]'
                          ('name [pattern: x y z]', 'name ')
                          ]
                         )
def test_remove_square_bracket_info (name, clean_name):
    assert f.remove_square_bracket_info(name) == clean_name


@pytest.mark.parametrize('name, result', 
                         [('EC 1.1.1.n11', True),
                          ('EC1.1.1.n11', False),
                          ('Aldoketoreductase 7', False)
                          ]
                         )

def test_is_ec(name, result):
    assert f.is_ec(name) == result

class TestCleanProteins:
    
    @pytest.fixture
    def exclude_ec_keep_cofactor(self):
        return {'na(m)e', 'name', 'name#', 'cofactor', 
                'name(cofactor', '3-alpha-(17-beta)-hydroxysteroid dehydrogenase', 
                'NAD(+)', 'another protein', 'name  [not_pattern a]',
                '3-alpha-(17-beta)-hydroxysteroid dehydrogenase (NAD(+))', 
                'name [not_pattern a] ', 'name [not_pattern a]'}
    
    @pytest.fixture
    def keep_ec_exclude_cofactor(self):
        return {'na(m)e', 'name', 'name#', 'EC 1.1.1.n11', 'EC1.1.1.n11', 
                'name(cofactor', '3-alpha-(17-beta)-hydroxysteroid dehydrogenase', 
                'NAD(+)', 'another protein', 'name  [not_pattern a]',
                '3-alpha-(17-beta)-hydroxysteroid dehydrogenase (NAD(+))', 
                'name [not_pattern a] ', 'name [not_pattern a]'}
    
    @pytest.fixture
    def exclude_ec_exclude_cofactor(self):
        return {'na(m)e', 'name', 'name#', 
                'name(cofactor', '3-alpha-(17-beta)-hydroxysteroid dehydrogenase', 
                'NAD(+)', 'another protein', 'name  [not_pattern a]',
                '3-alpha-(17-beta)-hydroxysteroid dehydrogenase (NAD(+))', 
                'name [not_pattern a] ', 'name [not_pattern a]'}
    
    @pytest.fixture
    def keep_ec_keep_cofactor(self):
        return {'na(m)e', 'name', 'name#', 'EC 1.1.1.n11', 'EC1.1.1.n11', 'cofactor', 
                'name(cofactor', '3-alpha-(17-beta)-hydroxysteroid dehydrogenase', 
                'NAD(+)', 'another protein', 'name  [not_pattern a]',
                '3-alpha-(17-beta)-hydroxysteroid dehydrogenase (NAD(+))', 
                'name [not_pattern a] ', 'name [not_pattern a]'}
    
    @pytest.mark.parametrize('exclude_ec, exclude_specific, expected_results_fixture', 
                             [(True, [], 'exclude_ec_keep_cofactor'),
                              (False, ['cofactor'], 'keep_ec_exclude_cofactor'),
                              (True, ['cofactor'], 'exclude_ec_exclude_cofactor'),
                              (False, [], 'keep_ec_keep_cofactor')
                              ]
                             )
    def test_normal(self, uniprot_results_df, request, exclude_ec, exclude_specific, expected_results_fixture):
        expected_output = request.getfixturevalue(expected_results_fixture)
        assert expected_output == f.clean_proteins(uniprot_results_df,
                                                   exclude_ec, 
                                                   exclude_specific)
        
class TestCleanGenes:
    
    @pytest.fixture
    def results_exclude_gene_6(self):
        return {'gene1', 'gene2', 'gene3'}
    
    @pytest.fixture
    def results(self):
        return {'gene1', 'gene2', 'gene3', 'gene6'}
    
    @pytest.mark.parametrize('exclude_specific, results_fixture',
                             [([], 'results'),
                              (['gene6'], 'results_exclude_gene_6'),
                              (['not_in_name'], 'results')
                              ]
                             )
    def test_normal(self, request, uniprot_results_df, exclude_specific, results_fixture):
        results = f.clean_genes(uniprot_results_df,
                                exclude_specific)
        expected_results = request.getfixturevalue(results_fixture)
        assert results == expected_results

class TestDfToQuery:
    
    @pytest.fixture
    def setup_df_to_query(self, monkeypatch):
        #clean genes and clean proteins return sets so you cannot be sure of the 
        #order - have to monkeypatch
        def patch_clean_genes(df, exclude_specific):
            #return a list so you can be sure of order
            return ['gene1', 'gene2#', 'gene3']
        def patch_clean_proteins(df, exclude_specific):
            return ['protein1', 'protein2#', 'protein3']
        monkeypatch.setattr('Bulkprot.functions.clean_genes', 
                            patch_clean_genes)
        monkeypatch.setattr('Bulkprot.functions.clean_proteins', 
                            patch_clean_proteins)
        

    
    @pytest.fixture
    def results(self):
        return '(gene:"gene1")OR(gene:"gene2%23")OR(gene:"gene3")'\
               'OR(protein_name:"protein1")OR(protein_name:"protein2%23")'\
               'OR(protein_name:"protein3")'
    @pytest.fixture
    def results_exclude_gene1_protein1(self):
        return '(gene:"gene2%23")OR(gene:"gene3")OR(protein_name:"protein2%23")'\
               'OR(protein_name:"protein3")'
               
    @pytest.mark.parametrize('exclude_specific, results_fixture',
                             [([], 'results'),
                              (['gene1', 'protein1'], 'results_exclude_gene1_protein1'),
                              (['not_in_name'], 'results')
                              ]
                             )
    def test_normal(self, request, setup_df_to_query, uniprot_results_df,
                    exclude_specific, results_fixture):
        expected_output = request.getfixturevalue(results_fixture)
        assert expected_output == f.df_to_query(uniprot_results_df, 
                                                exclude_specific)
        
class TestQueriesToTable:
    #monkeypatch request.get - return headings, raise connection error, or raise status code
    def setup_ok():
        pass
    def setup_connection_error()