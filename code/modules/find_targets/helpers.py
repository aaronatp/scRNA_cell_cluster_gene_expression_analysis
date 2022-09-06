from lookup_gene import *
from differential_expression import *
import numpy as np
import pandas as pd
import re
  
  
surface_genes = pd.read_csv('./Gene_expression_analysis/surface_genes.csv')
surface_genes.rename(columns={'To': 'Genes'}, inplace=True)
surface_genes.drop(columns=['Unnamed: 0', 'From'], inplace=True)  
  
brain_genes = pd.read_csv('./Downloads/genes_expressed_in_brain_proteinatlas.tsv', sep='\t')
human_genes = pd.read_csv('./Downloads/human_genes_proteinatlas.tsv', sep='\t')

# Helper
def _synonym_cleaner(string):
    try:
        return string.strip().split(',')
    except AttributeError:
        return ''

synonyms = human_genes['Gene synonym'].apply(_synonym_cleaner)
surfaceome_concat = ' '.join(i for i in surface_genes['Genes'])

# Helper
def _str_strip(string):
    return string.strip()

vec_strip = np.vectorize(_str_strip)

# Helper
def _match_surf(list_of_lists=synonyms):
    '''Construct extended_surfaceome - basic gene names for surface genes and also synonyms
        
        - Helpful for identifying GBM genes that are surface genes because you get more options to match
        individual GBM gene with'''
    extended_surfaceome = []
    for lst in list_of_lists:
        if any(match := re.search(gene, surfaceome_concat, flags=re.I) for gene in lst):
            print(f"Match is {match.group(0).strip()}")
            extended_surfaceome.extend(vec_strip(lst))
    
    return extended_surfaceome
  
# Helper
def _express_lookup_gene(gene):
    '''lookup_gene() shortcut fit for purposes of genes_expressed()'''
    return data.loc[f"{gene}"][data.loc[f"{gene}"] > 0]
