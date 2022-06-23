from lookup_gene import *
from differential_expression import *
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import re

import sys, os
  
  
brain_genes = pd.read_csv('/Users/aaronpresser/Downloads/genes_expressed_in_brain_proteinatlas.tsv', sep='\t')
human_genes = pd.read_csv('/Users/aaronpresser/Downloads/human_genes_proteinatlas.tsv', sep='\t')

# Helper
def synonym_cleaner(string):
    try:
        return string.strip().split(',')
    except AttributeError:
        return ''

synonyms = human_genes['Gene synonym'].apply(synonym_cleaner)


surface_genes = pd.read_csv('/Users/aaronpresser/files_from_quest/Filesfrom_Quest/Peng_RNA_seq/Gene_expression_analysis/surface_genes.csv')
surface_genes.rename(columns={'To': 'Genes'}, inplace=True)
surface_genes.drop(columns=['Unnamed: 0', 'From'], inplace=True)

surfaceome_concat = ' '.join(i for i in surface_genes['Genes'])

# Helper
def str_strip(string):
    return string.strip()

vec_strip = np.vectorize(str_strip)

# Helper
def match_surf(list_of_lists=synonyms):
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
def express_lookup_gene(gene):
    '''lookup_gene() shortcut fit for purposes of genes_expressed()'''
    return data.loc[f"{gene}"][data.loc[f"{gene}"] > 0]

# User
def genes_expressed(cluster, data=data_treat):
    '''Returns a list of all the genes recorded to be expressed in a particular cluster'''
    expressed = []
    for i in data.index:
        try:
            express_lookup_gene(gene=i)[cluster]
            print(f"{i} is expressed in {cluster}")
            expressed.append(i)
        except KeyError:
            print(f"\t{i} isn't expressed")
            continue
    return expressed

gbm_genes = genes_expressed(cluster='6_TREAT')

gbm_surface_genes = [gene for gene in gbm_genes if gene in match_surf(synonyms)]
