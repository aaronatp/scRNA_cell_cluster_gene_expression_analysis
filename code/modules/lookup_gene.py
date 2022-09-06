# data is the raw gene expression counts data
# cells_to_clusters tells you what cluster each cell belongs to

import numpy as np
import pandas as pd
import matplotlib.pyplot as plt

import sys, os

from collections import Counter

# Helper
class HiddenPrints:
    '''Hides print statements of code called within it'''
    def __enter__(self):
        self._original_stdout = sys.stdout
        sys.stdout = open(os.devnull, 'w')

    def __exit__(self, exc_type, exc_val, exc_tb):
        sys.stdout.close()
        sys.stdout = self._original_stdout   

# User    
def sorted_counter_dict(gene, data=data, proportion=True):
    """Returns a dict of clusters and proportion of cells that express 'gene'"""
    gene_names = data.loc[f'{gene}'][data.loc[f'{gene}'] > 0].index
    num_cells_in_each_cluster = cells_to_clusters['Cell_total']
    
    gene_dict = {}
    for i, j in zip(gene_names, num_cells_in_each_cluster):
        try:
            gene_dict[i] += 1
        except KeyError:
            gene_dict[i] = 1
    
    if proportion == True:
        gene_dict = {k: i/data_dict[k] for k,i in gene_dict.items()}
    else:
        gene_dict = {k: i for k,i in gene_dict.items()}
    return {k: v for k, v in sorted(gene_dict.items(), key=lambda item: item[1], reverse=True)}

# Helper
def _sort_dict_by_freq(store):
    '''Sorts a dictionary's items by the frequency of the values in it'''
    freq = {} 
    for value in store.values(): 
        if freq.get(value) == None: 
            freq[value] = 1 
        else: 
            freq[value] += 1 
            pass 
    pos = sorted(freq.items(), key=lambda item: item[1], reverse=True) 
    pos = list(map(lambda x: x[0], pos)) 
    store = dict(sorted(store.items(), key=lambda item: pos.index(item[1]))) 
    return store


def _gene_expression_helper(gene, sure, df):
    options = df.filter(like=f'{gene}', axis=0).index
    if len(options) == 0:
        return print("Try again")
    elif len(options) > 1 and sure==False:
        list_of_options = '\n'.join(i for i in options)
        print(f"There are multiple genes that match your search:")
        return print(list_of_options)
    elif len(options) == 1 or sure==True:
        sorted_dict = sorted_counter_dict(gene)
        return sorted_dict
    

def _gene_list_printing_statement(sorted_cluster_significances):
    most_common_clusters = list(sorted_cluster_significances.values()) 
    first = max(set(most_common_clusters), key=most_common_clusters.count)
    most_common_clusters = list(filter(lambda a: a != first, most_common_clusters))
    second = max(set(most_common_clusters), key=most_common_clusters.count)
    most_common_clusters = list(filter(lambda a: a != second, most_common_clusters))
    third = max(set(most_common_clusters), key=most_common_clusters.count)
    
    print(f"The genes and their most common clusters are:\n{sorted_cluster_significances}"
            f"The most common cluster is {first}, followed by {second}, followed by {third}")
    
    
def _gene_list_helper(gene=gene, sure=sure, df=df):
    counter==True and gene_list is not None:
    cluster_significances = {}
    for i in gene_list:
        try:
            cluster_significances[i] = list(sorted_counter_dict(gene=i).keys())[0]
        except KeyError:
            continue
    
    sorted_cluster_significances = _sort_dict_by_freq(cluster_significances)
    return _gene_list_printing_statement(sorted_cluster_significances)


# User
def lookup_gene(gene=None, sure=False, gene_list=None, df=data):
    """By default, will see whether input gene is unique,
        - if the input matches more than one gene in list, returns list of genes that match
        - if the input doesn't match any gene in list, tells the user to try again
        - if the input matches one gene in list, returns a pd.Series of clusters and expression levels
        """
    if (gene is not None) and (gene_list is not None):
        return print("Use either \"gene\" or \"gene_list\" but not both!")
    
    elif gene is not None:
        return _gene_expression_helper(gene=gene, sure=sure, df=df)
    
    elif gene_list is not None:
        return _gene_list_helper(gene=gene, sure=sure, df=df)
