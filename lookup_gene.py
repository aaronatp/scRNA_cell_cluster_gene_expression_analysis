import numpy as np
import pandas as pd
import matplotlib.pyplot as plt

import sys, os

from collections import Counter

# Read in and preprocess data
data = pd.read_csv('/Users/aaronpresser/files_from_quest/Filesfrom_Quest/Peng_RNA_seq/Gene_expression_analysis/seurat_integrated.csv')
cells_to_clusters = pd.read_csv('/Users/aaronpresser/files_from_quest/Filesfrom_Quest/Peng_RNA_seq/Gene_expression_analysis/cells_to_clusters.csv')

data.columns = ['Genes'] + list(cells_to_clusters['V1'])
data.set_index("Genes", inplace=True)

cell_idents = pd.read_csv('/Users/aaronpresser/files_from_quest/Filesfrom_Quest/Peng_RNA_seq/Gene_expression_analysis/integrated_cell_idents.csv')
cell_idents = list(cell_idents['V1'])
cell_clusters = list(data.columns)
new_cols = []
for a,b in zip(data.columns, cell_idents):
    a = str(a) + '_' + b
    new_cols.append(a)

data.columns = new_cols

data_ctrl = data.filter(regex='CTRL$', axis=1)
data_treat = data.filter(regex='TREAT$', axis=1)

# Create dict with cell clusters as keys and number of cells in each as values
data_dict = dict(sorted(dict(Counter(new_cols)).items()))
length = set(key.split('_')[0] for key in data_dict.keys())
cells_to_clusters['Cell_total'] = 0
for i in range(len(length)):
    cells_to_clusters['Cell_total'][str(i) + '_CTRL'] = data_dict[str(i) + '_CTRL']
    cells_to_clusters['Cell_total'][str(i) + '_TREAT'] = data_dict[str(i) + '_TREAT']

    
class HiddenPrints:
    def __enter__(self):
        self._original_stdout = sys.stdout
        sys.stdout = open(os.devnull, 'w')

    def __exit__(self, exc_type, exc_val, exc_tb):
        sys.stdout.close()
        sys.stdout = self._original_stdout    

        
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


def sort_dict_by_freq(store):
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


def lookup_gene(gene=None, sure=False, gene_list=None, df=data, counter=False, raw=False):
    """By default, will see whether input gene is unique,
        - if the input matches more than one gene in list, returns list of genes that match
        - if the input doesn't match any gene in list, tells the user to try again
        - if the input matches one gene in list, returns a pd.Series of clusters and expression levels
        
        If counter=True,
        - returns a dict of clusters ordered by how many cells are in each"""
    if gene is not None:
        options = df.filter(like=f'{gene}', axis=0).index
        if len(options) > 1 and sure==False:
            list_of_options = '\n'.join(i for i in options)
            print(f"There are multiple genes that match your search:")
            return list_of_options
        if len(options) == 0:
            return print("Try again")
        if sure == True and raw == True:
            return data.loc[f'{gene}'][data.loc[f'{gene}'] > 0]
        if len(options) == 1 or sure==True:
#             if proportion==True:
            sorted_dict = sorted_counter_dict(gene)
#             else:
#                 ### Currently setting 'proportion == False' returns 'AttributeError: 'bool' object has no attribute 'loc'
#                 sorted_dict = sorted_counter_dict(gene, proportion==False)
#             if __name__ == 'main':
#                 print(f"The cell clusters that express {gene}, sorted according to how many cells in the cluster express {gene}, are:")
            return sorted_dict
    
    if counter==False and gene_list is not None:
        return print("In order to input a gene list, set 'counter == True'")
    if counter==True and gene_list is not None:
        cluster_significances = {}
        for i in gene_list:
            try:
                cluster_significances[i] = list(sorted_counter_dict(gene=i).keys())[0]
            except KeyError:
                continue
        cluster_significances = sort_dict_by_freq(cluster_significances)
        
        most_common_clusters = list(cluster_significances.values()) 
        first = max(set(most_common_clusters), key=most_common_clusters.count)
        most_common_clusters = list(filter(lambda a: a != first, most_common_clusters))
        second = max(set(most_common_clusters), key=most_common_clusters.count)
        most_common_clusters = list(filter(lambda a: a != second, most_common_clusters))
        third = max(set(most_common_clusters), key=most_common_clusters.count)
        
        return print(f"The genes and their most common clusters are:\n{cluster_significances}\n\n\
The most common cluster is {first}, followed by {second}, followed by {third}")
    
