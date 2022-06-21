import lookup_gene, 
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import re

import sys, os


def expression_distribution(gene, cluster):
    '''Returns dict - keys are statistical concepts and dict-values are numerical values for those concepts'''
    distribution = {}
    try:
        expression_levels = list(lookup_gene(gene=gene, sure=True, raw=True)[cluster])
        distribution['Median'] = np.median(expression_levels)
        distribution['Mean'] = np.mean(expression_levels)
        distribution['Std'] = np.std(expression_levels)
        distribution['High'] = np.max(expression_levels)
        distribution['Low'] = np.min(expression_levels)
        distribution['Pct of cells'] = sorted_counter_dict(gene)[cluster]
        
        print(f"The expression distribution of {gene} in cluster {cluster}:")
        return distribution
    
    except KeyError:
        return print(f"{gene} doesn't appear to be found in cluster {cluster}. Perhaps try a different cluster")
    
    except TypeError:
        return print(f"{gene} appears to only be expressed in one cell in this cluster. "
                    "No distribution can be inferred from this data unfortunately")
      

def plot_expression_distribution(gene, cluster):
    '''Returns plot of distribution of gene transcripts levels for 'gene' in 'cluster''''
    with HiddenPrints():
        cell_number = data_dict[cluster]
        median = expression_distribution(gene=gene, cluster=cluster)['Median']
        low = expression_distribution(gene=gene, cluster=cluster)['Low']
    
    print(f"{gene} has {cell_number} cells, and has a median expression of {median} transcripts per cell. "
         f"There are at least {low} transcripts in each cell")
    
    a = np.hstack(list(lookup_gene(gene=gene, sure=True, raw=True)[cluster]))
    _ = plt.hist(a, bins='auto')
    plt.show()


def effect(dict1, dict2, dict_entry):
    return dict1[f"{dict_entry}"] - dict2[f"{dict_entry}"]

def compare_treat_ctrl(dictionary, cluster_num):
    '''Compares expression of genes in treat and control of same cluster
        - Returns dict: keys are genes, values are list of statistical information'''
    treat_effect = {}
    with HiddenPrints():
        for i in list(dictionary.keys()):
            ctrl = expression_distribution(gene=i, cluster=str(cluster_num) + '_CTRL')
            treat = expression_distribution(gene=i, cluster=str(cluster_num) + '_TREAT')
            treat_effect[i] = {'CTRL Median': ctrl['Median'], 'Median effect': effect(treat, ctrl, 'Median'), 'Pct of cells effect': effect(treat, ctrl, 'Pct of cells')}

    return treat_effect


def compare_to_rest(cluster1, cluster2, thresh, n_keep, group=None, df=data_treat):
    '''Returns dict - keys are genes that are differentially expressed between cluster1
        and cluster2, and are not expressed greater than 'thresh' in any other cluster,
        and keys are expression levels'''
    with HiddenPrints():
        diff_expressed = most_diff_exp(cluster1, cluster2, n_keep, df=data_treat)[0]
    
    clusters_keep = {}
    for i in diff_expressed:
        keep = dict(filter(lambda elem: elem[1] > thresh, lookup_gene(gene=i, sure=True).items()))
        if group is not None:
            keep = {k:v for k,v in keep.items() if k.endswith(group)}
        if len(keep) == 1:
            clusters_keep[i] = list(keep)[0]
    
    print(f"\nThe list of genes that are differentialy expressed between cluster {cluster1} and cluster "
          f"{cluster2}, that are not expressed in more than {thresh} of any other cluster, are:")
    return clusters_keep
  
  
brain_genes = pd.read_csv('/Users/aaronpresser/Downloads/genes_expressed_in_brain_proteinatlas.tsv', sep='\t')
human_genes = pd.read_csv('/Users/aaronpresser/Downloads/human_genes_proteinatlas.tsv', sep='\t')

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

def str_strip(string):
    return string.strip()

vec_strip = np.vectorize(str_strip)

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
        
extended_surfaceome = match_surf(synonyms)


def express_lookup_gene(gene):
    '''lookup_gene() shortcut fit for purposes of genes_expressed()'''
    return data.loc[f"{gene}"][data.loc[f"{gene}"] > 0]

def genes_expressed(cluster, data=data_treat):
    """Returns a list of all the genes recorded to be expressed in a particular cluster"""
    genes = []
    for i in data.index:
        try:
            express_lookup_gene(gene=i)[cluster]
            print(f"{i} is expressed in {cluster}")
            genes.append(i)
        except KeyError:
            print(f"{i} isn't expressed\n")
            continue
    return genes

gbm_genes = genes_expressed(cluster='6_TREAT')

gbm_surface_genes = [gene for gene in gbm_genes if gene in extended_surfaceome]

  
