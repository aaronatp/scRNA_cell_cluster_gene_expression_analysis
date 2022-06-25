from lookup_gene import *
from plotting_functions import *
import numpy as np
import pandas as pd

# User
def diff_expression(gene, cluster1, cluster2):
    '''Finds difference in proportion of cells that express gene between two clusters
            
        - Returns: 
            - -2 if neither cluster expresses gene
            - Difference between two proportions if gene is expressed in at least one cluster
            ''' 
    dict1 = lookup_gene(gene=gene, sure=True)
    counter = 0
    try:
        dict1[cluster1]
    except KeyError:
        dict1[cluster1] = 0
        counter += 1
    try:
        dict1[cluster2]
    except KeyError:
        dict1[cluster2] = 0
        counter += 1
    if counter == 2: 
        return -2
    else:
        return dict1[cluster1] - dict1[cluster2]
    
# User
def most_diff_exp(cluster1, cluster2, n_keep, df=data_treat):
    '''Finds most differentially expressed genes between two clsuters
        
        - Returns tuple: 
            - return[0] is dict whose keys are genes and values are differences in expression
            for gene between two clusters 
            - return[1] is list of genes that caused a ValueError
            '''
    most_diff_expressed = {}
    value_error = []
    values = [0]
    for i in df.index:
        try:
            with HiddenPrints():
                if diff_expression(i, cluster1, cluster2) == -2:
                    print(f"{i} broke")
                    continue
                else:
                    diff_exp = diff_expression(i, cluster1, cluster2)
                    print(f"Differential expression of {i} = {diff_exp}")
                    for a, b in enumerate(values):
                        if len(most_diff_expressed) == n_keep and diff_exp < b:
                            continue
                        else:
                            # Insert gene:diff_exp element into dict
                            most_diff_expressed = {k: v for k, v in (list(most_diff_expressed.items())[:a] + [(i, diff_exp)] + list(most_diff_expressed.items())[a:])}
                            # Sort dict according to size of values - should be redundant for above line but isn't - fix above line
                            most_diff_expressed = {k: v for k, v in sorted(most_diff_expressed.items(), key=lambda item: item[1], reverse=True)}

                            if len(most_diff_expressed) > n_keep: 
                                most_diff_expressed.pop(list(most_diff_expressed.keys())[-1])
                            values = most_diff_expressed.values()
        
        except ValueError:
            print(f"{i} caued a ValueError")
            value_error.append(i)
    
    print(f"\nThe most differentially expressed genes are:\n{most_diff_expressed}")
    print("\n\n\n")
    print(f"The list of genes that caused a value error are:\n{value_error}")
    
    return most_diff_expressed, value_error
  
# Save output of most_diff_exp('6_TREAT', '0_TREAT', 100)
with open('./Gene_expression_analysis/diff_exp_treat_6_0_gbm_myeloid.txt','w') as d:
    d.write(str(diff_exp_treat_6_0))
    
# User
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

# Helper
def effect(dict1, dict2, dict_entry):
    return dict1[f"{dict_entry}"] - dict2[f"{dict_entry}"]

# User
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
