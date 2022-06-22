from lookup_gene import *
from differential_expression import *
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt

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
      

def format_string(string):
    split_string = string.split('. ')
    counter = 0
    new_list = []
    for a, b in enumerate(split_string):
        new_list.append(b)
        counter += 1
        if counter == 1:
            new_list.append('\n')
            counter = 0
    
    output_string = ' '.join(i for i in new_list)  
    
    return output_string

    
def plotting_statement(gene, cluster_num, group_num):
    groups = ['_CTRL', '_TREAT']
    with HiddenPrints():
        cell_number = data_dict[str(cluster_num) + groups[group_num]]
        median = expression_distribution(gene=gene, cluster=str(cluster_num) + groups[group_num])['Median']
        low = expression_distribution(gene=gene, cluster=str(cluster_num) + groups[group_num])['Low']
    
    return format_string(f"{gene} is expressed {cell_number} cells in cluster {str(cluster_num) + groups[group_num]}. It has a median expression of {median} transcripts per cell. There are at least {low} transcripts in each cell\n")


def plot_expression_distribution(gene, cluster_num):
    '''Returns plot of distribution of gene transcripts levels for 'gene' in 'cluster' '''
    groups = ['_CTRL', '_TREAT']
    
    counter = 0
    try:
        data1 = np.hstack(list(lookup_gene(gene=gene, sure=True, raw=True)[str(cluster_num) + groups[0]]))
        counter += 1
        data2 = np.hstack(list(lookup_gene(gene=gene, sure=True, raw=True)[str(cluster_num) + groups[1]]))
        
    except KeyError:
        return print(f"Can't find {gene} in both clusters")
    
    except TypeError:
        return print(f"{gene} is only recorded to be expressed in one cell in {str(cluster_num) + groups[counter]} unfortunately")
    
    fig, ax = plt.subplots(1,2, sharex=True, figsize=(15,5))
    plt1 = ax[0].hist(data1, bins='rice')
    plt2 = ax[1].hist(data2, bins='rice')
    
    ax[0].text(0.5,-0.4, plotting_statement(gene, cluster_num, group_num=0), size=12, ha='center',
             transform=ax[0].transAxes)
    ax[1].text(0.5,-0.4, plotting_statement(gene, cluster_num, group_num=1), size=12, ha="center", 
             transform=ax[1].transAxes)
    
    return
