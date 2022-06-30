from lookup_gene import *
from differential_expression import *
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from matplotlib.ticker import MaxNLocator
from statsmodels.distributions.empirical_distribution import ECDF

# User
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
      
# Helper
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

# Helper
def plotting_statement(gene, cluster_num, group_num):
    groups = ['_CTRL', '_TREAT']
    with HiddenPrints():
        cell_number = data_dict[str(cluster_num) + groups[group_num]]
        median = expression_distribution(gene=gene, cluster=str(cluster_num) + groups[group_num])['Median']
        low = expression_distribution(gene=gene, cluster=str(cluster_num) + groups[group_num])['Low']
        pct_cells_expressing = expression_distribution(gene=gene, cluster=str(cluster_num) + groups[group_num])['Pct of cells']
    
    val = float("{:.2f}".format(pct_cells_expressing*100))
    return format_string(f"{gene} is expressed in {int(pct_cells_expressing * cell_number)} of the {cell_number} cells "
                         f"({val}%) in cluster {str(cluster_num) + groups[group_num]}. It has a median expression of "
                         f"{median} transcripts per cell. There are at least {low} transcripts in each cell\n")
# User
def plot_expression_distribution(gene, cluster_num, ecdf=False):
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
    
    if ecdf == True:
        data_groups = [data1, data2]
        markers = ['o', '.']
        labels = ['CTRL', 'TREAT']
        for data, marker, group in zip(data_groups, markers, labels):
            data = ECDF(data)
            plt.plot(data.x, data.y, label=f"{group}", marker=f"{marker}")
        
        plt.legend(loc="lower right")
        plt.xlabel(f'{gene} gene transcripts', fontsize=12)
        return
    
    fig, ax = plt.subplots(1,2, sharex=True, figsize=(15,5))
    plt1 = ax[0].hist(data1, bins='rice')
    plt2 = ax[1].hist(data2, bins='rice')
    
    ax[0].text(0.5,-0.45, plotting_statement(gene, cluster_num, group_num=0), size=12, ha='center',
             transform=ax[0].transAxes)
    ax[0].yaxis.set_major_locator(MaxNLocator(integer=True))
    ax[0].set_xlabel('Gene transcripts', fontsize=12)
    ax[0].set_ylabel('Number of cells', fontsize=12)
    
    ax[1].text(0.5,-0.45, plotting_statement(gene, cluster_num, group_num=1), size=12, ha="center", 
             transform=ax[1].transAxes)
    ax[1].yaxis.set_major_locator(MaxNLocator(integer=True))
    ax[1].set_xlabel('Gene transcripts', fontsize=12)
    
    return
