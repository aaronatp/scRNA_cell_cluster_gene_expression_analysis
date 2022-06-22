import lookup_gene
import numpy as np
import pandas as pd

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
with open('/Users/aaronpresser/files_from_quest/Filesfrom_Quest/Peng_RNA_seq/Gene_expression_analysis/diff_exp_treat_6_0_gbm_myeloid.txt','w') as d:
    d.write(str(diff_exp_treat_6_0))
    

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
