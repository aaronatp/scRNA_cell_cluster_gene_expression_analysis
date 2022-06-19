import lookup_file
import pandas as pd

def diff_expression(gene, cluster1, cluster2):
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
    

def compare_to_rest(cluster1, cluster2, thresh, n_keep, df=data):
    with HiddenPrints():
        diff_expressed = most_diff_exp(cluster1, cluster2, n_keep, df=data_treat[:100])[0]
    
    clusters_keep = {}
    for i in diff_expressed:
        keep = dict(filter(lambda elem: elem[1] > thresh, lookup_gene(gene=i, sure=True).items()))
        if len(keep) == 1:
            clusters_keep[i] = list(keep)[0]
    
    print(f"\nThe list of genes that are differentialy expressed between cluster 1 and cluster 2," 
            f"that are not expressed in more than {thresh} of any other cluster, are:")
    return clusters_keep

    
def effect(dict1, dict2, dict_entry):
    return dict1[f"{dict_entry}"] - dict2[f"{dict_entry}"]

def compare_treat_ctrl(dictionary, cluster_num):
    treat_effect = {}
    with HiddenPrints():
        for i in list(dictionary.keys()):
            ctrl = expression_distribution(gene=i, cluster=str(cluster_num) + '_CTRL')
            treat = expression_distribution(gene=i, cluster=str(cluster_num) + '_TREAT')
            treat_effect[i] = {'CTRL Median': ctrl['Median'], 'Median effect': effect(treat, ctrl, 'Median'), 'Pct of cells effect': effect(treat, ctrl, 'Pct of cells')}

    return treat_effect


def expression_distribution(gene, cluster):
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
    with HiddenPrints():
        cell_number = data_dict[cluster]
        median = expression_distribution(gene=gene, cluster=cluster)['Median']
        low = expression_distribution(gene=gene, cluster=cluster)['Low']
    
    print(f"{gene} has {cell_number} cells, and has a median expression of {median} transcripts per cell. "
         f"There are at least {low} transcripts in each cell")
    
    a = np.hstack(list(lookup_gene(gene=gene, sure=True, raw=True)[cluster]))
    _ = plt.hist(a, bins='auto')
    plt.show()
