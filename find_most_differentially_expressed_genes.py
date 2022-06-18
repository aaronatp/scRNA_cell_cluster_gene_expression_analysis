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
    

def most_diff_exp(cluster1, cluster2, n_keep, df=data):
    most_diff_expressed = {}
    value_error = []
    values = [0]
    for i in df.index:
        try:
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
  
  if __name__ == 'main':
    most_diff_exp(25, 6, 50)
    
    
