import helpers
import pandas as pd

synonyms = human_genes['Gene synonym'].apply(_synonym_cleaner)

def genes_expressed(cluster, data=data_treat):
    '''Returns a list of all the genes recorded to be expressed in a particular cluster'''
    expressed = []
    for i in data.index:
        try:
            _express_lookup_gene(gene=i)[cluster]
            print(f"{i} is expressed in {cluster}")
            expressed.append(i)
        except KeyError:
            print(f"\t{i} isn't expressed")
            continue
    return expressed

gbm_genes = genes_expressed(cluster='6_TREAT')

gbm_surface_genes = [gene for gene in gbm_genes if gene in _match_surf(synonyms)]
