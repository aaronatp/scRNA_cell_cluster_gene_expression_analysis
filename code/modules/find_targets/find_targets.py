import helpers
import pandas as pd

surface_genes = pd.read_csv('./Gene_expression_analysis/surface_genes.csv')
surface_genes.rename(columns={'To': 'Genes'}, inplace=True)
surface_genes.drop(columns=['Unnamed: 0', 'From'], inplace=True)  
  
brain_genes = pd.read_csv('./Downloads/genes_expressed_in_brain_proteinatlas.tsv', sep='\t')
human_genes = pd.read_csv('./Downloads/human_genes_proteinatlas.tsv', sep='\t')

synonyms = human_genes['Gene synonym'].apply(_synonym_cleaner)
surfaceome_concat = ' '.join(i for i in surface_genes['Genes'])


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
