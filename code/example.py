import numpy as np
import pandas as pd

from collections import Counter

# Read in and preprocess data
data = pd.read_csv('./Gene_expression_analysis/seurat_integrated.csv')
cells_to_clusters = pd.read_csv('./Gene_expression_analysis/cells_to_clusters.csv')

data.columns = ['Genes'] + list(cells_to_clusters['V1'])
data.set_index("Genes", inplace=True)

cell_idents = pd.read_csv('./Gene_expression_analysis/integrated_cell_idents.csv')
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
