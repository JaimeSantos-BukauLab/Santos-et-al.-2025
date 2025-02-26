import numpy as np
import pickle
import pandas as pd
import csv

#Load SeRP Confidence interval files for the two replicates
with open("/path/to/2_loCI_calc/CI_rep1", "rb") as handle:
	BTF3_1 = pickle.load(handle)
with open("/path/to/2_loCI_calc/CI_rep2", "rb") as handle:
	BTF3_2 = pickle.load(handle)

#Additional dictionaries
with open("/path/to/Sequence_human_proteins.pick", "rb") as handle:
	aa_seq = pickle.load(handle)
with open("/path/to/Localization_human_proteins.pick", "rb") as handle:
	loc = pickle.load(handle)

#Filter for genes with >10RPKMs
ref = pd.read_csv('/path/to/included_genes_10RPKM.csv', delimiter=';')
included_genes = {}
for index, row in ref.iterrows():
    included_genes[row["gene"]] = 'YES'

def find_consecutive_ranges(arr):
    ranges = []
    start = end = arr[0]
    for num in arr[1:]:
        if num == end + 1:
            end = num
        else:
            ranges.append((start, end))
            start = end = num
    ranges.append((start, end))

    return ranges


for g in included_genes:
                        if g not in aa_seq:
                                continue
                        if g not in included_genes.keys():
                                continue
#Filter by location (i.e: Cytoplasmatic and Nuclear)
                        if loc[g] != 'CP' and loc[g] != 'NC':
                                continue
                        t = np.mean([BTF3_1[g][3][:, 3], BTF3_2[g][3][:, 3]], axis=0)
                        x = np.where(t >= 0.6)[0]
                        if len(x) != 0:
                                consecutive_ranges = find_consecutive_ranges(x)
                                filtered_ranges = {}
                                filtered_range = [(start, end) for start, end in consecutive_ranges if (end - start) > 9]
                                if len(filtered_range) != 0:
                                   print (g,';',loc[g],';',filtered_range,';',len(t))
