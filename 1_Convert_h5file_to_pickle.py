#!/usr/local/bin/python3

#import atomium
import pickle
import os
import sys
import matplotlib.pyplot as plt
import h5py
import numpy as np
import re
from collections import Counter

path = "/path/to/h5files/"
out_path = path+"out_pickle/"
with open("/path/to/Sequence_human_proteins.pick", "rb") as handle:
	ref = pickle.load(handle)

for file_name in ['Selectome_rep1.h5', 'TotalTranslatome_rep1.h5']:
	f = h5py.File(path + file_name, 'r')
	CDS_reads = 0
	for gene in f:
		CDS_reads += np.sum(f[gene][1])
	print(CDS_reads)
	periodicity = [0, 0, 0]
	meta_start, meta_end = {}, {}
	g_count = 0
	out0, out1, out2, out3 = {}, {}, {}, {}
	for gene in f:
		raw_reads = np.zeros(int(f[gene].attrs['cds_length'])+3)
		raw_reads[f[gene][0]] = f[gene][1]

		rpm_reads = (raw_reads*1e6)/CDS_reads
		n = np.mean(raw_reads)
		if n > 0:
			g_count+=1
			for p, i in enumerate(rpm_reads):
				meta_start.setdefault(p, [])
				meta_start[p].append(i)
			for p, i in enumerate(rpm_reads[::-1]):
				meta_end.setdefault(p, [])
				meta_end[p].append(i)

		if len(raw_reads)%3 == 0:
			raw_reads_codon = np.sum(raw_reads.reshape(-1, 3), axis=1)
			rpm_reads_codon = np.sum(rpm_reads.reshape(-1, 3), axis=1)
		else:
			raw_reads_codon = np.sum(np.append(raw_reads, np.zeros(3-len(raw_reads)%3)).reshape(-1, 3), axis=1)
			rpm_reads_codon = np.sum(np.append(rpm_reads, np.zeros(3-len(raw_reads)%3)).reshape(-1, 3), axis=1)
		out0[gene] = raw_reads_codon[:-1]
		out1[gene] = rpm_reads_codon[:-1]

		gene_length = int(f[gene].attrs['cds_length']/3)
		sum_raw_read = np.sum(raw_reads)
		rpkm = (sum_raw_read * 1e9) / (CDS_reads * gene_length)
		rpm = (sum_raw_read * 1e6) / CDS_reads
		out2[gene] = [rpkm, rpm, sum_raw_read, gene_length]

		try:
			peri_h = raw_reads.reshape((-1, 3))
		except:
			continue
		peri_m = np.sum(peri_h, axis=0)
		periodicity[0]+=peri_m[0]
		periodicity[1]+=peri_m[1]
		periodicity[2]+=peri_m[2]

	####output statistics####

	out3['periodicity'] = [periodicity/sum(periodicity), periodicity]
	meta_start_n = {}
	for i, j in meta_start.items():
		if i < 500:
			meta_start_n[i] = np.mean(j)
	out3['meta_start'] = meta_start_n
	meta_end_n = {}
	for i, j in meta_end.items():
		if i < 500:
			meta_end_n[-i] = np.mean(j)
	out3['meta_end'] = meta_end_n
	out3['all_used_reads_in_CDS'] = CDS_reads
	out3['genes_in_metagene'] = g_count

	c, d = 0, 0
	for i, j in ref.items():
		if i in out1.keys():
			continue
		else:
                        if len(j) != 0:
                                out0[i] = np.zeros(len(j))
                                out1[i] = np.zeros(len(j))
                                out2[i] = [0, 0, 0, len(j)]

	for X, n in zip([out0, out1, out2, out3], ['_raw_codon_P15', '_rpm_codon_P15', '_rpkm', '_stats']):
	#for X, n in zip([out0, out1], ['_raw_codon_A12', '_rpm_codon_A12']):
		with open(out_path + file_name + n + ".pick", "wb") as handle:
			pickle.dump(X, handle, protocol=pickle.HIGHEST_PROTOCOL)

	print(file_name, 'done')
