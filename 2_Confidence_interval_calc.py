import numpy as np
import pickle
from statsmodels.stats.proportion import proportion_confint

####parameter####
window30 = 30
window45 = 45
#in codons

####iterate over files#####
path = '/path/to/h5files/out_pickle/'
input_files = [['Selectome_rep1', 'TotalTranslatome_rep1'], ['Selectome_rep2', 'TotalTranslatome_rep2']]
output_file_names = ['CI_rep1', 'CI_rep2']

for experiment, out_name in zip(input_files, output_file_names):
	master = {}
	with open(path + str(experiment[0]) + '.h5_raw_codon_P15.pick', 'rb') as f:
		total = pickle.load(f)
	with open(path + str(experiment[1]) + '.h5_raw_codon_P15.pick', 'rb') as f:
		ip = pickle.load(f)
	#count all reads
	outfile_CI, outfile_score = {}, {}
	sum_total = np.sum([np.sum(x) for x in total.values()])
	sum_ip = np.sum([np.sum(x) for x in ip.values()])
	normalization = sum_ip/sum_total
	common_keys= []
	for key in total.keys():
                if key in ip.keys():
                        common_keys.append(key)
	for p, gene in enumerate(common_keys):

		total_gene = total[gene]
		ip_gene = ip[gene]

		total_window30, ip_window30 = [], []
		for pos30 in np.arange(len(total_gene)):
			if pos30 <= np.floor(window30/2):
				start30 = 0
			else:
				start30 = int(pos30-np.floor(window30/2))
			if pos30 > len(total_gene)-np.ceil(window30/2):
				stop30 = len(total_gene)
			else:
				stop30 = int(pos30+np.ceil(window30/2))
			total_window30.append(np.sum(total_gene[start30:stop30]))
			ip_window30.append(np.sum(ip_gene[start30:stop30]))

		CI30 = proportion_confint(ip_window30, np.add(ip_window30, total_window30), alpha=0.05, method='agresti_coull')
		odds_CI_low30 = np.divide(CI30[0], np.subtract(1, CI30[0])) / normalization
		CI_add30 = CI30[1].copy()
		CI_add30[CI_add30 == 1] = 0.9
		odds_CI_high30 = np.divide(CI30[1], np.subtract(1, CI_add30)) / normalization
		odds_CI_mean30 = np.mean([odds_CI_high30, odds_CI_low30], axis=0)


		total_window45, ip_window45 = [], []
		for pos45 in np.arange(len(total_gene)):
			if pos45 <= np.floor(window45/2):
				start45 = 0
			else:
				start45 = int(pos45-np.floor(window45/2))
			if pos45 >= len(total_gene)+np.ceil(window45/2):
				stop45 = len(total_gene)
			else:
				stop45 = int(pos45+np.ceil(window45/2))
			total_window45.append(np.sum(total_gene[start45:stop45]))
			ip_window45.append(np.sum(ip_gene[start45:stop45]))

		CI45 = proportion_confint(ip_window45, np.add(ip_window45, total_window45), alpha=0.05, method='agresti_coull')
		odds_CI_low45 = np.divide(CI45[0], np.subtract(1, CI45[0])) / normalization
		CI_add45 = CI45[1].copy()
		CI_add45[CI_add45 == 1] = 0.9
		odds_CI_high45 = np.divide(CI45[1], np.subtract(1, CI_add45)) / normalization
		odds_CI_mean45 = np.mean([odds_CI_high45, odds_CI_low45], axis=0)

		binding_HC30 = np.where(odds_CI_low30 > 1.5)[0]
		binding_LC30 = np.where((odds_CI_low30 < 1.5) & (odds_CI_high30 > 1.5))[0]
		no_binding30 = np.where(odds_CI_high30 < 1.5)[0]
		binding30 = np.zeros(len(odds_CI_low30))
		binding30[binding_HC30] = 1
		binding30[binding_LC30] = 0.2
		binding30[no_binding30] = 0

		binding_HC45 = np.where(odds_CI_low45 > 1.5)[0]
		binding_LC45 = np.where((odds_CI_low45 < 1.5) & (odds_CI_high45 > 1.5))[0]
		no_binding45 = np.where(odds_CI_high45 < 1.5)[0]
		binding45 = np.zeros(len(odds_CI_low45))
		binding45[binding_HC45] = 1
		binding45[binding_LC45] = 0.2
		binding45[no_binding45] = 0

		if len(odds_CI_low45) > 40:
			max_loCI45 = np.amax(odds_CI_low45[30:-10])
			max_loCI_pos45 = np.argmax(odds_CI_low45[30:-10])+30
		else:
			max_loCI45 = 0
			max_loCI_pos45 = -1

		tag = []
		if len(np.where(CI45[1] == 1)[0]) != 0:
			tag.append('missing values')

		master[gene] = [tag, max_loCI45, max_loCI_pos45, \
					np.array([odds_CI_high30, odds_CI_mean30, odds_CI_low30, binding30]).T, \
					np.array([odds_CI_high45, odds_CI_mean45, odds_CI_low45, binding45]).T]
		if p%1000 == 0:
			print(p)

	with open(path + "2_loCI_calc/" + str(out_name) + '.pick', 'wb') as handle:
		pickle.dump(master, handle, protocol=pickle.HIGHEST_PROTOCOL)
