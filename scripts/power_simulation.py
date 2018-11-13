import argparse
import numpy as np
from scipy import stats
import statsmodels.api as sm
import statsmodels.formula.api as smf
import matplotlib.pyplot as plt

# simulate gene expression 
# 90% WT: 	Normal(20, 10)
# 10% MUT: 	5% Normal(20, 10)
#			5% Nomral(20, 10) - Normal(effect_size, effect_size/2)

def simulate(effect_size, N, SNPs):
	num_MUT_e = SNPs*N*MUT_perc*eff_perc
	num_MUT_n = SNPs*N*MUT_perc*(1-eff_perc)
	num_WT = N - SNPs*N*MUT_perc

	MUT_e = np.random.normal(WT_exp, var_frac*WT_exp, num_MUT_e) + \
		np.random.normal(effect_size, abs(var_frac*effect_size), num_MUT_e)
	MUT_n = np.random.normal(WT_exp, var_frac*WT_exp, num_MUT_n)
	
	MUT = np.append(MUT_e, MUT_n)
	WT = np.random.normal(WT_exp, var_frac*WT_exp, num_WT)

	y = np.append(MUT, WT)
	
	single_pvalue = []
	if num_MUT_e <= 5:
		# single
		for i in range(int(num_MUT_e)):
			x_single = np.append(np.ones(1), np.zeros(N - 1))
			x_single = sm.add_constant(x_single)
			results = sm.OLS(y, x_single).fit()
			single_pvalue.append(list(results.pvalues)[1])		# first pvalue is for intercept	
	else:
		# single
		num_occurr = N*MUT_perc		# num of occurrence for each variant
		for i in range(int(SNPs*eff_perc)):
			x_single = np.append(np.ones(num_occurr), np.zeros(N - num_occurr))
			x_single = sm.add_constant(x_single)
			results = sm.OLS(y, x_single).fit()
			single_pvalue.append(list(results.pvalues)[1])

	# grouped
	x_grouped = np.append(np.ones(len(MUT)), np.zeros(len(WT)))
	x_grouped = sm.add_constant(x_grouped)
	results = sm.OLS(y, x_grouped).fit()
	grouped_pvalue = list(results.pvalues)[1]

	return single_pvalue, grouped_pvalue


def run(N, num_iter, effect_size, SNPs):
	avg_sig_single = []
	avg_sig_grouped = []

	for i in range(len(N)):
		print(N[i])
		results_single = []
		results_grouped = []
		sig_single = 0
		sig_grouped = 0
		for j in range(num_iter):
			p_single, p_grouped = simulate(effect_size, N[i], SNPs)
			sig_single += (sum([p < 0.05 for p in p_single])/len(p_single))
			if p_grouped < 0.05:
				sig_grouped += 1
		avg_sig_single.append(sig_single/num_iter)
		avg_sig_grouped.append(sig_grouped/num_iter)
	
	print("Power (single):", avg_sig_single)
	print("Power (grouped):", avg_sig_grouped)
	return avg_sig_single, avg_sig_grouped


def plot_sample_size(log_N, single_avg_significant, grouped_avg_significant, N):
	# effect size = -20
	# Power (single): [0.50980000000000003, 0.5131, 0.50349999999999995, 0.51090000000000002, 0.51239999999999997, 0.72470000000000001, 0.96289999999999998, 0.99880000000000002, 1.0, 1.0]
	# Power (grouped): [0.3268, 0.5126, 0.6434, 0.7526, 0.8317, 0.9792, 1.0, 1.0, 1.0, 1.0]

	# effect size = -10
	# Power (single): [0.19620000000000001, 0.20300000000000001, 0.19159999999999999, 0.20200000000000001, 0.1966, 0.30769999999999997, 0.59719999999999995, 0.85440000000000005, 0.98719999999999997, 1.0]
	# Power (grouped): [0.1227, 0.1849, 0.2328, 0.3033, 0.358, 0.5877, 0.9298, 0.9981, 1.0, 1.0]

	plt.plot(log_N, single_avg_significant, 'o-', label="Single Variant")
	plt.plot(log_N, grouped_avg_significant, 'o-', label="Grouped Variants")
	plt.xticks(log_N, N)
	plt.xlabel("Sample Size")
	plt.ylabel("Power")
	plt.title("Power vs. Sample Size")
	plt.legend(loc='upper left')

	plt.show()


def plot_sample_size_2(log_N, avg_sig_single_20, avg_sig_grouped_20, avg_sig_single_10, avg_sig_grouped_10, N):
	#REF: https://matplotlib.org/gallery/showcase/bachelors_degrees_by_gender.html#sphx-glr-gallery-showcase-bachelors-degrees-by-gender-py
	# dark green, light green, dark blue, light blue --> group_10, group_20, single_10, single_20
	color_seq = ['#3D5B0A', '#8eba43', '#003C68', '#2988bc'] 
	fig, ax = plt.subplots(1, 1, figsize=(7, 7))
	# Remove the plot frame lines. They are unnecessary here.
	# ax.spines['top'].set_visible(False)
	# ax.spines['bottom'].set_visible(False)
	# ax.spines['right'].set_visible(False)
	# ax.spines['left'].set_visible(False)
	ax.set_ylim(0, 1.05)
	
	plt.plot(log_N, avg_sig_grouped_10, 'o-', label="Grouped-variants (effect size = -10)", color=color_seq[0])
	plt.plot(log_N, avg_sig_grouped_20, 'o-', label="Grouped-variants (effect size = -20)", color=color_seq[1])
	plt.plot(log_N, avg_sig_single_10, 'o-', label="Single-variant (effect size = -10)", color=color_seq[2])
	plt.plot(log_N, avg_sig_single_20, 'o-', label="Single-variant (effect size = -20)", color=color_seq[3])

	plt.xticks(np.log10([200, 400, 1000, 2000, 5000, 10000, 25000]), [200, 400, 1000, 2000, 5000, 10000, 25000], fontsize=12)
	plt.xlabel("Sample Size", fontsize=18, weight='bold', labelpad=15)
	plt.ylabel("Power", fontsize=18, weight='bold', labelpad=15)
	# plt.title("Power vs. Sample Size", fontsize=20, weight='bold', y=1.05)
	plt.legend(loc='best', frameon=False)
	plt.tight_layout()

	# plt.show()
	plt.savefig("/Users/dgl/eQTL/doc/figures/fig1/combined/myfig1b", dpi=300)

if __name__ == '__main__':
	# Fix number of SNPs = 10. 
	# Assume rare variants so each person will have one variant at most. 
	SNPs = 10
	eff_perc = 0.5

	# When simulate for single variant, do 1 SNP is sufficient (i.e. for each simulation
	# one test for single variant with 1 SNP; one test for grouped variants with 10 SNPs.

	# MUT_perc = frequency of each SNP present in a sample
	MUT_perc = 0.001
	WT_perc = 1 - MUT_perc

	WT_exp = 20
	effect_size_20 = -20
	effect_size_10 = -10
	effect_size = -10
	var_frac = 0.5
	# effect_size = np.arange(-5, -55, -5)

	## x-axis: sample_size, effect_size = -20
	# sample size N
	# N = np.array([200, 400, 600, 800, 1000, 2000, 3000, 4000, 5000, 10000, 20000, 50000])
	N = np.append(np.arange(200, 1000, step=200), np.arange(1000, 26000, step=1000))   
	# i.e. # effective MUT = [1, 2, 3, 4, 5, 10, 25, 50, 100, 250]
	log_N = np.log10(N)

	num_iter = 10000

	# Begin simulation
	# avg_sig_single, avg_sig_grouped = run(N, num_iter, effect_size, SNPs)
	# plot_sample_size(log_N, avg_sig_single, avg_sig_grouped, N)

	# avg_sig_single_20, avg_sig_grouped_20 = run(N, num_iter, effect_size_20, SNPs)
	# avg_sig_single_10, avg_sig_grouped_10 = run(N, num_iter, effect_size_10, SNPs)
	avg_sig_single_20 = [0.50739999999999996, 0.50770000000000004, 0.50960000000000005, 0.49990000000000001, 0.49780000000000002, 0.72050000000000003, 0.84030000000000005, 0.92379999999999995, 0.96109999999999995, 0.98089999999999999, 0.98970000000000002, 0.99650000000000005, 0.99739999999999995, 0.99860000000000004, 0.99970000000000003, 0.99980000000000002, 0.99990000000000001, 1.0, 1.0, 0.99990000000000001, 1.0, 0.99990000000000001, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0]
	avg_sig_grouped_20 = [0.3199, 0.5088, 0.6411, 0.7468, 0.8243, 0.9785, 0.9977, 0.9998, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0]
	avg_sig_single_10 = [0.1968, 0.20030000000000001, 0.1961, 0.20219999999999999, 0.20100000000000001, 0.31809999999999999, 0.41439999999999999, 0.50980000000000003, 0.59509999999999996, 0.65759999999999996, 0.72119999999999995, 0.77170000000000005, 0.82310000000000005, 0.86050000000000004, 0.88870000000000005, 0.90549999999999997, 0.92910000000000004, 0.94450000000000001, 0.95650000000000002, 0.96630000000000005, 0.97340000000000004, 0.97889999999999999, 0.98240000000000005, 0.98640000000000005, 0.99119999999999997, 0.99129999999999996, 0.99450000000000005, 0.99570000000000003, 0.99709999999999999]
	avg_sig_grouped_10 = [0.1238, 0.1807, 0.2392, 0.2986, 0.3524, 0.5968, 0.7587, 0.8669, 0.9329, 0.9642, 0.9813, 0.992, 0.9958, 0.9973, 0.9989, 0.9996, 0.9998, 0.9997, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0]

	plot_sample_size_2(log_N, avg_sig_single_20, avg_sig_grouped_20, avg_sig_single_10, avg_sig_grouped_10, N)






