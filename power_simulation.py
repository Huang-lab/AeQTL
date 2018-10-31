import argparse
import numpy as np
from scipy import stats
import statsmodels.api as sm
import statsmodels.formula.api as smf
import matplotlib.pyplot as plt

# simulate gene expression 
# 90% WT: 	Normal(20, 10)
# 10% MUT: 	5% Normal(20, 10)
#			5% Nomral(20, 10) - Normal(5, 2)

# linear regression
def sim_single(effect_size, N, SNPs):
	MUT_e = np.random.normal(WT_exp, var_frac*WT_exp, SNPs*N*MUT_perc*eff_perc) + np.random.normal(effect_size, abs(var_frac*effect_size), SNPs*N*MUT_perc*eff_perc)
	if len(MUT_e) == 0:
		MUT_n = np.random.normal(WT_exp, var_frac*WT_exp, SNPs*N*MUT_perc)
	else:
		MUT_n = np.random.normal(WT_exp, var_frac*WT_exp, SNPs*N*MUT_perc*(1-eff_perc))
	MUT = np.append(MUT_e, MUT_n)
	WT = np.random.normal(WT_exp, var_frac*WT_exp, N-SNPs*N*MUT_perc)

	y = np.append(MUT, WT)
	
	x = np.zeros(len(MUT)+len(WT))
	
	# x[0] = 1
	x[0:(N*MUT_perc)] = 1
	
	slope, intercept, r_value, p_value, std_err = stats.linregress(x,y)

	return p_value


def sim_grouped(effect_size, N, SNPs):
	# MUT_effective & MUT_noneffective
	MUT_e = np.random.normal(WT_exp, var_frac*WT_exp, SNPs*N*MUT_perc*eff_perc) + np.random.normal(effect_size, abs(var_frac*effect_size), SNPs*N*MUT_perc*eff_perc)
	
	#### CHANGE TO RNG ####
	if len(MUT_e) == 0:
		MUT_n = np.random.normal(WT_exp, var_frac*WT_exp, SNPs*N*MUT_perc)
	else:
		MUT_n = np.random.normal(WT_exp, var_frac*WT_exp, SNPs*N*MUT_perc*(1-eff_perc))
	#######################

	MUT = np.append(MUT_e, MUT_n)
	WT = np.random.normal(WT_exp, var_frac*WT_exp, N-SNPs*N*MUT_perc)
	MUT_x = np.ones(len(MUT))
	WT_x = np.zeros(len(WT))
	x = np.append(MUT_x, WT_x)
	y = np.append(MUT, WT)
	slope, intercept, r_value, p_value, std_err = stats.linregress(x,y)
	return p_value


def test_single():
	print("Testing single...")
	# single SNPs
	single_avg_pvalue = []
	single_avg_significant = [] 	# avg. percentage of significant pvalue = power
	for i in range(len(N)):
		print(N[i])
		results = []
		significant = 0
		for j in range(num_iter):
			pvalue = sim_single(effect_size, N[i], SNPs)
			if pvalue < 0.05:
				significant += 1
			results.append(pvalue)
		# print(sum(results)/len(results))
		# print(significant/len(results))
		single_avg_pvalue.append(sum(results)/len(results))
		single_avg_significant.append(significant/len(results))

	print("P-value (single):", single_avg_pvalue)
	print("Power (single):", single_avg_significant)
	return single_avg_pvalue, single_avg_significant


def test_grouped():
	print("Testing grouped...")
	# grouped SNPs
	avg_pvalue = []
	avg_significant = [] 	# avg. percentage of significant pvalue
	for i in range(len(N)):
		print(N[i])
		results = []
		significant = 0
		for j in range(num_iter):
			pvalue = sim_grouped(effect_size, N[i], SNPs)
			if pvalue < 0.05:
				significant += 1
			results.append(pvalue)
		avg_pvalue.append(sum(results)/len(results))
		avg_significant.append(significant/len(results))

	print("P-value (grouped):", avg_pvalue)
	print("Power (grouped):", avg_significant)
	return avg_pvalue, avg_significant


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


def plot_sample_size(log_N, single_avg_significant, grouped_avg_significant):
	# effect size = -20
	# Power (single): [0.50980000000000003, 0.5131, 0.50349999999999995, 0.51090000000000002, 0.51239999999999997, 0.72470000000000001, 0.96289999999999998, 0.99880000000000002, 1.0, 1.0]
	# Power (grouped): [0.3268, 0.5126, 0.6434, 0.7526, 0.8317, 0.9792, 1.0, 1.0, 1.0, 1.0]

	# effect size = -10
	# Power (single): [0.19620000000000001, 0.20300000000000001, 0.19159999999999999, 0.20200000000000001, 0.1966, 0.30769999999999997, 0.59719999999999995, 0.85440000000000005, 0.98719999999999997, 1.0]
	# Power (grouped): [0.1227, 0.1849, 0.2328, 0.3033, 0.358, 0.5877, 0.9298, 0.9981, 1.0, 1.0]

	plt.plot(log_N, single_avg_significant, 'o-', label="Single Variant")
	plt.plot(log_N, grouped_avg_significant, 'o-', label="Grouped Variants")
	plt.xlabel("Sample Size (log10)")
	plt.ylabel("Power")
	plt.title("Power vs. Sample Size")
	plt.legend(loc='upper left')

	plt.show()


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
	effect_size = -20
	effect_size = -10
	var_frac = 0.5
	# effect_size = np.arange(-5, -55, -5)

	## x-axis: sample_size, effect_size = -20
	# sample size N
	N = np.array([200, 400, 600, 800, 1000, 2000, 5000, 10000, 20000, 50000])
	# i.e. # effective MUT = [1, 2, 3, 4, 5, 10, 25, 50, 100, 250]
	log_N = np.log10(N)

	num_iter = 10000

	# Begin simulation
	avg_sig_single, avg_sig_grouped = run(N, num_iter, effect_size, SNPs)
	plot_sample_size(log_N, avg_sig_single, avg_sig_grouped)















