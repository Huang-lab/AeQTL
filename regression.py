import argparse
import numpy as np
from scipy import stats
import matplotlib.pyplot as plt
import statsmodels.api as sm
import statsmodels.formula.api as smf
import pandas as pd
import os

def argument_parser():
	parser = argparse.ArgumentParser()
	parser.add_argument('-b', '--bin_file_dir', required=True)
	parser.add_argument('-cn', '--covariate_numerical', nargs="+")
	parser.add_argument('-cc', '--covariate_categorical', nargs="+")
	parser.add_argument('-s', '--covariate_file_path')

	args = vars(parser.parse_args())

	bin_file_dir = args['bin_file_dir']
	covariate_numerical = args['covariate_numerical']
	covariate_categorical = args['covariate_categorical']
	covariate_file_path = args['covariate_file_path']

	return bin_file_dir, covariate_numerical, covariate_categorical, covariate_file_path


def add_covariate(bin_file, covariate_numerical, covariate_categorical, source_df):
	bin_file_cov = os.path.splitext(bin_file)[0]+"_cov.txt"

	with open(bin_file_cov, "w") as out:
		with open(bin_file, "r") as f:
			for line in f:
				if line.startswith("genotype"):
					out.write(line.strip("\n")+"\t"+"\t".join(covariate_numerical+covariate_categorical)+"\n")
					continue
				out.write(line.strip("\n"))
				sample = line.strip("\n").split("\t")[1]
				for cov in covariate_numerical:
					source = source_df.loc[source_df['sample_id'] == sample, cov].values
					if len(source) == 1: # make sure sample_id is unique
						out.write("\t"+str(source[0]))
					else:
						out.write("\tNA")
				for cov in covariate_categorical:
					source = source_df.loc[source_df['sample_id'] == sample, cov].values
					if len(source) == 1: # make sure sample_id is unique
						out.write("\t"+str(source[0]))
					else:
						out.write("\tNA")
				out.write("\n")

	return bin_file_cov


def multiregress(df, tested_gene, covariate_numerical, covariate_categorical):
	# Fit the model
	tested_gene = tested_gene.replace('-', '_')
	df.columns = [c.replace('-', '_') for c in df.columns]
	myformula = tested_gene + ' ~ genotype'
	for cov in covariate_numerical:
		myformula += ' + ' + cov
	for cov in covariate_categorical:
		df[cov+'_cat'] = pd.Categorical(df[cov]).codes
		myformula += " + " + cov+"_cat"
	model = smf.ols(formula = myformula, data=df).fit()

	# Return coefficients and p-values of all factors: intercept, genotype, and covariates
	return list(model.params), list(model.pvalues)

	

def main(bin_file_dir, covariate_numerical, covariate_categorical, covariate_file_path):
	source_df = pd.read_table(covariate_file_path, sep="\t", header=0)
	print("reading covariate file --done")

	if bin_file_dir[-1] == '/':
		summary_file = bin_file_dir[0:len(bin_file_dir)-1] + ".summary.txt"
	else:
		summary_file = bin_file_dir + ".summary.txt"

	with open(summary_file, 'w') as out:
		out.write("\t".join(["region","gene","coef_intercept","coef_genotype"] \
			+ ["coef_" + cov for cov in covariate_numerical] \
			+ ["coef_" + cov for cov in covariate_categorical] \
			+ ["pvalue_intercept", "pvalue_genotype"] \
			+ ["pvalue_" + cov for cov in covariate_numerical] \
			+ ["pvalue_" + cov for cov in covariate_categorical] ) \
			+ "\n")
		for bin_file in sorted(os.listdir(bin_file_dir)):
			region = os.path.splitext(bin_file)[0]
			bin_file_cov = add_covariate(bin_file_dir+bin_file, covariate_numerical, covariate_categorical, source_df)
			df = pd.read_table(bin_file_cov, header=0)
			for col in list(df):
				if col == list(covariate_numerical+covariate_categorical)[0]:
					break
				if col == "genotype" or col == "sample_id":
					continue
				gene = col
				out.write(region+"\t"+gene+"\t")
				coefs, pvalues = multiregress(df, gene, covariate_numerical, covariate_categorical)
				# Write coefficients and p-values (5 digits after decimal point) for each factor
				# REF: https://stackoverflow.com/questions/23418600/python-convert-a-list-of-float-to-string
				out.write("\t".join(['{:.5f}'.format(x) for x in coefs] \
					+ ['{:.5f}'.format(x) for x in pvalues]) + "\n")


if __name__ == '__main__':
	bin_file_dir, covariate_numerical, covariate_categorical, covariate_file_path = argument_parser()
	main(bin_file_dir, covariate_numerical, covariate_categorical, covariate_file_path)
