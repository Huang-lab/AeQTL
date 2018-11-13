'''
main wrapper script
'''

import pysam
import subprocess
import shlex
import argparse
import re
import os
import sys
import annotation
import regression

def argument_parser():
	parser = argparse.ArgumentParser()

	parser.add_argument('-v', '--vcf_path', required=True)
	parser.add_argument('-b', '--bed_path', required=True)
	parser.add_argument('-e', '--expression_matrix_path', required=True)
	parser.add_argument('-cn', '--covariate_numerical', nargs="+")
	parser.add_argument('-cc', '--covariate_categorical', nargs="+")
	parser.add_argument('-s', '--covariate_file_path')
	parser.add_argument('-o', '--output_dir', required=True)

	args = vars(parser.parse_args())

	vcf_path = args['vcf_path']
	bed_path = args['bed_path']
	expression_matrix_path = args['expression_matrix_path']
	covariate_numerical = args['covariate_numerical']
	covariate_categorical = args['covariate_categorical']
	covariate_file_path = args['covariate_file_path']
	output_dir = args['output_dir']

	if output_dir[-1] != '/':
		output_dir = output_dir + '/'

	if not os.path.isdir(output_dir):
		mkdircmd = 'mkdir ' + output_dir
		subprocess.call(shlex.split(mkdircmd))
	else:
		cmd1 = 'rm -rf ' + output_dir
		cmd2 = 'mkdir ' + output_dir
		subprocess.call(shlex.split(cmd1))
		subprocess.call(shlex.split(cmd2))

	return vcf_path, bed_path, expression_matrix_path, covariate_numerical, covariate_categorical, covariate_file_path, output_dir


# e.g.
# Expression matrix: 500 genes and 100 samples
# Genotype file: 300 variants and 95 samples
# Covariate matrix: 4 covariates and 80 samples; 3 covariates used
def file_info(vcf_path, expression_matrix_path, covariate_numerical, covariate_categorical, covariate_file_path):
	vcf = pysam.VariantFile(vcf_path)
	vcf_var = 0
	for record in vcf.fetch():
		vcf_sample = len(record.samples)
		vcf_var += 1

	expression_gene = 0
	with open(expression_matrix_path, 'r') as f:
		samples = f.readline().strip("\n").split("\t")[1:]
		expression_sample = len(samples)
		f.readline()
		for line in f:
			expression_gene += 1

	cov_used_num = len(covariate_numerical) + len(covariate_categorical)
	cov_sample = 0
	with open(covariate_file_path, 'r') as f:
		cov_num = len(f.readline().strip("\n").split("\t"))-1
		for line in f:
			cov_sample += 1

	return vcf_sample, vcf_var, expression_sample, expression_gene, cov_sample, cov_num, cov_used_num


def main():
	vcf_path, bed_path, expression_matrix_path, covariate_numerical, covariate_categorical, covariate_file_path, output_dir = argument_parser()

	# extract basic information from input files
	vcf_sample, vcf_var, expression_sample, expression_gene, cov_sample, cov_num, cov_used_num = file_info(vcf_path, expression_matrix_path, covariate_numerical, covariate_categorical, covariate_file_path)
	print("VCF file:", vcf_sample, "samples and", vcf_var, "variants")
	print("Expression matrix:", expression_sample, "samples and", expression_gene, "genes")
	print("Covariate file:", cov_sample, "samples and", cov_num, "covariates;", cov_used_num, "covariates used")

	# map samples to genes, write bin expression file for linear regression
	annotation.main(vcf_path, bed_path, expression_matrix_path, output_dir)

	# build linear regression model and perform t test, write output summary file
	print('buiding regression model...')
	regression.main(output_dir, covariate_numerical, covariate_categorical, covariate_file_path)

	print('done')
	return 0

if __name__=='__main__':
	main()
