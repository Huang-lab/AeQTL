import sys
import os
import subprocess
import shlex
from aeqtl import aeqtl
import argparse

def argument_parser():
	parser = argparse.ArgumentParser()

	parser.add_argument('-v', '--vcf_path', required=True, help="Path to VCF file")
	parser.add_argument('-b', '--bed_path', required=True, help="Path to BED file")
	parser.add_argument('-e', '--expression_matrix_path', required=True, help="Path to expression file")
	parser.add_argument('-cn', '--covariate_numerical', nargs="+", help="Name(s) of numerical covariate(s); space separated")
	parser.add_argument('-cc', '--covariate_categorical', nargs="+", help="Name(s) of categorical covariate(s); space separated")
	parser.add_argument('-s', '--covariate_file_path', help="Path to covariate file")
	parser.add_argument('-o', '--output_dir', required=True, help="Path to write output files")

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


def main():
	vcf_path, bed_path, expression_matrix_path, covariate_numerical, covariate_categorical, covariate_file_path, output_dir = argument_parser()
	AeQTL = aeqtl.aeqtl(vcf_path=vcf_path, bed_path=bed_path, expression_matrix_path=expression_matrix_path,
		covariate_numerical=covariate_numerical, covariate_categorical=covariate_categorical, 
		covariate_file_path=covariate_file_path, output_dir=output_dir)

	# extract basic information from input files
	vcf_sample, vcf_var, expression_sample, expression_gene, cov_sample, cov_num, cov_used_num = AeQTL.file_info()
	print("VCF file:", vcf_sample, "samples and", vcf_var, "variants")
	print("Expression matrix:", expression_sample, "samples and", expression_gene, "genes")
	print("Covariate file:", cov_sample, "samples and", cov_num, "covariates;", cov_used_num, "covariates used")

	# map samples to genes, write bin expression file for linear regression
	AeQTL.annotation()

	# build linear regression model and perform t test, write output summary file
	print('buiding regression model...')
	AeQTL.regression()

	print('done')
	return 0















