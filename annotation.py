import pysam
import argparse
import vcf
from bx_interval_tree.intersection import IntervalTree, Interval
import numpy as np
from scipy import stats
import matplotlib.pyplot as plt
import subprocess


def argument_parser():
	parser = argparse.ArgumentParser()
	parser.add_argument('-v', '--vcf_path', required=True)
	parser.add_argument('-b', '--bed_path', required=True)
	parser.add_argument('-e', '--expression_matrix_path', required=True)
	parser.add_argument('-o', '--output_dir', required=True)

	args = vars(parser.parse_args())

	vcf_path = args['vcf_path']
	bed_path = args['bed_path']
	expression_matrix_path = args['expression_matrix_path']
	output_dir = args['output_dir']

	# if output_dir[-1] != '/':
	# 	output_dir = output_dir + '/'

	# if not os.path.isdir(output_dir):
	# 	mkdircmd = 'mkdir ' + output_dir
	# 	subprocess.call(shlex.split(mkdircmd))

	return vcf_path, bed_path, expression_matrix_path, output_dir


def get_total_vcf_sample(vcf_path):
	all_samples = set()

	vcf = pysam.VariantFile(vcf_path)
	for record in vcf.fetch():
		for sample in record.samples:
			all_samples.add(sample)
		break

	return all_samples


def get_total_exp_sample(expression_matrix_path):
	all_samples = set()

	with open(expression_matrix_path, "r") as f:
		samples = f.readline().strip("\n").split("\t")[1:]
		for s in samples:
			all_samples.add(s)

	return all_samples


def get_matched_sample(vcf_path, expression_matrix_path):
	vcf_all_samples = get_total_vcf_sample(vcf_path)
	exp_all_samples = get_total_exp_sample(expression_matrix_path)
	
	total_vcf_s = len(vcf_all_samples)
	total_exp_s = len(exp_all_samples)

	matched_sample = set()

	for s in vcf_all_samples:
		if s in exp_all_samples:
			matched_sample.add(s)

	total_matched_s = len(matched_sample)

	return total_vcf_s, total_exp_s, total_matched_s


def build_interval_trees(bed_path):
	# Build a tree for each chromosome
	# trees[0-21]: chr1 - chr22
	# trees[22]: chrX
	# trees[23]: chrY
	trees = []
	for i in range(0, 24):
		trees.append(IntervalTree())

	# Generate RG_dict
	# key: regionName
	# value: list of genes to be tested
	RG_dict = {}

	with open(bed_path, "r") as f:
		for line in f:
			cols = line.strip("\n").split("\t")
			# Insert interval into the tree, Interval [start, end)
			chrom = cols[0]
			start = int(cols[1])
			end = int(cols[2])
			regionName = cols[3]
			testedGenes = cols[4].split(";")
			# Build interval trees
			if chrom == 'X':
				trees[22].insert_interval(Interval(start, end+1, value={'regionName': regionName}))
			elif chrom == 'Y':
				trees[23].insert_interval(Interval(start, end+1, value={'regionName': regionName}))
			else:
				try:					
					index = int(chrom) - 1
					trees[index].insert_interval(Interval(start, end+1, value={'regionName': regionName}))
				except ValueError:
					pass
			# Append gene names to dictionary for each region (key)
			for geneName in testedGenes:
				RG_dict.setdefault(regionName, {})[geneName] = 1

	print("build_interval_trees --done")
	return trees, RG_dict


def map_samples(trees, vcf_path):
	# key: regionName
	# value: list of sample names with variants
	RS_dict = {}

	vcf = pysam.VariantFile(vcf_path)
	for record in vcf.fetch():

		# Determine which tree to query 
		chrom = record.chrom
		if chrom == "X":
			chrom_index = 22
		elif chrom == "Y":
			chrom_index = 23
		else:
			try:
				chrom_index = int(chrom) - 1
			except ValueError:
				pass

		# Determine query interval
		start = int(record.pos)
		SVTYPE = record.info.get("SVTYPE")
		if SVTYPE == None:	# SNP: region = SNP base
			end = start+1
		elif SVTYPE == "INS":	# INS: region = flanking bases
			end = start+2
		elif SVTYPE == "DEL":	# DEL: region = deleted bases
			SVLEN = int(record.info.get("SVLEN"))
			end = start+abs(SVLEN)
		else:
			raise ValueError("Undefined SVTYPE: %s" % SVTYPE)

		# Find overlapping regions
		intervals = trees[chrom_index].find(start, end)
		regions = set()
		for i in intervals:
			regions.add(i.value['regionName'])

		# If variant overlaps with any regions, append sample names to dict of each region
		if (regions != set()):
			# Get sample names with this variant
			names = []
			samples = record.samples
			for i in range(0, len(samples)):
				if samples[i].get("GT") != (None, None) and samples[i].get("GT") != (0, 0):
					names.append(samples[i].name)
			# Append sample names to dictionary for each region (key)
			for regionName in regions:
				for name in names:
					RS_dict.setdefault(regionName, {})[name] = 1
	print("map_samples --done")
	return RS_dict


def bin_expression(RG_dict, RS_dict, expression_matrix_path, output_dir):
	# Write out expression profile for each region
	for region in RG_dict.keys():
		# if there is no variant in this region
		if region not in RS_dict.keys():
			print("WARNING: NO VARIANTS IN REGION \""+region+"\"")
			continue
		MUT_samples = RS_dict[region].keys()
		with open(output_dir+str(region)+".txt",'w') as out:
			out.write("genotype\tsample_id")
			exp_list = []
			with open(expression_matrix_path, "r") as f:
				for i, line in enumerate(f):
					if i == 0:
						samples = line.strip("\n").split("\t")[1:]
					else:
						geneName = line.strip("\n").split("\t")[0]
						if geneName in RG_dict[region].keys():
							exp_list.append(line.strip("\n").split("\t"))
							out.write("\t"+geneName)
			out.write("\n")
			for j in range(0, len(samples)):
				sample = samples[j]
				if sample in MUT_samples:
					out.write("1\t"+sample)	# Genotype = 1 for MUT
				else:
					out.write("0\t"+sample)	# Genotype = 0 for WT
				for k in range(0, len(exp_list)):
					out.write("\t"+exp_list[k][j+1]) # j+1 because first index is geneName
				out.write("\n")
	return 

def main(vcf_path, bed_path, expression_matrix_path, output_dir):
	total_vcf_s, total_exp_s, total_matched_s = get_matched_sample(vcf_path, expression_matrix_path)
	print("VCF sample total: %d, EXP sample total: %d, MATCHED sample total: %d" % (total_vcf_s, total_exp_s, total_matched_s))
	
	trees, RG_dict = build_interval_trees(bed_path)
	RS_dict = map_samples(trees, vcf_path)
	bin_expression(RG_dict, RS_dict, expression_matrix_path, output_dir)


if __name__ == '__main__':
	vcf_path, bed_path, expression_matrix_path, output_dir = argument_parser()
	main(vcf_path, bed_path, expression_matrix_path, output_dir)



