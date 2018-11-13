import os
import sys
import pysam
import vcf
from bx_interval_tree.intersection import IntervalTree, Interval
import statsmodels.api as sm
import statsmodels.formula.api as smf
import pandas as pd

class aeqtl(object):
	''' Example usage:
			AeQTL = aeqtl.aeqtl()
			AeQTL.fileinfo()
			AeQTL.annotation()
			AeQTL.regression()
	'''
	def __init__( self , **kwargs ):
		self.vcf_path = kwargs.get( 'vcf_path', "" )
		self.bed_path = kwargs.get( 'bed_path', "" )
		self.expression_matrix_path = kwargs.get( 'expression_matrix_path', "" )
		self.covariate_numerical = kwargs.get( 'covariate_numerical', [] )
		self.covariate_categorical = kwargs.get( 'covariate_categorical', [] )
		self.covariate_file_path = kwargs.get( 'covariate_file_path', "" )
		self.output_dir = kwargs.get( 'output_dir', "" )

	# e.g.
	# Expression matrix: 500 genes and 100 samples
	# Genotype file: 300 variants and 95 samples
	# Covariate matrix: 4 covariates and 80 samples; 3 covariates used
	def file_info(self):
		vcf = pysam.VariantFile(self.vcf_path)
		vcf_var = 0
		for record in vcf.fetch():
			vcf_sample = len(record.samples)
			vcf_var += 1

		expression_gene = 0
		with open(self.expression_matrix_path, 'r') as f:
			samples = f.readline().strip("\n").split("\t")[1:]
			expression_sample = len(samples)
			f.readline()
			for line in f:
				expression_gene += 1

		cov_used_num = len(self.covariate_numerical) + len(self.covariate_categorical)
		cov_sample = 0
		with open(self.covariate_file_path, 'r') as f:
			cov_num = len(f.readline().strip("\n").split("\t"))-1
			for line in f:
				cov_sample += 1

		return vcf_sample, vcf_var, expression_sample, expression_gene, cov_sample, cov_num, cov_used_num
	
	##################
	### ANNOTATION ###
	##################

	def get_total_vcf_sample(self):
		all_samples = set()
		vcf = pysam.VariantFile(self.vcf_path)
		for record in vcf.fetch():
			for sample in record.samples:
				all_samples.add(sample)
			break
		return all_samples


	def get_total_exp_sample(self):
		all_samples = set()
		with open(self.expression_matrix_path, "r") as f:
			samples = f.readline().strip("\n").split("\t")[1:]
			for s in samples:
				all_samples.add(s)
		return all_samples


	def get_matched_sample(self):
		vcf_all_samples = self.get_total_vcf_sample()
		exp_all_samples = self.get_total_exp_sample()
		total_vcf_s = len(vcf_all_samples)
		total_exp_s = len(exp_all_samples)
		matched_sample = set()
		for s in vcf_all_samples:
			if s in exp_all_samples:
				matched_sample.add(s)
		total_matched_s = len(matched_sample)
		return total_vcf_s, total_exp_s, total_matched_s


	def check_tested_genes(self):
		all_genes = []
		f = open(self.bed_path, 'r')
		line = f.readline()
		f.close()
		if len(line.strip("\n").split("\t")) < 4:
			sys.exit("Error: number of columns in BED file is smaller than 4")
		if len(line.strip("\n").split("\t")) == 4:
			with open(self.expression_matrix_path, 'r') as exp:
				for i, l in enumerate(exp):
					if i == 0:
						continue
					all_genes.append(l.strip("\n").split("\t")[0])
		return all_genes


	def build_interval_trees(self):
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

		# Check whether user has provided the column of tested genes in BED file
		all_genes = self.check_tested_genes()

		with open(self.bed_path, "r") as f:
			for line in f:
				cols = line.strip("\n").split("\t")
				# Insert interval into the tree, Interval [start, end)
				chrom = cols[0]
				start = int(cols[1])
				end = int(cols[2])
				regionName = cols[3]
				# If user has provided tested genes
				if len(all_genes) == 0:
					testedGenes = cols[4].split(";")
				# If user has not provided tested genes for ANY region, 
				# by default test each region against every gene in expression file
				else:
					testedGenes = all_genes
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


	def map_samples(self, trees):
		# key: regionName
		# value: list of sample names with variants
		RS_dict = {}

		vcf = pysam.VariantFile(self.vcf_path)
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


	def bin_expression(self, RG_dict, RS_dict):
		# Write out expression profile for each region
		for region in RG_dict.keys():
			# If there is no variant in this region
			if region not in RS_dict.keys():
				print("WARNING: NO VARIANTS IN REGION \""+region+"\"")
				continue
			MUT_samples = RS_dict[region].keys()
			with open(self.output_dir+str(region)+".txt",'w') as out:
				out.write("genotype\tsample_id")
				exp_list = []
				with open(self.expression_matrix_path, "r") as f:
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

	def annotation(self):
		total_vcf_s, total_exp_s, total_matched_s = self.get_matched_sample()
		print("VCF sample total: %d, EXP sample total: %d, MATCHED sample total: %d" % (total_vcf_s, total_exp_s, total_matched_s))
		
		trees, RG_dict = self.build_interval_trees()
		RS_dict = self.map_samples(trees)
		self.bin_expression(RG_dict, RS_dict)

	##################
	### REGRESSION ###
	##################

	def add_covariate(self, bin_file, source_df):
		bin_file_cov = os.path.splitext(bin_file)[0]+"_cov.txt"

		with open(bin_file_cov, "w") as out:
			with open(bin_file, "r") as f:
				for line in f:
					if line.startswith("genotype"):
						out.write(line.strip("\n")+"\t"+"\t".join(self.covariate_numerical+self.covariate_categorical)+"\n")
						continue
					out.write(line.strip("\n"))
					sample = line.strip("\n").split("\t")[1]
					for cov in self.covariate_numerical:
						source = source_df.loc[source_df['sample_id'] == sample, cov].values
						if len(source) == 1: # make sure sample_id is unique
							out.write("\t"+str(source[0]))
						else:
							out.write("\tNA")
					for cov in self.covariate_categorical:
						source = source_df.loc[source_df['sample_id'] == sample, cov].values
						if len(source) == 1: # make sure sample_id is unique
							out.write("\t"+str(source[0]))
						else:
							out.write("\tNA")
					out.write("\n")

		return bin_file_cov


	def multiregress(self, df, tested_gene):
		# Fit the model
		tested_gene = tested_gene.replace('-', '_')
		df.columns = [c.replace('-', '_') for c in df.columns]
		myformula = tested_gene + ' ~ genotype'
		for cov in self.covariate_numerical:
			myformula += ' + ' + cov
		for cov in self.covariate_categorical:
			df[cov+'_cat'] = pd.Categorical(df[cov]).codes
			myformula += " + " + cov+"_cat"
		model = smf.ols(formula = myformula, data=df).fit()

		# Return coefficients and p-values of all factors: intercept, genotype, and covariates
		return list(model.params), list(model.pvalues)

		
	def regression(self):
		source_df = pd.read_table(self.covariate_file_path, sep="\t", header=0)
		print("reading covariate file --done")

		if self.output_dir[-1] == '/':
			summary_file = self.output_dir[0:len(self.output_dir)-1] + ".summary.txt"
		else:
			summary_file = self.output_dir + ".summary.txt"

		with open(summary_file, 'w') as out:
			out.write("\t".join(["region","gene","coef_intercept","coef_genotype"] \
				+ ["coef_" + cov for cov in self.covariate_numerical] \
				+ ["coef_" + cov for cov in self.covariate_categorical] \
				+ ["pvalue_intercept", "pvalue_genotype"] \
				+ ["pvalue_" + cov for cov in self.covariate_numerical] \
				+ ["pvalue_" + cov for cov in self.covariate_categorical] ) \
				+ "\n")
			for bin_file in sorted(os.listdir(self.output_dir)):
				region = os.path.splitext(bin_file)[0]
				bin_file_cov = self.add_covariate(self.output_dir+bin_file, source_df)
				df = pd.read_table(bin_file_cov, header=0)
				for col in list(df):
					if col == list(self.covariate_numerical+self.covariate_categorical)[0]:
						break
					if col == "genotype" or col == "sample_id":
						continue
					gene = col
					out.write(region+"\t"+gene+"\t")
					coefs, pvalues = self.multiregress(df, gene)
					# Write coefficients and p-values (5 digits after decimal point) for each factor
					# REF: https://stackoverflow.com/questions/23418600/python-convert-a-list-of-float-to-string
					out.write("\t".join(['{:.5f}'.format(x) for x in coefs] \
						+ ['{:.5f}'.format(x) for x in pvalues]) + "\n")
