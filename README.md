# AeQTL

eQTL analysis using region-based aggregation of rare variants. 

## Requirements

- python 3.5
- pip
- bx_interval_tree (see installation instructions below)
- git (optional)

## Installation

First, install IntervalTree from bx-python. We strongly recommend using a standalone package called *bx_interval_tree* which is smaller and easier to compile than bx-python. 

```
git clone https://github.com/ccwang002/bx_interval_tree
cd bx_interval_tree
python setup.py install
cd ..
```

Continue to install AeQTL by choosing one of the options below.

#### (1) From PyPI

The easiest way to install AeQTL is from PyPI.

	pip install aeqtl

#### (2) From source code

Alternatively, download the source code of AeQTL

	git clone https://github.com/Huang-lab/AeQTL

Then install AeQTL

	cd AeQTL
	pip install .

#### Optional (but recommended)

Append the path to AeQTL to your PATH environment variable

	export PATH=/path/to/AeQTL/bin:$PATH

## Run

	aeqtl -v <vcf file> -b <bed file> -e <expression file> \
    	  -cn <numerical covariates> -cc <categorical covariates> -s <covariate file> \
          -o <output directory>

## Input data format

*Note: demo input files with compatible format can be found in the "demo" folder*

#### VCF file

A standard multi-sample VCF file with file extension .vcf (or .vcf.gz). Sample IDs in VCF file, expression file, and covariate file should match exactly. 

#### BED file

A BED file (tab separated) with at least four columns and without header. The format of the file should follow:

	<chromosome>	<start>		<end>		<region_name>		<tested_genes>

An example row:

	chr17			41197693	41197821	BRCA1				BRCA1;SLC25A39;HEXIM2

The first four columns are required. The fifth column is a list of genes separated by ";". If the fifth column (tested_genes) is not provided, AeQTL by default will test each region with every gene from the expression file. 

#### Expression file

A matrix-format, tab separated .tsv file with gene expression from RNA-seq. The first row (header) of the file should follow:

	gene_id		<sample_id_1>		<sample_id_2>		<sample_id_3>		...

and the first column of the file should follow:

	gene_id
	<gene_1>
	<gene_2>
	...

#### Covariate file

A tab separated .tsv file with column names corresponding to covariates. A column of sample IDs with column name "sample_id" is required. Covariates entered in AeQTL and their corresponding column names must match exactly. However, the covariate file can contain other unused columns as well. If entering a categorical covariate, please make sure each category has the same value throughout the file (i.e. avoid instances such as having both "FEMALE" and "female" in the same column).


## Output data format

A tab separated .tsv file of summary statistics (up to 5 digits after the decimal point). Each row is an eQTL test between a region and a gene. The file contains the following fields:

- region
- gene
- coef_intercept
- coef_genotype
- coef_\<covariate> *(for each covariate)*
- pvalue_intercept
- pvalue_genotype
- pvalue_\<covariate> *(for each covariate)*








