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

## Usage

	aeqtl -v <VCF file, cannot be used together with MAF file> \
	      -m <MAF file, cannot be used together with VCF file> \
	      -b <BED file> \
	      -e <Expression file> \
	      -cn <Numerical covariates, separated by ","> \
	      -cc <Categorical covariates, separated by ","> \
	      -s <Covariate file> \
	      -ts <Threshold of mutated sample number per region, default=1> \
	      -tv <Threshold of variant number per region, default=1> \
          -o <output directory>

## Input data format

*Note: demo input files with compatible format can be found in the folder `demo/myfiles`*

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

A tab separated .tsv file of summary statistics. Each row is an eQTL test between a region and a gene. The file contains the following fields:

- region
- gene
- coef_intercept
- coef_genotype
- coef_\<covariate> *(for each covariate)*
- pvalue_intercept
- pvalue_genotype
- pvalue_\<covariate> *(for each covariate)*
- adj_pvalue_intercept
- adj_pvalue_genotype
- adj_pvalue_\<covariate> *(for each covariate)*


## Example using TCGA somatic data

This is a demo analysis of TCGA somatic truncations on chromosome 17 using AeQTL. The whole procedure from data pre-processing to generating the final summary file will be included. 

See folder `demo/tcga_somatic` for input files.

#### 1. Download public TCGA MAF file

Download the TCGA MC3 Public MAF from <https://api.gdc.cancer.gov/data/1c8cfe5f-e52d-41ba-94da-f15ea1337efc> (Ellrott et al., 2018) and decompress the file.
    
    gunzip mc3.v0.2.8.PUBLIC.maf.gz

Extract truncations on chromosome 17. Here we define truncations as variants whose `Variant_Classification` is one of the following: `Frame_Shift_Del`, `Frame_Shift_Ins`, `Nonsense_Mutation`, `Splice_Site`, and `Translation_Start_Site`.
    
    awk '{ if ($5 == 17) { print } }' mc3.v0.2.8.PUBLIC.maf > mc3.v0.2.8.PUBLIC.chr17.maf
    egrep 'Frame_Shift_Del|Frame_Shift_Ins|Nonsense_Mutation|Splice_Site|Translation_Start_Site' mc3.v0.2.8.PUBLIC.chr17.maf > mc3.v0.2.8.PUBLIC.chr17.truncations.maf
    head -n 1 mc3.v0.2.8.PUBLIC.maf | cat - mc3.v0.2.8.PUBLIC.chr17.truncations.maf > tmp && mv tmp mc3.v0.2.8.PUBLIC.chr17.truncations.maf

#### 2. Data pre-processing
*(Note: this part is processed in **R**)*

- **Filter expression file**
    
    Read file (a subset data is used for demoing purpose).

        exp <- read.table("BRCA_RSEM_hugo.matrix.demo.txt", header=TRUE, sep="\t")
    
    Subset file to only include tumor samples.

        columns_to_keep <- c(1, which(substr(colnames(exp),14,14)=="0"))
        exp_tumor_only <- subset(exp, select=columns_to_keep)
    
    Since low expression levels will likely present technical noises, we can also filter out genes with low expression medians across samples. Here we use log(2) as the threshold.
    
        medians <- apply(exp_tumor_only[,2:ncol(exp_tumor_only)], 1, median)
        exp_tumor_only_filtered <- exp_tumor_only[which(log(medians)>2),]
    
    Write out the filtered file.
    
        write.table(exp_tumor_only_filtered,file="BRCA_RSEM_hugo.matrix.demo.tumor.only.log2.txt", quote=FALSE, sep="\t", row.names=FALSE)
    
- **Re-format sample ID in MAF and expression files**

    Sample IDs used in both MAF and expression files are formatted using TCGA Sample Barcode (<https://docs.gdc.cancer.gov/Encyclopedia/pages/TCGA_Barcode/>). Since we match samples on a participant level, we will use the first 12 characters of the barcode as sample ID.
    
        maf <- read.table("mc3.v0.2.8.PUBLIC.chr17.truncations.maf", header=TRUE)
        maf$Tumor_Sample_Barcode <- substr(maf$Tumor_Sample_Barcode, 1, 12)
        exp <- read.table("BRCA_RSEM_hugo.matrix.demo.tumor.only.log2.txt", header=TRUE)
        write.table(maf, file="mc3.v0.2.8.PUBLIC.chr17.truncations.maf", quote=FALSE, sep="\t", row.names=FALSE)
        write.table(exp, file="BRCA_RSEM_hugo.matrix.demo.tumor.only.log2.txt", quote=FALSE, sep="\t", row.names=FALSE)

#### 3. Run AeQTL

Run AeQTL with 6 covariates, including 3 numerical covariates (age, PC1, PC2) and 3 categorical covariates (gender, race, tumor stage).  

    aeqtl -m mc3.v0.2.8.PUBLIC.chr17.truncations.maf \
          -b allCDS_1BasedStart_2bpFlanks.chr17.bed \
          -e BRCA_RSEM_hugo.matrix.demo.tumor.only.log2.txt \
          -cn age PC1 PC2 \
          -cc gender race ajcc_pathologic_tumor_stage \
          -s PanCan_ClinicalData.txt \
          -o tcga_somatic_truncations_output
          
The folder `demo/tcga_somatic/tcga_somatic_truncations_output/` contains all the intermediate files, and see summary file at `demo/tcga_somatic/tcga_somatic_truncations_output.summary.txt`.

