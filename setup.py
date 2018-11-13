#REF: https://packaging.python.org/tutorials/packaging-projects/
import setuptools

with open("README.md", "r") as fh:
	long_description = fh.read()

setuptools.setup(
	name="aeqtl",
	version="0.1.0",
	author="Guanlan Dong",
	author_email="guanlan.dong@wustl.edu",
	description="eQTL analysis using region-based aggregation of rare variants.",
	long_description=long_description,
	long_description_content_type="text/markdown",
	url="https://github.com/Huang-lab/AeQTL",
	classifiers=[
		"Programming Language :: Python :: 3",
		"License :: OSI Approved :: MIT License",
		"Operating System :: OS Independent",
	],
	entry_points={
	'console_scripts': [
			'aeqtl=aeqtl:main',
		],
	},
	package_dir = { 
		'aeqtl':'aeqtl', 
	},
	packages = [ \
		'aeqtl' ,
	],
	install_requires=[
        'pysam==0.13',
        'PyVCF',
        'Cython', # for bx_interval_tree
        'statsmodels',
        'scipy',
        'pandas',
        'bx_interval_tree'
    ],
)
