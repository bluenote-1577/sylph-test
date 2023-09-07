# sylph-test - scripts for recreating results from sylph paper

In the scripts/* folder you will find a set of python scripts for recreating all figures from our paper. All intermediate results are already present in this directory except for the following:

1. If you want to rerun the MWAS, look at the `analysis.ipynb` notebook. You will require sylphs results for all Parkinson's samples against UHGG, which can be found at https://zenodo.org/record/8151279.
2. If you want to rerun the scripts/manhat.py script, you will need to unzip the gzipped files in this directory.

For example, `python scripts/chng_between_sample.py` will regenerate the figures for the skin microbiome fungi result in our paper. 

3. If you want to look at the results for the 45 sample CPE K. pneumoniae tracking, look at `kang_analysis.ipynb`.
## Requirements

1. Python 3 with the following libraries:
* matplotlib  
* numpy
* seaborn
* natsort
* scipy
* statsmodels
2. jupyter notebook if you want to rerun the MWAS results

