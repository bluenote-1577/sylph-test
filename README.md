# sylph-test - scripts for recreating results from sylph paper

In the scripts/* folder you will find a set of python scripts for recreating results from Fig.1, Fig. 2, Fig.3, Fig. 4, Fig. 6 from our paper. The notebooks will show how Figs. 2 and 5 were obtained. All intermediate results are already present in this directory except for the following:

For example, `python scripts/chng_between_sample.py` will regenerate the figures for the skin microbiome fungi result in our paper. 

For files that are gzipped, you may have to ungzip. 

### Fig. 1

* Look at `scripts/synthetic_pois_plot.py`

### Fig. 2

* Look at `profile_analysis.ipynb`.

### Fig. 3

* `cd` into `real_gut_results_v0.5` or `real_gut_results_v0.5/hadza` and look at the `scripts` python scripts. 

### Fig. 4

* Look at `scripts/diagonal_ani_nn.py` and `scripts/mock_community_plot.py`

### Fig. 5

* Look at `analysis.ipynb` for the MWAS statistical calculation and `scripts/manhat.py`.
* Look at `scripts/chng_between_sample.py` for the Chng et al analysis

### Fig. 6
* Look at `real_profile_analysis.ipynb`.

1. Python 3 with the following libraries:
* matplotlib  
* numpy
* seaborn
* natsort
* scipy
* statsmodels
2. jupyter notebook if you want to rerun the MWAS results

