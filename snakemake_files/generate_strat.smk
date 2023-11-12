import pandas as pd
import numpy as np
import glob

number_of_bases = 750_000_000
anis = [95,96,97,98,99]
lists = expand("strat_{ani}", ani = anis)
iters = [0,1,2,3,4,5,6,7,8,9]
#iters = [0]

rule all:
  input:
    expand("strat/species{ani}_{it}.tsv", ani=anis, it=iters),

rule generate_all:
  input:
    "strat_{ani}.tsv"
  output:
    "strat/species{ani}_{it}_R1.fq.gz",
    "strat/species{ani}_{it}_R2.fq.gz",
    "strat/species{ani}_{it}.tsv"
  params:
    "strat",
    number_of_bases,
    True
  script:
    "sample_and_generate.py"
