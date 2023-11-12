import pandas as pd
import numpy as np
import glob

number_of_bases = 3_000_000_0
species_95_list = "vir_spec95_list.txt"
species_90_list = "vir_spec90_list.txt"
#iters = [0,1,2,3,4,5,6,7,8,9,10]
iters = [0,1,2,3,4,5,6,7,8,9,10]
#iters = [0]

rule all:
  input:
    expand("vir_species95/species95_{it}.tsv", it = iters)

rule generate_species95:
  input:
    species_95_list,
  output:
    "vir_species95/species95_{iter}_R1.fq.gz",
    "vir_species95/species95_{iter}_R2.fq.gz",
    "vir_species95/species95_{iter}.tsv"
  params:
    "vir_species95",
    number_of_bases,
    True

  script:
    "sample_and_generate_vir.py"

