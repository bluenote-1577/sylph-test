import pandas as pd
import numpy as np
import glob

number_of_bases = 3_000_000_000
species_98_list = "species_98_list.tsv"
species_96_list = "species_96_list.tsv"
species_90_list = "species_90_list.tsv"
iters = [0,1,2,3,4,5,6,7,8,9]
#iters = [0]

rule all:
  input:
    expand("species98/species98_{it}_nano.tsv", it = iters),
    expand("species96/species96_{it}_nano.tsv", it = iters),
    expand("species90/species90_{it}_nano.tsv", it = iters)

rule generate_species98:
  input:
    species_98_list
  output:
    "species98/species98_{iter}_nano.fq.gz",
    "species98/species98_{iter}_nano.tsv"
  params:
    "species98",
    number_of_bases,
    True

  script:
    "sample_and_generate.py"

rule generate_species96:
  input:
    species_96_list
  output:
    "species96/species96_{iter}_nano.fq.gz",
    "species96/species96_{iter}_nano.tsv"
  params:
    "species96",
    number_of_bases,
    True
  script:
    "sample_and_generate.py"

rule generate_species90:
  input:
    species_96_list,
    species_90_list,
  output:
    "species90/species90_{iter}_nano.fq.gz",
    "species90/species90_{iter}_nano.tsv"
  params:
    "species90",
    number_of_bases,
    True

  script:
    "sample_and_generate.py"

