import glob
import numpy as np

input_folder = "/home/jshaw/scratch/2023_prita_test/atopic_chng/"

reads_1 = sorted(glob.glob(input_folder+"*_1*"))
reads_2 = sorted(glob.glob(input_folder+"*_2*"))
output_folder = "/home/jshaw/scratch/2023_zip_test/chng_flashed/"
reads_out = [output_folder + x.replace('_1','').split('/')[-1] for x in reads_1]


rule all:
    input:
        reads_1,
        reads_2
    output:
        reads_out
    run:
        for i in range(len(reads_1)):
            r1 = reads_1[i]
            r2 = reads_2[i]
            rout = reads_out[i]
            shell(f"flash {r1} {r2} -t 10 --min-overlap 20 --max-overlap 100; cat out*.fastq > {rout}")
