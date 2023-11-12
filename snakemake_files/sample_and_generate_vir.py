import numpy as np
import subprocess
import os
import random
num_samples = 60

def call(s):
    print(s)
    subprocess.run(s.split())

def count_bases_in_fasta(file_path):
    total_bases = 0

    with open(file_path, 'r') as fasta_file:
        for line in fasta_file:
            line = line.strip()  # Remove leading/trailing whitespace
            if not line.startswith('>'):  # Skip header lines
                total_bases += len(line)

    return total_bases

def illumina_sample(reference, abundance, outpath, num_bases_total):
    bases = count_bases_in_fasta(reference)
    num_reads = int(num_bases_total * abundance / 300)
    call(f"wgsim -e 0.01 -1 150 -2 150 -r 0 -N {num_reads} {reference} {outpath}_R1.temp_short {outpath}_R2.temp_short")
    return bases

def sample(input_file, x, input_file2 = None, y= None):
    toret = []
    with open(input_file, 'r') as f:
        lines = f.readlines()
        random.shuffle(lines)
        toret = toret + lines[:x]

    if input_file2 != None:
        with open(input_file2, 'r') as f:
            lines = f.readlines()
            random.shuffle(lines)
            toret = toret + lines[:y]
    toret = [x.strip().split('\t') for x in toret]

    return toret

if len(snakemake.input) == 1:
    samples = sample(snakemake.input[0], num_samples)
else:
    samples = sample(snakemake.input[0], int(num_samples/4), snakemake.input[1], int(num_samples * 3 / 4))

abundances = np.random.lognormal(1, 2, num_samples)
print(abundances)
total_sum = sum(abundances)
abundances = [x/total_sum for x in abundances]
print(abundances)
params = snakemake.params
bases = []

for i in range(num_samples):
    file = "/mnt/disks/1tb/" + samples[i][0]
    call(f"gzip -d {file}")
    file_id = file.rstrip().split('/')[-1]
    ab = abundances[i]
    r = illumina_sample(file[0:-3], ab, params[0] + "/" + file_id, params[1])
    bases.append(r)
    call(f"gzip {file[0:-3]}")

out1 = snakemake.output[0]
out2 = snakemake.output[1]
out3 = snakemake.output[2]

with open(out3, 'w') as f:
    for i in range(num_samples):
        ref = samples[i][0]
        file = samples[i][-1].rstrip()
        ab = abundances[i]
        r = bases[i]
        cov  = params[1]  * ab / r
        f.write(f"{file}\t{ab}\t{r}\t{cov}\t{ref}\n")

os.system(f"cat /mnt/disks/1tb/simulation_reads/" + f"{params[0]}/*R1.temp_short >" + f" /mnt/disks/1tb/simulation_reads/" +f"{out1[0:-3]}")
os.system(f"cat /mnt/disks/1tb/simulation_reads/" + f"{params[0]}/*R2.temp_short >" + f" /mnt/disks/1tb/simulation_reads/" +f"{out2[0:-3]}")
call(f"pigz -p 20  /mnt/disks/1tb/simulation_reads/" + f"{out1[0:-3]}")
call(f"pigz -p 20 /mnt/disks/1tb/simulation_reads/" + f"{out2[0:-3]}")
os.system(f"rm /mnt/disks/1tb/simulation_reads/" + f"{params[0]}/*.temp_short")
os.system(f"rm /mnt/disks/1tb/simulation_reads/" + f"{params[0]}/*temp_long*")
