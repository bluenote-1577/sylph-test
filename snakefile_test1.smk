import glob
import numpy as np

#input_folder = "./pipeline_test_reads/"
#input_folder = "./mock_community_reads/"
#input_folder = "./gut_reads/"
#input_reads_path = glob.glob(input_folder + "*.fastq.gz")
#input_read_file = [x.split('/')[-1] for x in input_reads_path]

iters = [str(x) for x in range(20)]
cs = [str(x) for x in [100, 500, 1000]]
ks = [str(x) for x in [21, 31]]

input_genome = "./refs/MN-03.fa"
input_genome_95 = "./refs/k.africana.fa"
input_genome_85 = "./refs/k.aerogenes.fa"
output_folder = "./synthetic_simulated_reads/"
output_results_folder = "./synthetic_simulated_results_jul3/"
output_results_folder_95 = "./synthetic_simulated_results_95_jul3/"
output_results_folder_85 = "./synthetic_simulated_results_85_jul3/"
abundances= np.geomspace(0.01,3.0, num = 10)
output_reads = expand(output_folder + "{x}-{y}-synthetic.fq", x = iters, y = abundances)
output_results = expand(output_results_folder + "{k}-{c}-{x}-{y}.tsv",k = ks, c = cs,  x = iters, y = abundances)
output_results_95 = expand(output_results_folder_95 + "{k}-{c}-{x}-{y}.tsv",k = ks, c = cs,  x = iters, y = abundances)
output_results_85 = expand(output_results_folder_85 + "{k}-{c}-{x}-{y}.tsv",k = ks, c = cs,  x = iters, y = abundances)

rule subsamp:
    input:
        output_results,
        output_results_95,
        output_results_85,

rule generate_reads:
    input:
        input_genome
    output:
        output_reads
    run:
        for it in iters:
            for ab in abundances:
                shell(f"~/scratch/software/art_bin_GreatSmokyMountains/art_illumina -ss HS25 -i refs/MN-03.fa -l 150 -f {ab} -o synthetic_simulated_reads/{it}-{ab}-synthetic -s 1 -na -qs 120")

rule contain:
    input:
        output_reads
    output:
        output_results
    run:
        for it in iters:
            for ab in abundances:
                for k in ks:
                    for c in cs:
                        shell(f"sylph contain {input_genome} synthetic_simulated_reads/{it}-{ab}-synthetic.fq -c {c} -k {k} -m 0.75 > {output_results_folder}/{k}-{c}-{it}-{ab}.tsv -t3")


rule contain_95:
    input:
        output_reads
    output:
        output_results_95
    run:
        for it in iters:
            for ab in abundances:
                for k in ks:
                    for c in cs:
                        shell(f"sylph contain {input_genome_95} synthetic_simulated_reads/{it}-{ab}-synthetic.fq -c {c} -k {k} -m 0.75 > {output_results_folder_95}/{k}-{c}-{it}-{ab}.tsv -t3")


rule contain_85:
    input:
        output_reads
    output:
        output_results_85
    run:
        for it in iters:
            for ab in abundances:
                for k in ks:
                    for c in cs:
                        shell(f"sylph contain {input_genome_85} synthetic_simulated_reads/{it}-{ab}-synthetic.fq -c {c} -k {k} -m 0.75 > {output_results_folder_85}/{k}-{c}-{it}-{ab}.tsv -t3")


