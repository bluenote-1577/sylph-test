import glob
import numpy as np

#input_folder = "./pipeline_test_reads/"
#input_folder = "./mock_community_reads/"
input_folder = "./gut_reads/"
output_folder = "./output_pipeline_test_reads/"
input_reads_path = glob.glob(input_folder + "*.fastq*")
input_read_file = [x.split('/')[-1] for x in input_reads_path]
iters = [str(x) for x in range(4)]
abundances= np.geomspace(0.01,0.5, num = 10)
output_reads = expand(output_folder + "{x}-{y}-{z}", x = iters, y = abundances, z= input_read_file)

rule subsamp:
    input:
        output_reads

rule do_subsamp:
    input:
        input_reads_path
    output:
        output_reads
    run:
        for it in iters:
            for ab in abundances:
                for (i,read) in enumerate(input_reads_path):
                    #shell(f"seqtk sample {read} {ab} -s {it} | gzip > {output_folder}{it}-{ab}-{input_read_file[i]}")
                    shell(f"seqtk sample {read} {ab} -s {it} > {output_folder}{it}-{ab}-{input_read_file[i]}")
