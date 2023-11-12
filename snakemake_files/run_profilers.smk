import pandas as pd
import glob 

samples = pd.read_table("samples.tsv").set_index("samples", drop=False)
#samples = pd.read_table("mock_samples.tsv").set_index("samples", drop=False)
print(samples)

reads_f = [f"{sample}_R1.fq.gz" for sample in samples.index]
reads_r = [f"{sample}_R2.fq.gz" for sample in samples.index]

ganon_db = "/mnt/disks/1tb/benchmark/r89/r89-ganon"
kraken_db = "/mnt/disks/1tb/benchmark/r89/r89-kraken"
kmcp_db = "/mnt/disks/1tb/benchmark/r89/r89-kmcp.kmcp"
r89 = "/mnt/disks/1tb/benchmark/r89"
shell("mkdir -p ganon_out/species90")
shell("mkdir -p ganon_out/species96")
shell("mkdir -p ganon_out/species98")
shell("mkdir -p kraken_out/species90")
shell("mkdir -p kraken_out/species96")
shell("mkdir -p kraken_out/species98")
shell("mkdir -p kmcp_out/species90")
shell("mkdir -p kmcp_out/species96")
shell("mkdir -p kmcp_out/species98")
shell("mkdir -p benchmark/species90")
shell("mkdir -p benchmark/species96")
shell("mkdir -p benchmark/species98")



#kmcp search -d r89/r89-kmcp.kmcp -o kmcp-test-result.tsv.gz -1 ../simulation_reads/species98/species98_6_R1.fq.gz -2 ../simulation_reads/species98/species98_6_R2.fq.gz -j 50^C
#jimshawster@tokyo-240gb /m/d/1/simulation_reads> kmcp profile -X r89/ -T r89/r89-taxid.map -m 3 kmcp-test-result.tsv.gz  -o kmcp-test.profile -C kmcp-test.c.profile -s sample^C                  (taxprof)
#jimshawster@tokyo-240gb /m/d/1/simulation_reads> ganon classify -d r89/r89-ganon -p ../simulation_reads/species98/species98_6_R1.fq.gz ../simulation_reads/species98/species98_6_R2.fq.gz -o ganon-test -t 60 -a^C

rule metaphlan:
    input:
      expand("metaphlan_out/{sample}.metaphlan", sample=samples.index)

rule motus:
    input:
      expand("motus_out/{sample}.motus", sample=samples.index)

rule sylph:
  input:
    "sylph_out/sylph_output.tsv"

rule kraken:
    input:
      expand("kraken_out/{sample}.kraken", sample=samples.index)

rule kmcp:
    input:
      expand("kmcp_out/{sample}.kmcp", sample=samples.index)

rule ganon:
    input:
      expand("ganon_out/{sample}.ganon.tre", sample=samples.index)


rule sylph_all:
  input:
    reads_f,
    reads_r
  output:
    "sylph_out/sylph_output.tsv"
  params:
    reads_f,
    reads_r
  benchmark:
    "benchmarks/sylph.benchmark"
  run:
    shell("sylph sketch -1 {params[0]} -2 {params[1]} -d sylph_sketches -t 50")
    shell("sylph profile sylph_sketches/* /mnt/disks/1tb/benchmark/r89/r89-c200.sylqueries -t 50 > {output}")

rule braken_all:
  input:
    "{sample}_R1.fq.gz",
    "{sample}_R2.fq.gz",
  output:
    "kraken_out/{sample}.kraken"
  params:
    kraken_db,
  benchmark:
    "benchmarks/{sample}.kraken.benchmark"
  run:
    db = params[0]
    rf = input[0]
    rb = input[1]
    shell("kraken2 --paired --db {db} --threads 50 {rf} {rb} --gzip-compressed --report kraken_out/{wildcards.sample}.kreport > temp")
    shell("bracken -d {db} -i kraken_out/{wildcards.sample}.kreport -o {output} -r 150")


rule kmcp_all:
  input:
    "{sample}_R1.fq.gz",
    "{sample}_R2.fq.gz",
  output:
    "kmcp_out/{sample}.kmcp"
  params:
    kmcp_db,
    r89
  benchmark:
    "benchmarks/{sample}.kmcp.benchmark"
  run:
    db = params[0]
    rf = input[0]
    rb = input[1]
    shell("kmcp search -d {params[0]} -o {wildcards.sample}.kmcp.tsv.gz -1 {rf} -2 {rb} -j 50")
    shell("kmcp profile -X {params[1]} -T {params[1]}/r89-taxid.map -m 3 {wildcards.sample}.kmcp.tsv.gz -o {output} -C {output}.profile -s {wildcards.sample}")

rule ganon_all:
  input:
    "{sample}_R1.fq.gz",
    "{sample}_R2.fq.gz",
  output:
    "ganon_out/{sample}.ganon.tre"
  params:
    ganon_db,
    "ganon_out/{sample}.ganon"
  benchmark:
    "benchmarks/{sample}.ganon.benchmark",
  run:
    db = params[0]
    rf = input[0]
    rb = input[1]
    shell("ganon classify -d {params[0]} -p {rf} {rb} -o {params[1]} -t 50 ")

rule metaphlan_all:
  input:
    "{sample}_R1.fq.gz",
    "{sample}_R2.fq.gz",
  output:
    "metaphlan_out/{sample}.metaphlan"
  benchmark:
    "benchmarks/{sample}.metaphlan.benchmark"
  run:
    shell("metaphlan {input[0]},{input[1]} --bowtie2out metaphlan_out/{wildcards.sample}.bowtie2.bz2 --nproc 50 --input_type fastq -o {output}")

rule motus_all:
  input:
    "{sample}_R1.fq.gz",
    "{sample}_R2.fq.gz",
  output:
    "motus_out/{sample}.motus"
  benchmark:
    "benchmarks/{sample}.motus.benchmark"
  run:
      rf = input[0]
      rb = input[1]
      shell("motus profile -f {rf} -r {rb} -t 50 -A -o {output}")

