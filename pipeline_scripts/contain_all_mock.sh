#sylph contain nanopore_mock2_1gb.fq.prs reference/all_genomes_mock2/* --ci -c100 > results/mock2_nano_c100
#sylph contain nanopore_mock2_1gb.fq.prs reference/all_genomes_mock2/* --ci -c1000 > results/mock2_nano_c1000
#sylph contain pac_mock2_1gb.fq.prs reference/all_genomes_mock2/* --ci -c100 > results/mock2_pac_c100
#sylph contain pac_mock2_1gb.fq.prs reference/all_genomes_mock2/* --ci -c1000 > results/mock2_pac_c1000
#sylph contain ill_mock2_1gb.fq.prs reference/all_genomes_mock2/* --ci -c100 > results/mock2_ill_c100
#sylph contain ill_mock2_1gb.fq.prs reference/all_genomes_mock2/* --ci -c1000 > results/mock2_ill_c1000

#sylph contain nanopore_mock2_1gb.fq.prs reference/all_genomes_mock2/* --ci -c100 > results_gap/mock2_nano_c100
#sylph contain nanopore_mock2_1gb.fq.prs reference/all_genomes_mock2/* --ci -c1000 > results_gap/mock2_nano_c1000
#sylph contain pac_mock2_1gb.fq.prs reference/all_genomes_mock2/* --ci -c100 > results_gap/mock2_pac_c100
#sylph contain pac_mock2_1gb.fq.prs reference/all_genomes_mock2/* --ci -c1000 > results_gap/mock2_pac_c1000
#sylph contain ill_mock2_1gb.fq.prs reference/all_genomes_mock2/* --ci -c100 > results_gap/mock2_ill_c100
#sylph contain ill_mock2_1gb.fq.prs reference/all_genomes_mock2/* --ci -c1000 > results_gap/mock2_ill_c1000

sylph contain gtdb-r214-c100_gap.pgs pac_mock2_1gb.fq.prs --ci -t 20 -m 0.85 > gtdb-on-reads/gtdb-on-pac-c100.tsv
sylph contain gtdb-r214-c1000_gap.pgs pac_mock2_1gb.fq.prs --ci -t 20 -m 0.85 > gtdb-on-reads/gtdb-on-pac-c1000.tsv
sylph contain gtdb-r214-c100_gap.pgs nanopore_mock2_1gb.fq.prs --ci -t 20 -m 0.85 > gtdb-on-reads/gtdb-on-nano-c100.tsv
sylph contain gtdb-r214-c1000_gap.pgs nanopore_mock2_1gb.fq.prs --ci -t 20 -m 0.85 > gtdb-on-reads/gtdb-on-nano-c1000.tsv
sylph contain gtdb-r214-c100_gap.pgs  ill_mock2_1gb.fq.prs --ci -t 20 -m 0.85 > gtdb-on-reads/gtdb-on-ill-c100.tsv
sylph contain gtdb-r214-c1000_gap.pgs ill_mock2_1gb.fq.prs  --ci -t 20 -m 0.85 > gtdb-on-reads/gtdb-on-ill-c1000.tsv


##mash
#mash screen mock2_mash.msh final_reads/nanopore_mock2_1gb.fq -p 10 > results/mash_nano
#mash screen mock2_mash.msh final_reads/ill_mock2_1gb.fq -p 10 > results/mash_ill
#mash screen mock2_mash.msh final_reads/pac_mock2_1gb.fq -p 10 > results/mash_pac
#
##sourmash
#sourmash search  nanopore_mock2_1gb.fq.sig sourmash_mock2/* --max-contain -t 0.0 -o results/sour_nano
#sourmash search  ill_mock2_1gb.fq.sig sourmash_mock2/* --max-contain -t 0.0 -o results/sour_ill
#sourmash search  pac_mock2_1gb.fq.sig sourmash_mock2/* --max-contain -t 0.0 -o results/sour_pac
#

