#!/bin/sh

./pht.pl --in Teste/sequences/ --cpu 3 --quality_cutoff 33 --quality_percent_cutoff 90 --primers Teste/primers_pcos.csv --len_cutoff 20 --weigth 10 --ref_seqs Teste/refencias_pcos.fasta --depara Teste/nomes_JL.csv --cov 90 --ident_cutoff 90
