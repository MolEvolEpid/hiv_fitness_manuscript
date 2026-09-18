cat HIV1_CON_2021_pol_DNA.fasta hiv-db_gap_squeeze-pol.fasta > pol.fasta
cat HIV1_CON_2021_env_DNA.fasta hiv-db_gap_squeeze-env.fasta > env.fasta

mafft pol.fasta > pol_aligned.fasta
mafft env.fasta > env_aligned.fasta
