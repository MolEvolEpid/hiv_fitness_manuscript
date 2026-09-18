library(ape)

web_aln <- read.dna('HIV1_FLT_2022_6225-7786_DNA.fasta', format = 'fasta')
web_aln <- del.colgapsonly(web_aln[grepl('^B\\.', rownames(web_aln)),])
con_aln <- read.dna('HIV1_CON_2021_6225-7786_DNA.fasta', format = 'fasta')
con_aln <- del.colgapsonly(con_aln[grepl('^CON_B', rownames(con_aln)),])

write.FASTA(web_aln, 'HIV1_FLT_2022_6225-7786_DNA_B.fasta')
write.FASTA(con_aln, 'HIV1_CON_2021_6225-7786_DNA_B.fasta')

system("mafft --add HIV1_CON_2021_6225-7786_DNA_B.fasta --reorder HIV1_FLT_2022_6225-7786_DNA_B.fasta > HIV1_FLT_CON_2022_6225-7786_DNA_B.fasta")
