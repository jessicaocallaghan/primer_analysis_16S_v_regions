#blast command to make database from 97_otus file
makeblastdb -in 97_otus.fasta -dbtype nucl
#run blastn to the 97 percent OTUs, exclude at 10E-5, take top hit
blastn -db 97_otus.fasta -query seqs.fna -outfmt 6 -num_alignments 1 -evalue .00001 -out blastn_97.out
#run parser.py
python parser.py blastn_97.out blastn_97_with_tax.out
