# Concatenate random site fasta files together
cat random_hits_*.fa > random_hits_cat.fa

# Cleanup the final random site file
## Convert fasta file into a two-column tsv - assumes that the fasta file is single line!
## Remove lines whose sequences contain N's
## Remove lines with empty sequence fields
## Convert back into fasta file

awk 'BEGIN{RS=">";OFS="\t"}NR>1{print $1,$2}' random_hits_cat.fa | grep -v 'N' | awk  '$2 != ""' | awk '{print ">"$1"\n"$2}' > cleaned_random_hits.fa

# To better mimic fragment library sizes for NGS, filter the fragments for fragments between 150 bp and 800 bp

awk 'BEGIN{RS=">"}{print $1"\t"$2;}' cleaned_random_hits.fa | awk '{if (length($2) <= 900) {print}}' | awk '{if (length($2) >= 150) {print}}' | awk '{print ">"$1"\n"$2}' > cleaned_filtered.fa

# Generate list of sequences for size distribution analysis

awk 'BEGIN{RS=">";OFS="\t"}NR>1{print $1,$2}' cleaned_filtered.fa | grep -v 'N' | awk  '$2 != ""' | awk '{OFS = "\t"}{print $2}' > random_sequences.txt

# For longer reads, STAR struggles.
# To overcome this, you can re-compile STAR using the 'STARlong' parameters.
# Alternatively, trimming the reads is a viable alternative for our purposes.
# Unique mapping should still be possible for the vast majority of reads.

awk 'BEGIN{RS=">";OFS="\t"}NR>1{print $1,$2}' cleaned_filtered.fa | cut -c -300 | awk '{print ">"$1"\n"$2}' > cleaned_filtered_trimmed.fa
