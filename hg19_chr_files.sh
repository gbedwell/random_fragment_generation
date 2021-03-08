wget http://hgdownload.cse.ucsc.edu/goldenPath/hg19/bigZips/chromFa.tar.gz
tar zvfx chromFa.tar.gz

#If a single concatenated fasta file is ever desired, run:
#cat *.fa > hg19.fa
#NOTE: a single concatenated file will NOT work with the fragment generation python script
