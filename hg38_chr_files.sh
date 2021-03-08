wget https://hgdownload-test.gi.ucsc.edu/goldenPath/hg38/bigZips/hg38.chromFa.tar.gz
tar zvfx hg38.chromFa.tar.gz

#If a single concatenated fasta file is ever desired, run:
#cat chroms/chr*.fa > hg38.fa
#NOTE: a single concatenated file will NOT work with the fragment generation python script
