#!/usr/local/bin/python
######################################################################################################
# This script generates a required number of random fragments from the human genome.
# It can accomidate a mix of restriction enzymes (e.g. NheI/AvrII/SpeI/BamHI or MseI/BglII) for fragmentation.
# It can also fragment by theoretical random fragmentation (e.g. sonication).
#
# It is best to use this script to generate a single dataset of random integration site fragments to match wet lab experiment site numbers.
# For large random fragment datasets, it is best to use the appropriate iterative script that changes the random seed every N sites.
#
# The genome is defined in base_counts.txt.
# base_counts is a tsv file with four columns:
# 1: All chromosome names for the genome build (e.g. http://hgdownload.cse.ucsc.edu/goldenPath/hg19/chromosomes/ for hg19).
# 2: Chromosome sizes (e.g. https://hgdownload.cse.ucsc.edu/goldenPath/hg19/bigZips/hg19.chrom.sizes for hg19).
# 3: Relative (global) beginning of each chromosome within entire human genome.
# 4: Relative (global) end of each chromosome within entire human genome.
#
# The last row, col [4] is easily calculated as sum of all chromosome lengths
# The last row col [3] is calculated as [4] - [2].
# The other columns follow logically.
#
#
# Created December 2006, 2014 by Peter Cherepanov Peter.Cherepanov@crick.ac.uk
# Modified May 2020 by Greg Bedwell gregoryjbedwell@gmail.com
#
# Written for python 2.7. Migration to python 3 is underway.
######################################################################################################

import os, string, random, re, sys
from string import maketrans, translate, join
from itertools import islice
from random import randint

print "\n\n"
print "*"*100
print "Simulates random fragmentation of the human genome"
print "*"*100
print "\n\n"
print "Requires:\n"
print "\t(1) list of chromosomes\n"
print "\t(2) corresponding chromosome fasta files\n"
print "\t(3) chromosome bp counts file (base_counts.txt)\n\n" #see information and links above
print "Program writes out a joined fasta file (random_hits.fa) containing the generated fragments\n\n" #random_hits.txt currently commented out
print "-"*100
print "\n\n"
BUILD = raw_input("Enter the the genome build to use (hg19 or hg38)\n")
print "\n"

if (BUILD != 'hg19' and BUILD != 'hg38'):
    sys.exit('ERROR: Please enter a valid genome build identifier\n\n')

PATH = raw_input("Enter the path to the chromosome fasta files\n")
print "\n"
N = int(raw_input("Enter the number of random fragments to generate:\n")) #number of required random integration sites
print "\n"
FRAG = raw_input("Enter the genome fragmentation method:\n\nMB for MseI/BglII digestion, NASB for NheI/AvrII/SpeI/BamHI digestion, Random for random fragmentation\n")

if (FRAG != 'MB' and FRAG != 'NASB' and FRAG != 'Random' ):
    sys.exit('ERROR: Please enter a valid genome fragmentation method\n\n')


print "\n"
DISTANCE = int(raw_input("Enter the maximum allowed distance to cleavage site\n")) #maximum distance to cleavage site. I use 20000 for all fragmentation methods
print "\n"

if FRAG == 'Random':
    mu = int(raw_input("Enter the mean fragment length for random fragmentation\n")) #I use mu = 400
    print "\n"
    sigma = int(raw_input("Enter the standard deviation of the fragment length for random fragmentation\n")) #I use sigma = 50.
    print "\n"

print "-"*100
path = PATH    #location of the  chromosome fasta files and base_counts file
trans = maketrans("acgttumrwsykvhdbnACGTTUMRWSYKVHDBN", "TGCAAAKYWSRMBDHVNTGCAAAKYWSRMBDHVN")  #reverse complement rules

def reversecomp(sequence):
    sequence_list = list(sequence)
    sequence_list.reverse()
    reversed_sequence = join(sequence_list, '')
    reverse_complement_sequence = translate(reversed_sequence, trans)
    return reverse_complement_sequence

#random_site_list = open('random_hits.txt', 'w')
random_site_fasta = open('random_hits.fa', 'w')
random.seed()    #set the seed based on the current system time

if BUILD == 'hg19':
    i = 0
    while i < N:
        current_site = random.randint(1,3137161264) #set range = the number of bp in genome build. currently defined for hg19.
        #print current_site
        counts_file = open(path + '/base_counts_hg19.txt', 'r') #base_counts file should be named
        for line in counts_file.readlines():
            line_clean = line.replace('\n', '')
            foo = line_clean.split('\t')
            if (current_site > int(foo[2]) and current_site < int(foo[3])): #if current site > chr start (foo[2]) and < chr end (foo[3])
                bingo_chromosome = foo[0] #define chromosome containing integration site
                bingo_site = current_site - int(foo[2])
                #print bingo_chromosome, bingo_site, '\n'
                chromosome = open(path + '/' + foo[0] + '.fa', 'r') #open chromosome file containing bingo_site
                #print chromosome.readline()

                #here adding random orientation of the provirus....
                strand = random.randrange(2) #range can take on two possible values - 0 and 1.
                if strand == 0: #strand == 0 is + sense
                    chromosome_fragment = chromosome.read().replace('\n', '')[bingo_site - 1:bingo_site + DISTANCE - 1]  #limit = DISTANCE
                    #print "strand=", strand
                    #print chromosome_fragment
                else:
                    chromosome_fragment = chromosome.read().replace('\n', '')[bingo_site - DISTANCE - 1:bingo_site - 1]
                    #print "strand=", strand, "reversing..."
                    #print chromosome_fragment, '\n'
                    chromosome_fragment = reversecomp(chromosome_fragment)
                    #print chromosome_fragment, '\n'
                #print restre.search(chromosome_fragment)

    #fragment DNA by method of choice

                if FRAG == 'MB':
                    restre = re.compile('(TTAA|AGATCT)', re.IGNORECASE)  #define restriction site sequences, case-insensitive!
                    if  restre.search(chromosome_fragment):
                        #print restre.findall(chromosome_fragment)[0]
                        site = chromosome_fragment.find(restre.findall(chromosome_fragment)[0])
                        #print "cropped sequence=", str(site + 1), "bp"
                        if site > 14: #We want minimum 18 bp, including 2bp of RE site. This is for mapping purposes. Most fragments will be much larger.
                            #print chromosome_fragment[0:site + 2]   #will take first 2 bp of the restriction site.
                            entry1 = bingo_chromosome.replace('.fa', '') + '\t' + str(bingo_site) + '\t' + 'random_site_number_' + str(i) + '\t' + chromosome_fragment[0:site + 2] + '\n'
                    #random_site_list.write(entry1)
                            entry2 = '>' + bingo_chromosome.replace('.fa', '') + '_' + 'strand_' + str(strand) + '_' + str(bingo_site) + '_' + 'random_site_number_' + str(i) + '\n' + chromosome_fragment[0:site + 2] +'\n'
                            random_site_fasta.write(entry2)
                            i = i + 1
                            print 'generated site number  ', i

                if FRAG == 'NASB':
                    restre = re.compile('(GCTAGC|CCTAGG|ACTAGT|GGATCC)', re.IGNORECASE)  #pattern for searching the REsites, case-insensitive!
                    if  restre.search(chromosome_fragment):
                        #print restre.findall(chromosome_fragment)[0]
                        site = chromosome_fragment.find(restre.findall(chromosome_fragment)[0])
                        #print "cropped sequence=", str(site + 1), "bp"
                        if site > 14:  #We want minimum 18 bp, including 2bp of RE site. This is for mapping purposes. Most fragments will be much larger.
                            #print chromosome_fragment[0:site + 2]   #will take first 2 bp of the restriction site.
                            entry1 = bingo_chromosome.replace('.fa', '') + '\t' + str(bingo_site) + '\t' + 'random_site_number_' + str(i) + '\t' + chromosome_fragment[0:site + 2] + '\n'
                            #random_site_list.write(entry1)
                            entry2 = '>' + bingo_chromosome.replace('.fa', '') + '_' + 'strand_' + str(strand) + '_' + str(bingo_site) + '_' + 'random_site_number_' + str(i) + '\n' + chromosome_fragment[0:site + 2] +'\n'
                            random_site_fasta.write(entry2)
                            i = i + 1
                            print 'generated site number  ', i

                if FRAG == 'Random':
                    site = int(round(random.gauss(mu,sigma), 1))
                    #site = randint(rand_min, rand_max)
                    #print chromosome_fragment[0:site + 0]
                    #We didn't include a size limit here, as the model distribution strongly disfavors small fragments.
                    #This is borne out is the analysis of fragment lengths.
                    #We don't have that same control over the RE digested fragments.
                    entry1 = bingo_chromosome.replace('.fa', '') + '\t' + str(bingo_site) + '\t' + 'random_site_number_' + str(i) + '\t' + chromosome_fragment[0:site + 0] + '\n'
                    #random_site_list.write(entry1)
                    entry2 = '>' + bingo_chromosome.replace('.fa', '') + '_' + 'strand_' + str(strand) + '_' + str(bingo_site) + '_' + 'random_site_number_' + str(i) + '\n' + chromosome_fragment[0:site + 0] +'\n'
                    random_site_fasta.write(entry2)
                    i = i + 1
                    print 'generated site number  ', i

if BUILD == 'hg38':
    i = 0
    while i < N:
        current_site = random.randint(1,3209286105) #set range = the number of bp in genome build
        #print current_site
        counts_file = open(path + '/base_counts_hg38.txt', 'r') #base_counts file should be named
        for line in counts_file.readlines():
            line_clean = line.replace('\n', '')
            foo = line_clean.split('\t')
            if (current_site > int(foo[2]) and current_site < int(foo[3])): #if current site > chr start (foo[2]) and < chr end (foo[3])
                bingo_chromosome = foo[0] #define chromosome containing integration site
                bingo_site = current_site - int(foo[2])
                #print bingo_chromosome, bingo_site, '\n'
                chromosome = open(path + '/' + foo[0] + '.fa', 'r') #open chromosome file containing bingo_site
                #print chromosome.readline()

                #here adding random orientation of the provirus....
                strand = random.randrange(2) #range can take on two possible values - 0 and 1.
                if strand == 0: #strand == 0 is + sense
                    chromosome_fragment = chromosome.read().replace('\n', '')[bingo_site - 1:bingo_site + DISTANCE - 1]  #limit = DISTANCE
                    #print "strand=", strand
                    #print chromosome_fragment
                else:
                    chromosome_fragment = chromosome.read().replace('\n', '')[bingo_site - DISTANCE - 1:bingo_site - 1]
                    #print "strand=", strand, "reversing..."
                    #print chromosome_fragment, '\n'
                    chromosome_fragment = reversecomp(chromosome_fragment)
                    #print chromosome_fragment, '\n'
                #print restre.search(chromosome_fragment)

    #fragment DNA by method of choice

                if FRAG == 'MB':
                    restre = re.compile('(TTAA|AGATCT)', re.IGNORECASE)  #define restriction site sequences, case-insensitive!
                    if  restre.search(chromosome_fragment):
                        #print restre.findall(chromosome_fragment)[0]
                        site = chromosome_fragment.find(restre.findall(chromosome_fragment)[0])
                        #print "cropped sequence=", str(site + 1), "bp"
                        if site > 14: #We want minimum 18 bp, including 2bp of RE site. This is for mapping purposes. Most fragments will be much larger.
                            #print chromosome_fragment[0:site + 2]   #will take first 2 bp of the restriction site.
                            entry1 = bingo_chromosome.replace('.fa', '') + '\t' + str(bingo_site) + '\t' + 'random_site_number_' + str(i) + '\t' + chromosome_fragment[0:site + 2] + '\n'
                    #random_site_list.write(entry1)
                            entry2 = '>' + bingo_chromosome.replace('.fa', '') + '_' + 'strand_' + str(strand) + '_' + str(bingo_site) + '_' + 'random_site_number_' + str(i) + '\n' + chromosome_fragment[0:site + 2] +'\n'
                            random_site_fasta.write(entry2)
                            i = i + 1
                            print 'generated site number  ', i

                if FRAG == 'NASB':
                    restre = re.compile('(GCTAGC|CCTAGG|ACTAGT|GGATCC)', re.IGNORECASE)  #pattern for searching the REsites, case-insensitive!
                    if  restre.search(chromosome_fragment):
                        #print restre.findall(chromosome_fragment)[0]
                        site = chromosome_fragment.find(restre.findall(chromosome_fragment)[0])
                        #print "cropped sequence=", str(site + 1), "bp"
                        if site > 14:  #We want minimum 18 bp, including 2bp of RE site. This is for mapping purposes. Most fragments will be much larger.
                            #print chromosome_fragment[0:site + 2]   #will take first 2 bp of the restriction site.
                            entry1 = bingo_chromosome.replace('.fa', '') + '\t' + str(bingo_site) + '\t' + 'random_site_number_' + str(i) + '\t' + chromosome_fragment[0:site + 2] + '\n'
                            #random_site_list.write(entry1)
                            entry2 = '>' + bingo_chromosome.replace('.fa', '') + '_' + 'strand_' + str(strand) + '_' + str(bingo_site) + '_' + 'random_site_number_' + str(i) + '\n' + chromosome_fragment[0:site + 2] +'\n'
                            random_site_fasta.write(entry2)
                            i = i + 1
                            print 'generated site number  ', i

                if FRAG == 'Random':
                    site = int(round(random.gauss(mu,sigma), 1))
                    #site = randint(rand_min, rand_max)
                    #print chromosome_fragment[0:site + 0]
                    #We didn't include a size limit here, as the model distribution strongly disfavors small fragments.
                    #This is borne out is the analysis of fragment lengths.
                    #We don't have that same control over the RE digested fragments.
                    entry1 = bingo_chromosome.replace('.fa', '') + '\t' + str(bingo_site) + '\t' + 'random_site_number_' + str(i) + '\t' + chromosome_fragment[0:site + 0] + '\n'
                    #random_site_list.write(entry1)
                    entry2 = '>' + bingo_chromosome.replace('.fa', '') + '_' + 'strand_' + str(strand) + '_' + str(bingo_site) + '_' + 'random_site_number_' + str(i) + '\n' + chromosome_fragment[0:site + 0] +'\n'
                    random_site_fasta.write(entry2)
                    i = i + 1
                    print 'generated site number  ', i

#random_site_list.close()
random_site_fasta.close()
