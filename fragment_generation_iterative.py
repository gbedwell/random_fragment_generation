#!/usr/local/bin/python2

###########################################################################################################################
# Adapted from fragment_generation.py by Greg Bedwell (gregoryjbedwell@gmail.com), June 2020.
#
# This script is designed for iterative random fragment generation.
# It is written to output a new fragment file every N number of fragments. Each new file will have a new random seed (based on system time).
#
###########################################################################################################################

import os, string, random, re, sys
from string import maketrans, translate, join
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
print "Program writes out a joined fasta file (random_sites.fa) containing the generated fragments\n\n" #random_sites.txt currently commented out
print "-"*100
print "\n\n"
BUILD = raw_input("Enter the the genome build to use (hg19 or hg38)\n")
print "\n"

if (BUILD != 'hg19' and BUILD != 'hg38'):
    sys.exit('ERROR: Please enter a valid genome build identifier\n\n')

PATH = raw_input("Enter the path to the chromosome fasta files\n")
print "\n"
ITER = int(raw_input("Enter the desired number of iterations:\n"))
print "\n"
DIR = raw_input("Enter the name of the directory in which to place output files\n")
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

if not os.path.exists(DIR):
    os.makedirs(DIR)


for iter in range(1, ITER+1):
    path = PATH
    trans = maketrans("acgttumrwsykvhdbnACGTTUMRWSYKVHDBN", "TGCAAAKYWSRMBDHVNTGCAAAKYWSRMBDHVN")
    def reversecomp(sequence):
        sequence_list = list(sequence)
        sequence_list.reverse()
        reversed_sequence = join(sequence_list, '')
        reverse_complement_sequence = translate(reversed_sequence, trans)
        return reverse_complement_sequence
    random_site_fasta = open(DIR + '/random_hits_'+ str(iter)+ '.fa', 'w')
    random.seed()

    if BUILD == 'hg19':
        i = 0
        while i < N:
            current_site = random.randint(1,3137161264)
            counts_file = open(path + '/base_counts_hg19.txt', 'r')
            for line in counts_file.readlines():
                line_clean = line.replace('\n', '')
                foo = line_clean.split('\t')
                if (current_site > int(foo[2]) and current_site < int(foo[3])):
                    bingo_chromosome = foo[0]
                    bingo_site = current_site - int(foo[2])
                    chromosome = open(path + '/' + foo[0] + '.fa', 'r')
                    strand = random.randrange(2)
                    if strand == 0:
                        chromosome_fragment = chromosome.read().replace('\n', '')[bingo_site - 1:bingo_site + DISTANCE - 1]
                    else:
                        chromosome_fragment = chromosome.read().replace('\n', '')[bingo_site - DISTANCE - 1:bingo_site - 1]
                        chromosome_fragment = reversecomp(chromosome_fragment)


                    if FRAG == 'MB':
                        restre = re.compile('(TTAA|AGATCT)', re.IGNORECASE)  #define restriction site sequences, case-insensitive!
                        if  restre.search(chromosome_fragment):
                            #print restre.findall(chromosome_fragment)[0]
                            site = chromosome_fragment.find(restre.findall(chromosome_fragment)[0])
                            #print "cropped sequence=", str(site + 1), "bp"
                            if site > 14: #We want minimum 18 bp, including 2bp of RE site. This is for mapping purposes. Most fragments will be much larger. #Python counts from 0
                                #print chromosome_fragment[0:site + 2]   #will take first 2 bp of the restriction site.
                                entry1 = bingo_chromosome.replace('.fa', '') + '\t' + str(bingo_site) + '\t' + 'random_site_number_' + str(i) + '\t' + chromosome_fragment[0:site + 2] + '\n'
                        #random_site_list.write(entry1)
                                entry2 = '>' + bingo_chromosome.replace('.fa', '') + '_' + 'strand_' + str(strand) + '_' + str(bingo_site) + '_' + 'random_site_number_' + str(i) + '\n' + chromosome_fragment[0:site + 2] +'\n'
                                random_site_fasta.write(entry2)
                                i = i + 1
                                print 'generated site number  ', i

                    if FRAG == 'NASB':
                        restre = re.compile('(GCTAGC|CCTAGG|ACTAGT|GGATCC)', re.IGNORECASE)  #define restriction site sequences, case-insensitive!
                        if  restre.search(chromosome_fragment):
                            #print restre.findall(chromosome_fragment)[0]
                            site = chromosome_fragment.find(restre.findall(chromosome_fragment)[0])
                            #print "cropped sequence=", str(site + 1), "bp"
                            if site > 14: #We want minimum 18 bp, including 2bp of RE site. This is for mapping purposes. Most fragments will be much larger. #Python counts from 0
                                #print chromosome_fragment[0:site + 2]   #will take first 2 bp of the restriction site.
                                entry1 = bingo_chromosome.replace('.fa', '') + '\t' + str(bingo_site) + '\t' + 'random_site_number_' + str(i) + '\t' + chromosome_fragment[0:site + 2] + '\n'
                        #random_site_list.write(entry1)
                                entry2 = '>' + bingo_chromosome.replace('.fa', '') + '_' + 'strand_' + str(strand) + '_' + str(bingo_site) + '_' + 'random_site_number_' + str(i) + '\n' + chromosome_fragment[0:site + 2] +'\n'
                                random_site_fasta.write(entry2)
                                i = i + 1
                                print 'generated site number  ', i

                    if FRAG == 'Random':
                        site = int(round(random.gauss(mu,sigma), 1))
                        entry1 = bingo_chromosome.replace('.fa', '') + '\t' + str(bingo_site) + '\t' + 'random_site_number_' + str(i) + '\t' + chromosome_fragment[0:site + 0] + '\n'
                        entry2 = '>' + bingo_chromosome.replace('.fa', '') + '_' + 'strand_' + str(strand) + '_' + str(bingo_site) + '_' + 'random_site_number_' + str(i) + '\n' + chromosome_fragment[0:site + 0] +'\n'
                        random_site_fasta.write(entry2)
                        i = i + 1
                        print 'generated site number  ', i

    if BUILD == 'hg38':
        i = 0
        while i < N:
            current_site = random.randint(1,3209286105)
            counts_file = open(path + '/base_counts_hg38.txt', 'r')
            for line in counts_file.readlines():
                line_clean = line.replace('\n', '')
                foo = line_clean.split('\t')
                if (current_site > int(foo[2]) and current_site < int(foo[3])):
                    bingo_chromosome = foo[0]
                    bingo_site = current_site - int(foo[2])
                    chromosome = open(path + '/' + foo[0] + '.fa', 'r')
                    strand = random.randrange(2)
                    if strand == 0:
                        chromosome_fragment = chromosome.read().replace('\n', '')[bingo_site - 1:bingo_site + DISTANCE - 1]
                    else:
                        chromosome_fragment = chromosome.read().replace('\n', '')[bingo_site - DISTANCE - 1:bingo_site - 1]
                        chromosome_fragment = reversecomp(chromosome_fragment)


                    if FRAG == 'MB':
                        restre = re.compile('(TTAA|AGATCT)', re.IGNORECASE)  #define restriction site sequences, case-insensitive!
                        if  restre.search(chromosome_fragment):
                            #print restre.findall(chromosome_fragment)[0]
                            site = chromosome_fragment.find(restre.findall(chromosome_fragment)[0])
                            #print "cropped sequence=", str(site + 1), "bp"
                            if site > 14: #We want minimum 18 bp, including 2bp of RE site. This is for mapping purposes. Most fragments will be much larger. #Python counts from 0
                                #print chromosome_fragment[0:site + 2]   #will take first 2 bp of the restriction site.
                                entry1 = bingo_chromosome.replace('.fa', '') + '\t' + str(bingo_site) + '\t' + 'random_site_number_' + str(i) + '\t' + chromosome_fragment[0:site + 2] + '\n'
                        #random_site_list.write(entry1)
                                entry2 = '>' + bingo_chromosome.replace('.fa', '') + '_' + 'strand_' + str(strand) + '_' + str(bingo_site) + '_' + 'random_site_number_' + str(i) + '\n' + chromosome_fragment[0:site + 2] +'\n'
                                random_site_fasta.write(entry2)
                                i = i + 1
                                print 'generated site number  ', i

                    if FRAG == 'NASB':
                        restre = re.compile('(GCTAGC|CCTAGG|ACTAGT|GGATCC)', re.IGNORECASE)  #define restriction site sequences, case-insensitive!
                        if  restre.search(chromosome_fragment):
                            #print restre.findall(chromosome_fragment)[0]
                            site = chromosome_fragment.find(restre.findall(chromosome_fragment)[0])
                            #print "cropped sequence=", str(site + 1), "bp"
                            if site > 14: #We want minimum 18 bp, including 2bp of RE site. This is for mapping purposes. Most fragments will be much larger. #Python counts from 0
                                #print chromosome_fragment[0:site + 2]   #will take first 2 bp of the restriction site.
                                entry1 = bingo_chromosome.replace('.fa', '') + '\t' + str(bingo_site) + '\t' + 'random_site_number_' + str(i) + '\t' + chromosome_fragment[0:site + 2] + '\n'
                        #random_site_list.write(entry1)
                                entry2 = '>' + bingo_chromosome.replace('.fa', '') + '_' + 'strand_' + str(strand) + '_' + str(bingo_site) + '_' + 'random_site_number_' + str(i) + '\n' + chromosome_fragment[0:site + 2] +'\n'
                                random_site_fasta.write(entry2)
                                i = i + 1
                                print 'generated site number  ', i

                    if FRAG == 'Random':
                        site = int(round(random.gauss(mu,sigma), 1))
                        entry1 = bingo_chromosome.replace('.fa', '') + '\t' + str(bingo_site) + '\t' + 'random_site_number_' + str(i) + '\t' + chromosome_fragment[0:site + 0] + '\n'
                        entry2 = '>' + bingo_chromosome.replace('.fa', '') + '_' + 'strand_' + str(strand) + '_' + str(bingo_site) + '_' + 'random_site_number_' + str(i) + '\n' + chromosome_fragment[0:site + 0] +'\n'
                        random_site_fasta.write(entry2)
                        i = i + 1
                        print 'generated site number  ', i


random_site_fasta.close()
