#!/usr/bin/python

#################################################################################################################################################################################
####################### Authorship information: The basis of this script was written by Peter Cherepanov (peter.cherepanov@crick.ac.uk) 2006, 2014, 2021. #######################
######################################### It was modified into its current from by Greg Bedwell (gregoryjbedwell@gmail.com) 2020, 2021. #########################################
################################################### Please feel free to contact Greg Bedwell with any questions or concerns. ####################################################
#################################################################################################################################################################################

import os, re, sys
from string import maketrans, translate, join
import itertools
import random
import argparse


parser = argparse.ArgumentParser(description= "This script generates a defined number of random DNA fragments from the human genome. It is written to accomidate either pattern-based fragmentation (e.g. via restriction endonucleases) or random fragmentation. For the script to run properly, the user must provide a 'base_text.txt' file. This file should be a 4-column tsv file where column 1) is the names of all chromosomes present in the genome build, 2) is the absolute size of each chromosome (in bp), 3) is the relative start position of each chromosome relative to the total genome size, and 4) is the relative end position of each chromosome relative to the total genome size.")

parser.add_argument("-build", "--build", help="The genome build being used. Currently, hg19, hg38, and CHM13_v1.1 are implemented as options.")
parser.add_argument("-path", "--path", help="The path to the chromosome fasta files and the base_counts file.")
parser.add_argument("-bc", "--bc_filename", help="The name of the base_counts.txt file.")
parser.add_argument("-N", "--N", help="The number of random fragments to generate.")
parser.add_argument("-frag", "--frag", help="The fragmentation method. Currently accepts M for MseI, MB for MseI/BglII digestion, NASB for NheI/AvrII/SpeI/BamHI digestion, or random for random fragmentation.")
parser.add_argument("-dir", "--dir", help="The directory to save the output to.")
parser.add_argument("-distance", "--distance", help="The maximum distance allowed to the cleavage site.")
parser.add_argument("-mu", "--mu", help="Only used with random fragmentation. Defines the median fragment size.")
parser.add_argument("-sigma", "--sigma", help="Only used with random fragmentation. Defines the standard deviation of the median fragment size.")
parser.add_argument("-file_num", "--file_num", help="Numbers the output files. Useful for generating multiple random fragment files in parallel.")
parser.add_argument("-output", "--output", help="The output filetype. Can be fasta or fastq.")


args = parser.parse_args()

if args.build is None or args.path is None or args.bc_filename is None or args.N is None or args.frag is None or args.distance is None or args.dir is None or args.file_num is None or args.output is None:
    sys.exit("\nYou must provide values for 1) build, 2) path, 3) bc, 4) N, 5) frag, 6) distance, 7) dir, 8) file_num, and 9) output to run this script!\nOne or more of these values is currently not defined!\n\nSee help (-h or --help) for variable definitions.\n")
else:
    print("\n")
    print('build = ' + args.build)
    print('path = ' + args.path)
    print('bc_filename = ' + args.bc_filename)
    print('N = ' + args.N)
    print('frag = ' + args.frag)
    print('distance = ' + args.distance)
    print('directory = ' + args.dir)
    print('output = ' + args.output)

    if args.frag == str('random'):
        if args.mu and args.sigma:
            print('mu = ' + args.mu)
            print('sigma = ' + args.sigma)
            print("\n")
        else:
            sys.exit("\nYou must provide values for mu and sigma when frag = random!\n")


build, frag = args.build, args.frag

if (build != 'hg19' and build != 'hg38' and build != 'CHM13_v1.1'):
    sys.exit('\nERROR: Please enter a valid genome build identifier!\nOptions are: hg19, hg38, CHM13_v1.1\n')


if (frag != 'M' and frag != 'MB' and frag != 'NASB' and frag != 'random' ):
    sys.exit('\nERROR: Please enter a valid genome fragmentation method\n')


if not os.path.exists(args.dir):
    os.makedirs(args.dir)



file_str1, file_str2 = "_".join([str(args.build), str(args.frag), 'R1', str(args.file_num)]), "_".join([str(args.build), str(args.frag), 'R2', str(args.file_num)])

trans = maketrans("acgttumrwsykvhdbnACGTTUMRWSYKVHDBN", "TGCAAAKYWSRMBDHVNTGCAAAKYWSRMBDHVN")  #reverse complement rules

def reversecomp(sequence):
    sequence_list = list(sequence)
    sequence_list.reverse()
    reversed_sequence = join(sequence_list, '')
    reverse_complement_sequence = translate(reversed_sequence, trans)
    return reverse_complement_sequence

if args.output == str('fasta'):

    random_site_fasta1 = open(str(args.dir) + '/random_fragments_' + str(file_str1) + '.fa', 'w')
    random_site_fasta2 = open(str(args.dir) + '/random_fragments_' + str(file_str2) + '.fa', 'w')

if args.output == str('fastq'):

    random_site_fastq1 = open(str(args.dir) + '/random_fragments_' + str(file_str1) + '.fq', 'w')
    random_site_fastq2 = open(str(args.dir) + '/random_fragments_' + str(file_str2) + '.fq', 'w')

random.seed()    #set the seed based on the current system time

i = 0
while i < int(args.N):
    if args.build == str('hg19'):
        current_site = random.randint(1,3137161264)

    if args.build == str('hg38'):
        current_site = random.randint(1,3209286105)

    if args.build == str('CHM13_v1.1'):
        current_site = random.randint(1,3056899953)

    counts_file = open(str(args.path) + '/' + str(args.bc_filename), 'r') #base_counts file should be named
    for line in counts_file.readlines():
        line_clean = line.replace('\n', '')
        foo = line_clean.split('\t')

        if (current_site > int(foo[2]) and current_site < int(foo[3])): #if current site > chr start (foo[2]) and < chr end (foo[3])
            bingo_chromosome = foo[0] #define chromosome containing integration site
            bingo_site = current_site - int(foo[2])
                        #print bingo_chromosome, bingo_site, '\n'
                        #chromosome = open(path + '/' + foo[0] + '.fa', 'r') #open chromosome file containing bingo_site
                        #print chromosome.readline()
            with open(str(args.path) + '/' + foo[0] + '.fa', 'r') as f:
                strand = random.randrange(2)
                if strand == 0:
                    chromosome_fragment = "".join(line.strip() for line in itertools.islice(f, 1, None))[bingo_site - 1:bingo_site + int(args.distance) - 1]
                else:
                    chromosome_fragment = "".join(line.strip() for line in itertools.islice(f, 1, None))[bingo_site - int(args.distance) - 1:bingo_site - 1]

            #fragment DNA by method of choice

                if args.frag == str('M'):
                    restre = re.compile('(TTAA)', re.IGNORECASE)
                if args.frag == str('MB'):
                    restre = re.compile('(TTAA|AGATCT)', re.IGNORECASE)
                if args.frag == str('NASB'):
                    restre = re.compile('(GCTAGC|CCTAGG|ACTAGT|GGATCC)', re.IGNORECASE)

                if args.frag == str('random'):
                    restre = None

                if restre is None:

                    site = int(round(random.gauss(int(args.mu),int(args.sigma)), 1))

                    if site > 19 and site < 899: # We want fragments between 20 bp and 900 bp

                        seq1 = chromosome_fragment[0:min(site+0,150)]
                        seq2 = chromosome_fragment[max(site+0-150,0):site+0]
                        seq2 = reversecomp(seq2)
                            #site = randint(rand_min, rand_max)
                            #print chromosome_fragment[0:site + 0]
                            #We didn't include a size limit here, as the model distribution strongly disfavors small fragments.
                            #This is borne out is the analysis of fragment lengths.
                            #We don't have that same control over the RE digested fragments.
                            #entry1 = bingo_chromosome.replace('.fa', '') + '\t' + str(bingo_site) + '\t' + 'random_site_number_' + str(i) + '\t' + chromosome_fragment[0:site + 0] + '\n'
                            #random_site_list.write(entry1)

                        entry1_name = "_".join([bingo_chromosome.replace('.fa', ''), 'strand', str(strand), str(bingo_site), 'random_site_number', str(i+1)])

                        entry2_name = "_".join([bingo_chromosome.replace('.fa', ''), 'strand', str(strand), str(bingo_site), 'random_site_number', str(i+1)])

                        if args.output == str('fasta'):

                            entry1 = '>' + bingo_chromosome.replace('.fa', '') + '_' + 'strand_' + str(strand) + '_' + str(bingo_site) + '_' + 'random_site_number_' + str(i+1) + '\n' + seq1 + '\n'
                            random_site_fasta1.write(entry1)

                            entry2 = '>' + bingo_chromosome.replace('.fa', '') + '_' + 'strand_' + str(strand) + '_' + str(bingo_site) + '_' + 'random_site_number_' + str(i+1) + '\n' + seq2 + '\n'
                            random_site_fasta2.write(entry2)

                            i = i + 1

                            print 'generated site number  ', i

                        if args.output == str('fastq'):

                            Q1 = 'J'
                            Q2 = ['E','F','G']
                            Q3 = ['A','B']

                            samp1 = ''.join(random.sample(Q1*85,85))

                            tmp2 = random.sample(Q1*45,45)+[random.choice(Q2) for _ in range(10)]
                            samp2_1 = ''.join(random.sample(tmp2, len(tmp2)))
                            tmp2 = random.sample(Q1*45,45)+[random.choice(Q2) for _ in range(10)]
                            samp2_2 = ''.join(random.sample(tmp2, len(tmp2)))

                            tmp3 = random.sample(Q1*5,5)+[random.choice(Q3) for _ in range(5)]
                            samp3_1 = ''.join(random.sample(Q1*5,5)+[random.choice(Q3) for _ in range(5)])
                            tmp3 = random.sample(Q1*5,5)+[random.choice(Q3) for _ in range(5)]
                            samp3_2 = ''.join(random.sample(Q1*5,5)+[random.choice(Q3) for _ in range(5)])

                            fakeq1 = ''.join([samp1, samp2_1, samp3_1])
                            fakeq2 = ''.join([samp1, samp2_2, samp3_2])

                            entry1 = '@' + bingo_chromosome.replace('.fa', '') + '_' + 'strand_' + str(strand) + '_' + str(bingo_site) + '_random_site_number_' + str(i+1) + '\n' + seq1 +'\n' +'+\n'+fakeq1[0:len(seq1)]+'\n'
                            random_site_fastq1.write(entry1)

                            entry2 = '@' + bingo_chromosome.replace('.fa', '') + '_' + 'strand_' + str(strand) + '_' + str(bingo_site) + '_random_site_number_' + str(i+1) + '\n' + seq2 +'\n' +'+\n'+fakeq2[0:len(seq2)]+'\n'
                            random_site_fastq2.write(entry2)

                            i = i + 1

                            print 'generated site number  ', i

                else:
                    if restre.search(chromosome_fragment):

                        site = chromosome_fragment.find(restre.findall(chromosome_fragment)[0])

                        if site > 18 and site < 899: # Gives a minimum fragment length of 19 bp + 1 bp of the restriction site added in the following lines

                            seq1 = chromosome_fragment[0:min(site+1,150)]
                            seq2 = chromosome_fragment[max(site+1-150,0):site+1]
                            seq2 = reversecomp(seq2)

                            if args.output == str('fasta'):

                                entry1 = '>' + bingo_chromosome.replace('.fa', '') + '_' + 'strand_' + str(strand) + '_' + str(bingo_site) + '_random_site_number_' + str(i+1) + '\n' + seq1 + '\n'
                                random_site_fasta1.write(entry1)

                                entry2 = '>' + bingo_chromosome.replace('.fa', '') + '_' + 'strand_' + str(strand) + '_' + str(bingo_site) + '_random_site_number_' + str(i+1) + '\n' + seq2 + '\n'
                                random_site_fasta2.write(entry2)

                                i = i + 1

                                print 'generated site number  ', i



                            if args.output == str('fastq'):

                                Q1 = 'J'
                                Q2 = ['E','F','G']
                                Q3 = ['A','B']

                                samp1 = ''.join(random.sample(Q1*85,85))

                                tmp2 = random.sample(Q1*45,45)+[random.choice(Q2) for _ in range(10)]
                                samp2_1 = ''.join(random.sample(tmp2, len(tmp2)))
                                tmp2 = random.sample(Q1*45,45)+[random.choice(Q2) for _ in range(10)]
                                samp2_2 = ''.join(random.sample(tmp2, len(tmp2)))

                                tmp3 = random.sample(Q1*5,5)+[random.choice(Q3) for _ in range(5)]
                                samp3_1 = ''.join(random.sample(Q1*5,5)+[random.choice(Q3) for _ in range(5)])
                                tmp3 = random.sample(Q1*5,5)+[random.choice(Q3) for _ in range(5)]
                                samp3_2 = ''.join(random.sample(Q1*5,5)+[random.choice(Q3) for _ in range(5)])

                                fakeq1 = ''.join([samp1, samp2_1, samp3_1])
                                fakeq2 = ''.join([samp1, samp2_2, samp3_2])

                                entry1 = '@' + bingo_chromosome.replace('.fa', '') + '_' + 'strand_' + str(strand) + '_' + str(bingo_site) + '_random_site_number_' + str(i+1) + '\n' + seq1 +'\n' +'+\n'+fakeq1[0:len(seq1)]+'\n'
                                random_site_fastq1.write(entry1)

                                entry2 = '@' + bingo_chromosome.replace('.fa', '') + '_' + 'strand_' + str(strand) + '_' + str(bingo_site) + '_random_site_number_' + str(i+1) + '\n' + seq2 +'\n' +'+\n'+fakeq2[0:len(seq2)]+'\n'
                                random_site_fastq2.write(entry2)

                                i = i + 1

                                print 'generated site number  ', i

if args.output == str('fasta'):
    random_site_fasta1.close()
    random_site_fasta2.close()

if args.output == str('fastq'):
    random_site_fastq1.close()
    random_site_fastq2.close()
