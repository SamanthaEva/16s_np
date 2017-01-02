#!/usr/bin/env python
import sys
import os
import pyfasta
"""
Takes a fastafile and changes the fasta header (i.e. ">scaffold 0 3 34342 311")
to >$index_NN.
Author: Moa Hammarstrom
""" 

# Make sure the user entered the command line arguments correctly
if len(sys.argv) != 4:
    sys.stderr.write("USAGE IS: %s <fasta input> <index> <fasta output>\n" % sys.argv[0])
    sys.exit(1)

# Parse arguments
fileIn = sys.argv[1]
index = sys.argv[2]
fileOut = sys.argv[3]

fasta_in = pyfasta.Fasta(fileIn)
fasta_out = open(fileOut, "wb")

seq_names = fasta_in.keys()
nr_seqs = len(seq_names)
counter = 1
for seq in seq_names:
    fasta_out.write(">" + str(index) + "_" + str(counter) + "\n")
    fasta_out.write(str(fasta_in[seq]))
    if counter < nr_seqs:
        fasta_out.write("\n")
        counter+=1

fasta_out.close()
