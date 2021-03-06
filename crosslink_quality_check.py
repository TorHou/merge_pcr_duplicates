from matplotlib.backends.backend_pdf import PdfPages
from scipy import stats

import os
import matplotlib.pyplot as plt
import numpy
import pandas
import logging
import argparse

####################
##   ARGS INPUT   ##
####################

tool_description = """
This tool checks the reproducibility of the crosslinking sites between samples in fasta format.
Ideally the trend should follow a diagonal line with the highest reproducible motif in the right
upper most corner.
By default output is written to source file location.
Example usage:
crosslink_quality_check.py file1.fasta file2.fasta kmer_length -o output.pdf
"""

# parse command line arguments
parser = argparse.ArgumentParser(description=tool_description,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
# positional arguments
parser.add_argument(
    "file1",
    help="Path to first fasta file.")
parser.add_argument(
    'file2',
    nargs='+',
    help="Path to the other fasta files. Specify at least one more file for file1 to be compared to.")
parser.add_argument(
    'kmer_length',
    type=int,
    help="Length of the kmers. Keep in mind that the sequences should be long enough.")
# optional arguments
parser.add_argument(
    "-o", "--outfile",
    help="Write results to this file.")
parser.add_argument(
    "-d", "--debug",
    help="Print lots of debugging information",
    action="store_true")

args = parser.parse_args()
if args.debug:
    logging.basicConfig(level=logging.DEBUG, format="%(asctime)s - %(filename)s - %(levelname)s - %(message)s")
else:
    logging.basicConfig(format="%(filename)s - %(levelname)s - %(message)s")
logging.info("Parsed arguments:")
logging.info("  file1: '{}'".format(args.file1))
logging.info("  file2: '{}'".format(args.file2))
logging.info("  kmer_length: '{}'".format(args.kmer_length))
if args.outfile:
    logging.info("  outfile: enabled writing to file")
    logging.info("  outfile: '{}'".format(args.outfile))
logging.info("")

if args.kmer_length <= 1:
    raise Exception("[ERROR] kmer length too short. Your kmer length makes no sense.")

###################
##   READ DATA   ##
###################

print("[START]")
print("[NOTE] Read data")

# your read data in bed6 format
files = [args.file1]
files.extend(args.file2)

# read first line of first file to get length of the sequences
tmp = open(files[0])
firstline = tmp.readline()

# we assume that the middle of the sequence is the crosslink nucleotide
cl_nucleotide = int(len(firstline)/2)

# get starting and end point for the kmers
start_iter = cl_nucleotide - args.kmer_length + 1
end_iter = cl_nucleotide + args.kmer_length + 1

tmp.close()

print("[NOTE] finish")

######################
##   PROCESS DATA   ##
######################

print("[NOTE] Process data")

# create a dictornary for the kmers
kmer_dict = dict()

f = 0
for file in files:
    with open(file) as openfileobject:
        for seq in openfileobject:
            # check length of kmer
            if len(seq) < (args.kmer_length*2 - 1):
                raise Exception("[ERROR] kmer length too long, in other words, sequence is too short.")

            # go over your sequence and generate kmers of length args.kmer_length
            # put kmers into dictonary
            for i in range(start_iter, cl_nucleotide+1):
                kmer = seq[i:i+args.kmer_length]
                # if kmer already exists, then increment for file f
                if kmer in kmer_dict:
                    kmer_dict[kmer][f] += 1
                # if kmer does not exist, then create a new vector of size len(files)
                else:
                    init = numpy.zeros(len(files))
                    init[f] = 1
                    kmer_dict[kmer] = init

    openfileobject.close()
    f +=1

sum_vector = numpy.zeros(len(files))

kmer_values_vectors = kmer_dict.values()

# get the total number of kmers for each files (sample)
for values_vector in kmer_values_vectors:
    sum_vector += values_vector

# calculate the realtive abundance of each kmer for each file (sample)
for kmer in kmer_dict:
    kmer_dict[kmer] = kmer_dict[kmer]/sum_vector

print("[NOTE] finish")

##############
##   PLOT   ##
##############

print("[NOTE] Make plots")

plotpath = os.path.dirname(os.path.abspath(__file__)) + '/'

# create a plot and list of motifs with their sorted relative abundance for each pair of files
for i in range(0,len(files)):
    for j in range(i+1, len(files)):

        outfile_name = ""
        if args.outfile:
            outfile_name = args.outfile
        else:
            outfile_name = plotpath + 'Crosslink_Kmer_Quality_Check_' + str(i) + '_' + str(j) + '.pdf'

        pp = PdfPages(outfile_name)

        df = pandas.DataFrame(kmer_dict).T

        # do linear regression for the two files
        slope, intercept, r_value, p_value, std_err = stats.linregress(df[i], df[j])

        plt.plot(df[i], df[j],  ls='', marker='.', ms=10.0)
        plt.ylabel(files[i])
        plt.xlabel(files[j])

        max_x = max(df[i])
        max_y = max(df[j])

        plt.text(max_x/3, max_y/2 + max_y*.2, "R" + r'$^2 =$' + " " + str(r_value), fontsize=15)

        pp.savefig()
        pp.close()

        outfile_name = 'test-data/reproducible_motifs_' + str(i) + '_' + str(j) + '.csv'
        df_sorted = df[[i,j]]
        df_sorted = df_sorted.sort_values([i, j], ascending=[False, False])
        df_sorted.columns = [files[i], files[j]]
        df_sorted.to_csv(outfile_name, sep='\t')

print("[FINISH]")

