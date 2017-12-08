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
crosslink_quality_check.py exp_rep_1.fasta exp_rep_2.fasta controls.fasta kmer_length -o output.pdf
"""

# parse command line arguments
parser = argparse.ArgumentParser(description=tool_description,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
# positional arguments
parser.add_argument(
    "exp_rep_1",
    help="Path to first experiment replicate fasta file.")
parser.add_argument(
    "exp_rep_2",
    help="Path to second experiment replicate fasta file.")
parser.add_argument(
    'controls',
    nargs='+',
    help="Path to the control fasta files. Specify at least one more files.")
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
logging.info("  exp_rep_1: '{}'".format(args.exp_rep_1))
logging.info("  exp_rep_2: '{}'".format(args.exp_rep_2))
logging.info("  controls: '{}'".format(args.controls))
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
files = [args.exp_rep_1, args.exp_rep_2]
files.extend(args.controls)

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

df = pandas.DataFrame(kmer_dict).T

# sort the dictionary (now dataframe) accoridng two the first two files (here replicates of the experiemtn)
df_sorted = df.sort_values([0, 1], ascending=[False, False])
df_sorted.columns = files
df_sorted.to_csv(plotpath + 'reproducible_motifs.csv', sep='\t')

# change colun names back to integer for convience
df_sorted.columns = [x for x in range(len(files))]

# Find the n most reproducible motifs in the two replicates of your experiment
n = 10
top_n_motifs = ["red" for x in range(n)]
rest_of_points = ["blue" for x in range(len(df_sorted[0]) - n)]
colors_for_scatterplot = top_n_motifs + rest_of_points

# create a plot and list of motifs with their sorted relative abundance for each pair of files
p = 1
for i in range(0,len(files)):
    for j in range(i+1, len(files)):

        outfile_name_plot = ""
        outfile_name_motif_table = ""
        if args.outfile:
            outfile_name_plot = args.outfile + '_' + str(p) + '.pdf'
        else:
            outfile_name_plot = plotpath + 'Crosslink_Kmer_Quality_Check_' + str(i) + '_' + str(j) + '_' + str(p) + '.pdf'

        p  += 1
        pp = PdfPages(outfile_name_plot)

        # do linear regression for the two files
        slope, intercept, r_value, p_value, std_err = stats.linregress(df_sorted[i], df_sorted[j])

        plt.scatter(df_sorted[i], df_sorted[j], c=colors_for_scatterplot, s=2)
        plt.ylabel(files[i])
        plt.xlabel(files[j])

        max_x = max(df_sorted[i])
        max_y = max(df_sorted[j])

        plt.title("R" + r'$^2 =$' + " " + str(r_value), fontsize=15)

        pp.savefig()
        pp.close()
        plt.close()

print("[FINISH]")

