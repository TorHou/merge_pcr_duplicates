from __future__ import division
from collections import Counter
from matplotlib.backends.backend_pdf import PdfPages

import pandas
import matplotlib.pyplot as plt
import os
import argparse
import logging

####################
##   ARGS INPUT   ##
####################

tool_description = """
Creates a plot for quality control which depicts the over all reproducibility between
different samples. It counts the total number of cDNA that covers a particular position
between all samples and finds the positions that are shared between samples. The amount
of cDNA that covers a reproduced position (y-axis) is plotted against the total amount of
cDNA (x-axis).
By default output is written to source file location.
Example usage:
chromosome_quality_check.py file1.bed file2.bed --out output.pdf
"""

# parse command line arguments
parser = argparse.ArgumentParser(description=tool_description,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
# positional arguments
parser.add_argument(
    "file1",
    help="Path to first bed file.")
parser.add_argument(
    'file2',
    nargs='+',
    help="Path to the other bed files. Specify at least one more file for file1 to be compared to.")
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
if args.outfile:
    logging.info("  outfile: enabled writing to file")
    logging.info("  outfile: '{}'".format(args.outfile))
logging.info("")

###################
##   READ DATA   ##
###################

print("[START]")
print("[NOTE] Read data")

# your read data in bed6 format
files = [args.file1]
files.extend(args.file2)

groupcount = {}
main_groupcount = []
in_both = [0 for i in range(5)]
all = [0 for i in range(5)]

# go over all files
# group reads by chromosome start and strand (i.e., see with positions are covered)
for filename in files:
    df = pandas.read_table(filename, sep='\t', names=['chrom', 'start', 'stop', 'name', 'score', 'strand'])
    grouped = df.groupby(['chrom', 'start', 'strand'])
    groupcount = grouped.size().to_dict()
    main_groupcount.append(groupcount)

print("[NOTE] finish")

######################
##   PROCESS DATA   ##
######################

def count_values_for_set(s):
    "This counts the occurrence of the values in a set"
    sorted_set = sorted(s.values())
    return Counter(sorted_set);

print("[NOTE] Process data")

common_sets = set(main_groupcount[0].keys())

# reads that are the same between the sets
for i in range(0,len(main_groupcount)):
    sample = main_groupcount[i]
    common_sets = common_sets & set(sample.keys())
    counted_set = count_values_for_set(sample)
    # add to the counts the number of regions of each sample that are x-times duplicated
    for x in counted_set.keys():
        all[min(x - 1, 4)] += counted_set[x]

# add to the intersection-counts for each sample the number of reads
# that are x-times duplicated
for region in common_sets:
    for sample in main_groupcount:
        in_both[min(sample[region] - 1, 4)] += 1

print("In both files: ", in_both)
print("All :", all)

Ratio_both = ([(x / y) for x, y in zip(in_both, all)])

print("[NOTE] finish")

##############
##   PLOT   ##
##############

print("[NOTE] Make plots")

plotpath = os.path.dirname(os.path.abspath(__file__)) + '/'

pp = PdfPages(plotpath + 'Chromosome_Quality_Check.pdf')

col_label = ['container_1', 'container_2', 'container_3', 'container_4', 'container_5+']
chart = pandas.DataFrame(dict(zip(col_label, [Ratio_both])))
plt.plot(chart, ls='--', marker='.', ms=10.0)
plt.ylabel('Portion of cDNA for Reproduced Positions (%)')
plt.xlabel('cDNA Count per Position')
plt.axis([0, 4, 0, 1])
plt.xticks([0,1,2,3,4], ['1','2','3','4','>=5'])
pp.savefig()
pp.close()

print("[FINISH]")