
from matplotlib.backends.backend_pdf import PdfPages

import pysam
import os
import pandas
import numpy
import argparse
import logging
import matplotlib.pyplot as plt

####################
##   ARGS INPUT   ##
####################

tool_description = """
The tool calculates the cross convolution between the forward and reverse strand.
For a typical ChIPseq experiment, the cross-convolution should have a distinct peak
for a shift equal to the fragment size. Another peak for the read-length might also
be observable. 
By default output is written to source file location.
Example usage:
cross_convolution_analysis.py reads.bam genome_table.tsv shift -o output.pdf
"""

# parse command line arguments
parser = argparse.ArgumentParser(description=tool_description,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
# positional arguments
parser.add_argument(
    "reads",
    help="Path to read file in bam format. Remark: the file needs to be sorted.")
parser.add_argument(
    'genome_table',
    help="Path to the table of the reference genome (e.g. hg18) "
         "listing the chromosomes and length (e.g. chr1 \t 247249719).")
parser.add_argument(
    'shift',
    type=int,
    help="Size of the shift. You will shift the reverse strand over the forward strand" +
         "in oder to determine the highest common read counts per position.")
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
logging.info("  reads: '{}'".format(args.reads))
logging.info("  genome_table: '{}'".format(args.genome_table))
logging.info("  shift: '{}'".format(args.shift))
if args.outfile:
    logging.info("  outfile: enabled writing to file")
    logging.info("  outfile: '{}'".format(args.outfile))
logging.info("")

if args.shift <= 0:
    raise Exception("[ERROR] Shift must be a positive integer bigger than 0.")

###############################
##   READ AND PROCESS DATA   ##
###############################

print('[START]')

# read in the two files provided by the user
genome_table = pandas.read_table(args.genome_table, sep='\t', names=['chrom', 'length']).set_index('chrom').T.to_dict(orient='list')

pure_bam = pysam.AlignmentFile(args.reads)

# position-wise vector for the forward and reverse strand
chr_vec_forward = numpy.zeros(1)
chr_vec_reverse = numpy.zeros(1)

# vector holding the positions that the read covers
pos_vec = []

buff_chrom = ''
error_check_chr_list = []

# array holding the overall result later used for the plot
conv_array = numpy.zeros(args.shift + 1)

# go over all reads
for read in pure_bam:

    # get the name of the chromosome
    chrom = pure_bam.get_reference_name(read.rname)

    if buff_chrom != chrom:

        if chrom in error_check_chr_list:
            raise Exception("[ERROR] bam file not sorted. Please sort the bam file.")

        if chrom not in genome_table:
            raise Exception("[ERROR] Chromosome could not be found in the genome table." +
                            "Please check your reference genome.")

        # if we have look at all reads for a chromosome then calculated the convolution
        # between forward and reverse strand
        if buff_chrom != '':
            print('[NOTE] finish')
            print('[NOTE] Convolve for ' + buff_chrom)
            conv_array += numpy.convolve(chr_vec_forward[::-1], chr_vec_reverse, mode='valid')
            print('[NOTE] finish')

        # init the two position-wise vectors (add additionally position for the forward strand
        # with the size of the shift in order to calculate the convolution only for the shifted positions)
        chr_vec_forward = numpy.zeros(genome_table[chrom][0] + args.shift)
        chr_vec_reverse = numpy.zeros(genome_table[chrom][0])
        buff_chrom = chrom
        error_check_chr_list.extend(chrom)
        print('[NOTE] Start reading for ' + chrom)

    # if the sam format bit-flag = 0 then this is the forward strand
    if read.flag == 0:

        # just in case, check if the read is too long, i.e. going over the chromosome
        # get all position that the read covers
        # (for the forward strand add the shift length which is necessary for the convolution)
        if read.pos+len(read.seq) > genome_table[chrom][0]:
            pos_vec = [int(i+args.shift) for i in range(read.pos, genome_table[chrom][0])]
        else:
            pos_vec = [int(i+args.shift) for i in range(read.pos, read.pos+len(read.seq))]

        # increment each position on the forward strand that the read covers
        chr_vec_forward[pos_vec] += 1

    # if the sam-format bit-flag = 16 then this is the reverse strand
    elif read.flag == 16:

        if read.pos + len(read.seq) > genome_table[chrom][0]:
            pos_vec = [i for i in range(read.pos, genome_table[chrom][0])]
        else:
            pos_vec = [i for i in range(read.pos, read.pos+len(read.seq))]

        chr_vec_reverse[pos_vec] += 1   # reverse

##############
##   PLOT   ##
##############

print('[NOTE] Create plot')

plotpath = os.path.dirname(os.path.abspath(__file__)) + '/'

pp = PdfPages(plotpath + 'Cross-Convolution.pdf')

plt.plot(conv_array, ls='--', marker='.', ms=10.0)
plt.xlabel('Shift Distance')
plt.ylabel('Convoluted Read Counts')
pp.savefig()
pp.close()

print('[FINISH]')
