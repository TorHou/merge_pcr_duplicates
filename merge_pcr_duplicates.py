import pysam
from collections import Counter
import argparse
from Bio import SeqIO

bam_data = []
score = Counter()
endResult = []
fastq_data = {}
bam_filter = []
merge_data = []



def get_bam_filter(line, chrom_bam):
    if line.is_reverse:
        strand = '-'
    else:
        strand = '+'
    if 'N' not in fastq_data[line.query_name]:
        bam_filter.append((line.reference_start,
                           line.reference_start + line.reference_length,
                           line.query_name, strand,
                           chrom_bam, fastq_data[line.query_name]))
    return bam_filter


def bam_reader(bam_file, xs_tag):
    global merge_data
    bam = pysam.AlignmentFile(bam_file, "rb")
    filtered_data = []
    chromosomes = bam.references
    for chrom in chromosomes:
        data = bam.fetch(chrom, multiple_iterators=True)

        for entry in data:
            try:
                chrom_bam = bam.get_reference_name(entry.reference_id)
            except:
                continue
            if (entry.is_unmapped == False and chrom_bam == chrom):
                if xs_tag and entry.has_tag("XS"):
                    continue
                else:
                    filtered_data = get_bam_filter(entry, chrom_bam)
        chromosome_info(filtered_data)
        print_results(merge_data)
        merge_data = []
    bam.close()

def fastq_reader(fastq_file):
    fastq_dt = {}
    for fastq_data in SeqIO.parse(fastq_file, "fastq"):
        fastq_dt[str(fastq_data.id)] = str(fastq_data.seq)
    return fastq_dt


def chromosome_counter(i, bam_data):
    if args.ends:
        merge_checks = (bam_data[i][4], bam_data[i][0], bam_data[i][1], bam_data[i][5], bam_data[i][3])
    else:
        merge_checks = (bam_data[i][4], bam_data[i][0], bam_data[i][5], bam_data[i][3])
    score[merge_checks] += 1
    if score[merge_checks] == 1:
        merge_data.append((bam_data[i][4], bam_data[i][0], bam_data[i][1], bam_data[i][2], bam_data[i][3]))
    
    #append_merge = merge_data.append
    #if score[merge_checks] == 1:
        #append_merge((bam_data[i][4], bam_data[i][0], bam_data[i][1], bam_data[i][2], bam_data[i][3]))

def chromosome_info(bam_data):
    chr_info = []
    for i in range(len(bam_data)):
        rec_id = bam_data[i][2]
        if rec_id in fastq_data:
            chromosome_counter(i, bam_data)
        else:
            print(rec_id + " ID not found in fastq file")

def print_results(endResult):
    if endResult != [] :
        with open(args.output_file, "w") as f:
            if args.ends:
                for entry in endResult:
                    wr = ("%s\t%s\t%s\t%s\t%s\t%s" % (entry[0], entry[1], entry[2], entry[3],
                                                      score[(entry[0], entry[1], entry[2], fastq_data[entry[3]],
                                                             entry[4])], entry[4]))
                    f.write(wr + '\n')
            else:
                for entry in endResult:
                    wr = ("%s\t%s\t%s\t%s\t%s\t%s" % (entry[0], entry[1], entry[2], entry[3],
                                                      score[(entry[0], entry[1], fastq_data[entry[3]],
                                                             entry[4])], entry[4]))
                    f.write(wr + '\n')

tool_description = """
Merge PCR duplicates according to random barcode library.
Input:
* bam file containing fastq read-id and other details
* fastq library of random barcodes
Output:
* bed6 file with random barcode in name field and number of PCR duplicates as
  score, sorted by fields chrom, start, stop, strand, name
Example usage:
python3 pysam_test.py example1.bam example2.fa
"""

epilog = """
Author: Fayyaz Hussain
Status: Testing
"""


# parse command line arguments
parser = argparse.ArgumentParser(description=tool_description,
                                 epilog=epilog, formatter_class=argparse.RawDescriptionHelpFormatter)
parser.add_argument("bam_file", help="Path to bam file containing alignments.", metavar='BAM_File')
parser.add_argument("fastq_file", help="Path to fastq barcode library.", metavar='FASTQ_File')
parser.add_argument("-o", "--output_file", required=True, help="Write results to this file.",
                    metavar='Output_File')
parser.add_argument("--filter_by_nh_tag", action='store_true', default=False, help="Mapped reads with XS tag will be excluded.")
parser.add_argument("-e", "--ends", action='store_true', default=False, help="Consider sequence end coordinates when merging the PCR duplicates. By default reads with the same chromosome, strand, start and barcode are merged")
args = parser.parse_args()
args.filter_by_nh_tag
fastq_data = fastq_reader(args.fastq_file)
bam_reader(args.bam_file, args.filter_by_nh_tag)

