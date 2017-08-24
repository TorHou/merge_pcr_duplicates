import pysam
from collections import Counter
import argparse
from Bio import SeqIO

bam_data = []
score = Counter()
endResult = []
merge_data = []
fastq_data = {}
bam_filter = []


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

def bam_reader(bam_file, end, tag):
    bam = pysam.AlignmentFile(bam_file, "rb")
    filtered_data = []
    chromosomes = bam.references
    for chrom_all in chromosomes:
        data = bam.fetch(multiple_iterators=True, until_eof=True)

        for line in data:
            try:
                chrom_bam = bam.get_reference_name(line.reference_id)
            except:
                continue
            if tag == False:
                if (line.is_unmapped == False and line.has_tag("XS") == False and chrom_bam == chrom_all):
                    filtered_data = get_bam_filter(line, chrom_bam)

            elif tag == True:
                if (line.is_unmapped == False and chrom_bam == chrom_all):
                    filtered_data = get_bam_filter(line, chrom_bam)

        bam_data = chromosome_info(filtered_data, end)
        printing(bam_data, end)
    bam.close()


def fastq_reader(fastq_file):
    fastq_dt = {}
    update_fastq = fastq_dt.update
    for fastq_data in SeqIO.parse(fastq_file, "fastq"):
        update_fastq({str(fastq_data.id): str(fastq_data.seq)})
    return fastq_dt


def chromosome_counter(i, end, bam_data):
    if end == True:
        merge_checks = (bam_data[i][4], bam_data[i][0], bam_data[i][1], bam_data[i][5], bam_data[i][3])
    else:
        merge_checks = (bam_data[i][4], bam_data[i][0], bam_data[i][5], bam_data[i][3])
    score.update([merge_checks])
    append_merge = merge_data.append
    if score[merge_checks] == 1:
        append_merge((bam_data[i][4], bam_data[i][0], bam_data[i][1], bam_data[i][2], bam_data[i][3]))

    return merge_data

def chromosome_info(bam_data, end):
    chr_info = []
    for i in range(len(bam_data)):
        rec_id = bam_data[i][2]
        if rec_id in fastq_data:
            chr_info = chromosome_counter(i, end, bam_data)
        else:
            print(rec_id + " ID not found in fastq file")

    return chr_info

def printing(endResult, end):
    if endResult != [] :
        with open(args.output_file, "w") as f:
            if end == True:
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
Barcodes containing uncalled base 'N' are removed.
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
parser = argparse.ArgumentParser(prog="Chromosomes Information", description=tool_description,
                                 epilog=epilog, formatter_class=argparse.RawDescriptionHelpFormatter)
parser.add_argument("bam_file", help="Path to bam file containing alignments.", metavar='BAM_File')
parser.add_argument("fastq_file", help="Path to fastq barcode library.", metavar='FASTQ_File')
parser.add_argument("-o", "--output_file", required=True, help="Write results to this file.",
                    metavar='Output_File')
parser.add_argument("-t", "--tag", action='store_true', default=False, help="If XS tag is to be excluded.")
parser.add_argument("-e", "--end", action='store_false', default=True, help="If sequence end needs to be considered.")
args = parser.parse_args()
tag = args.tag
end = args.end
fastq_data = fastq_reader(args.fastq_file)
bam_reader(args.bam_file, end, tag)

