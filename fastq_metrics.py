import argparse
from Bio import SeqIO
from itertools import izip_longest
import os

def fastq_reader(input_list_fastq):

    reads_list = []

    for fastq_line in open(input_list_fastq, "r"):
        pair_reads = fastq_line.strip().split("\t")
        reads_list.append(pair_reads)

    return reads_list

def fastq_metrics(read_list):

    number_pair_reads = []

    for read_path in read_list:
        # Iterates over the forward and the reverse read
        readiter_fw = SeqIO.parse(open(read_path[0]), "fastq")
        readiter_rv = SeqIO.parse(open(read_path[1]), "fastq")

        temp_number_reads_fw = []
        temp_number_reads_rv = []

        temp_number_read_without_NNN_fw = []
        temp_number_read_without_NNN_rv = []

        for rec1, rec2 in izip_longest(readiter_fw, readiter_rv):
            # Obtain the list of ids and count them (lectures number)
            temp_number_reads_fw.append([rec1.id])
            temp_number_reads_rv.append([rec2.id])

            if rec1.seq.startswith("N"):
                pass
            else:
                temp_number_read_without_NNN_fw.append([rec1.id])

            if rec2.seq.startswith("N"):
                pass
            else:
                temp_number_read_without_NNN_rv.append([rec2.id])

        number_pair_reads.append([len(temp_number_reads_fw), len(temp_number_reads_rv),
                                  len(temp_number_read_without_NNN_fw), len(temp_number_read_without_NNN_rv)])


    return number_pair_reads

def calculate_fastq_metrics(input_fastq):

    # Create output directory
    output_dir = "./fastq_metrics"
    fastq_statistics_filename = output_dir + "/" + "fastq_statistics.txt"


    if not os.path.exists("./fastq_metrics"):
        os.makedirs(output_dir)

    fastq_statistics = open(fastq_statistics_filename, "w")

    fastq_statistics.write("#Reads_FW\t#Reads_RV\t#Reads_filterN_FW\t#Reads_filterN_RV\n")

    input_fastq_list = fastq_reader(input_fastq)

    number_reads_list = fastq_metrics(input_fastq_list)

    for n_reads in number_reads_list:
        numbers_reads = str(n_reads[0]) + "\t" + str(n_reads[1]) + "\t" + str(n_reads[2]) + "\t" + str(n_reads[3]) + "\n"
        fastq_statistics.write(numbers_reads)

def main():

    program_description = "WKO"

    parser = argparse.ArgumentParser(description=program_description)

    parser.add_argument("-i", "--input", type=str, help="WKO", required=True)

    args = parser.parse_args()

    fastq_list = args.input

    calculate_fastq_metrics(fastq_list)

if __name__ == "__main__":
    main()