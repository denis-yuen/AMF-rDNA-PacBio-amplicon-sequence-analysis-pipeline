import argparse
import os
from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord

def main():
    parser = argparse.ArgumentParser(description="""This script will reverse complement the sequences in a fasta file. Example command: python reverse_complement_fasta_file.py my_file.fasta""")

    ## required parameters for all types of data
    parser.add_argument("input_file", type=str, help="Specify a fasta file for your analysis.  Example: myInputFile.fasta")
    parser.add_argument("output_file", type=str, help="Specify name for output file for your analysis.  Example: myOutputFile.fasta")

    ## parse the arguments
    args = parser.parse_args()
    input_file = args.input_file
    output_file = args.output_file

## Check if the input file exists
    if not os.path.exists(input_file):
        print("File" + input_file + " was not found.")
        exit()

    file_out = open(output_file, "w")

    for seq_record in SeqIO.parse(input_file, "fasta"):
        this_id=seq_record.id
        this_seq=str(seq_record.seq.reverse_complement())

        rc_record = SeqRecord( Seq(this_seq), id=this_id )

        SeqIO.write(rc_record, file_out, 'fasta')

if __name__ == "__main__":
    main()
