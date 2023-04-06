import warnings
warnings.filterwarnings('ignore')
import pickle
import twobitreader
from plastid import Transcript
from Bio import SeqIO
from Bio.Seq import Seq
import argparse

parser = argparse.ArgumentParser()
parser.add_argument("-f", "--fasta", required=True)
parser.add_argument("-g", "--genome", required=True)
parser.add_argument("-t", "--transcripts", required=True)
parser.add_argument("-o", "--output", required=True)
parser.add_argument("-b", "--ignore-start-codons", required=True, type=float)
parser.add_argument("-e", "--ignore-end-codons", required=True, type=float)
parser.add_argument("-i", "--include-tag", nargs="*")
parser.add_argument("-z", "--offset", required=True, type=int)
parser.add_argument("-l", "--read-length", required=True, type=int)


args = vars(parser.parse_args())

two_bit_genome = twobitreader.TwoBitFile(args["genome"])


transcript_tag = args["include_tag"] 
if transcript_tag is not None:
    transcript_tag = set(transcript_tag)

with open(args["transcripts"], "rb") as infile:
    transcripts = pickle.load(infile)

transcriptome_clip_5 = int(
    round(3 * args["ignore_start_codons"])) - args["offset"]

transcriptome_clip_3 = int(
    round(3 * args["ignore_end_codons"])) +  args["offset"] -  args["read_length"]


seq_list = []

for tra in transcripts:
    if transcript_tag is not None:
        try:
            attrbs = tra.attr["tag"].split(",")
            if len(transcript_tag.intersection(attrbs)) < 1:
                continue
        except KeyError:
            continue

    cds = tra.get_cds()
    u5 = tra.get_utr5()
    u3 = tra.get_utr3()
    if len(cds) > 0:
        cds_seq = cds.get_sequence(two_bit_genome)
        if transcriptome_clip_5 < 0:
            if len(u5) > 0:
                u5_seq = u5.get_sequence(two_bit_genome)[
                    transcriptome_clip_5:]
                cds_seq = u5_seq + cds_seq
        else:
            cds_seq = cds_seq[transcriptome_clip_5:]
        if transcriptome_clip_3 < 0:
            if len(u3) > 0:
                u3_seq = u3.get_sequence(two_bit_genome)[
                    : -transcriptome_clip_3]
                cds_seq = cds_seq + u3_seq
        elif transcriptome_clip_3 > 0:
            cds_seq = cds_seq[: -transcriptome_clip_3]

        seq_list.append(
            SeqIO.SeqRecord(Seq(cds_seq), 
                            tra.attr["transcript_id"],
                            "", ""))

with open(args["output"], "w") as outfile:
    SeqIO.write(seq_list, outfile, "fasta")

