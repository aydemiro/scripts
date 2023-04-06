from plastid import GTF2_TranscriptAssembler
import warnings
warnings.filterwarnings('ignore')
import argparse
import pickle

parser = argparse.ArgumentParser()
parser.add_argument("-g", "--gtf", required=True)
parser.add_argument("-o", "--output", required=True)

args = vars(parser.parse_args())

gtf = args["gtf"]
output = args["output"]

transcripts = list(GTF2_TranscriptAssembler(gtf))

with open(output, "wb") as outfile:
    pickle.dump(transcripts, outfile)

