"""
Author: Ozkan Aydemir
Description: Simple script to create transcript ID to gene ID table from a gtf file.
For use on ghpcc cluster with plastid loaded either as a module or conda env.
"""

import pandas as pd
from plastid import GTF2_Reader
import bz2, gzip
import warnings
warnings.filterwarnings('ignore')
import argparse

# Parse command line arguments
parser = argparse.ArgumentParser()
parser.add_argument("-g", "--gtf", help="Input gtf file.",
                    required=True)
parser.add_argument("-o", "--output", help="Output file",
                    required=True)
parser.add_argument("-z", "--compression", choices=[None, "bz2", "gzip"],
                    default=None)


args = vars(parser.parse_args())
gtf_file = args["gtf"]
compression = args["compression"]
output_file = args["output"]

# Read gtf file using GTF2_Reader from plastid
if compression is None:
    anno = GTF2_Reader(gtf_file)
elif compression == "gzip":
    anno = GTF2_Reader(gzip.open(gtf_file, "rt"))
elif compression == "bz2":
    anno = GTF2_Reader(bz2.open(gtf_file, "rt"))
else:
    raise OSError(
        "Provided compression '{}' is not supported.".format(
            compression))
                    
# Populate a transcript to gene dictionary from the gtf
tx = {}
for feat in anno:
    if feat.attr["type"] == "transcript":
        tx[feat.attr["transcript_id"]] = {
            "gene_id": feat.attr["gene_id"]}

# create a dataframe from the dict
tx_df = pd.DataFrame.from_dict(
    tx, orient="index").reset_index().rename(
        columns={"index": "transcript_id"})

# save dataframe to file
tx_df[["transcript_id", "gene_id"]].to_csv(
    output_file, index=False)

