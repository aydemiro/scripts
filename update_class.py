import argparse
import pandas as pd


def update_class(class_file, count_file, output_file, count_source):
    dt = {'RTS_stage': str, 'within_CAGE_peak': str, 'polyA_motif_found': str}
    clf = pd.read_table(class_file, dtype=dt)
    cts = pd.read_table(count_file)

    if count_source == "tama":
        cts = cts.rename(columns={"merge_trans_id": "isoform",
                         "trans_read_count": "read_count"})
        cts["num_samples"] = cts["source_line"].apply(
            lambda a: len(a.split(",")))
        cts = cts[["isoform", "read_count", "num_samples", "source_line"]]

    elif count_source == "isoquant":
        cts = cts.rename(columns={"#feature_id": "isoform"})
        cts = cts.set_index("isoform")
        cts = pd.concat([(cts > 0).sum(axis=1), cts.sum(axis=1)], axis=1)
        cts = cts.reset_index()
        cts = cts.rename(columns={0: "num_samples", 1: "read_count"})
        cts = cts[["isoform", "read_count", "num_samples"]]

    clf = clf.merge(cts, how="left")
    clf["read_count"].fillna(0, inplace=True)
    clf["num_samples"].fillna(0, inplace=True)

    trans_5_thr = 50
    trans_3_thr = 50
    gene_5_thr = 50
    gene_3_thr = 50
    a_thr = 60

    cage = clf["within_CAGE_peak"] == "TRUE"
    polya = clf["polyA_motif_found"] == "TRUE"
    trans5 = abs(clf["diff_to_TSS"]) <= trans_5_thr
    trans3 = abs(clf["diff_to_TTS"]) <= trans_3_thr
    gene5 = abs(clf["diff_to_gene_TSS"]) <= gene_5_thr
    gene3 = abs(clf["diff_to_gene_TTS"]) <= gene_3_thr
    perc_a = clf["perc_A_downstream_TTS"] <= a_thr

    five_prime_trans_filter = trans5 | cage
    three_prime_trans_filter = trans3 | (polya & perc_a)
    five_prime_gene_filter = gene5 | cage
    three_prime_gene_filter = gene3 | (polya & perc_a)

    clf["five_prime_trans"] = "FALSE"
    clf.loc[five_prime_trans_filter, "five_prime_trans"] = "TRUE"

    clf["three_prime_trans"] = "FALSE"
    clf.loc[three_prime_trans_filter, "three_prime_trans"] = "TRUE"

    clf["three_prime_gene"] = "FALSE"
    clf.loc[three_prime_gene_filter, "three_prime_gene"] = "TRUE"

    clf["five_prime_gene"] = "FALSE"
    clf.loc[five_prime_gene_filter, "five_prime_gene"] = "TRUE"

    clf.to_csv(output_file, sep="\t", header=True, index=False)

    return


if __name__ == "__main__":
    # Read input arguments
    parser = argparse.ArgumentParser(
        description="""Update SQANTI3 classification file.""")

    parser.add_argument("-c", "--class-file", required=True)
    parser.add_argument("-t", "--count-file", required=True)
    parser.add_argument("-o", "--output-file", required=True)
    parser.add_argument("-s", "--count-source", required=True,
                        choices=["isoquant", "tama"])

    args = vars(parser.parse_args())

    update_class(**args)

