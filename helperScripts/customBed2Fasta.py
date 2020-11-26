#!/usr/bin/env python
import pandas
import argparse

scriptName = "customBed2Fasta"
Version = 1.2


def bed_to_fasta(file, out):
    """Convert .bed file to .fasta file.

    Loops over the rows of input bed file to extract the sequence
    alone from the seq column and save it to a new fasta file.
    """
    df = pandas.read_csv(
        str(file), sep="\t", header=None, names=["chr", "pos1", "pos2", "seq", "Tm"]
    )
    probe_inx = [">probe_%s_%s" % (str(i), str(out)) for i in range(1, len(df) + 1)]

    new_inx = pandas.Series(probe_inx)
    out_df = pandas.concat([new_inx, df["seq"]], axis=1, sort=False)

    with open(str(out) + ".fasta", "w+") as f:
        for index, row in out_df.iterrows():
            f.write(str(row[0]) + "\n")
            f.write(str(row["seq"]) + "\n")


def _parse_args():
    parser = argparse.ArgumentParser(
        description="Extract sequences from .bed file and save them as .fasta."
    )
    parser.add_argument(
        "-f",
        "--File",
        required=True,
        help=(
            "Input bed file, from which the single sequences will be extracted "
            "and saved as a fasta file. [required]"
        ),
    )
    parser.add_argument(
        "-o", "--out", required=True, help="Output fasta file name. [required]"
    )
    args = parser.parse_args()
    return args


def main():
    # Parse user arguments
    args = _parse_args()
    inFile = args.File
    outSuffix = args.out

    # Call bed_to_fasta to extract sequences
    bed_to_fasta(inFile, outSuffix)


if __name__ == "__main__":
    main()
