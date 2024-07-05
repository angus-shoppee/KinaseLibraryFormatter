
from argparse import ArgumentParser


def get_arg_parser() -> ArgumentParser:

    parser = ArgumentParser()

    parser.add_argument(
        "dataset_path",
        type=str,
        help="Path to the dataset to be processed"
    )

    parser.add_argument(
        "-s",
        "--sequence",
        type=str,
        required=True,
        help="Name of the column containing phospho-peptide sequence, including phosphorylation probability scores"
    )

    parser.add_argument(
        "-x",
        "--position",
        type=str,
        required=True,
        help="Name of the column containing position of phospho-site in peptide sequence"
    )

    parser.add_argument(
        "-f",
        "--fold-change",
        type=str,
        required=True,
        help="Name of the column containing fold change"
    )

    parser.add_argument(
        "-p",
        "--p-value",
        type=str,
        required=True,
        help="Name of the column containing p values"
    )

    parser.add_argument(
        "-g",
        "--gene-name",
        type=str,
        required=False,
        help="(If using the --protein-info flag) Name of the column containing gene names"
    )

    parser.add_argument(
        "-i",
        "--protein-id",
        type=str,
        required=False,
        help="(If using the --protein-info flag) Name of the column containing protein IDs"
    )

    parser.add_argument(
        "-o",
        "--output",
        type=str,
        required=True,
        help="Output file path"
    )

    parser.add_argument(
        "--header",
        action="store_true",
        help="(Optional) Include header row from output file"
    )

    parser.add_argument(
        "--protein-info",
        action="store_true",
        help="(Optional) Include protein ID and gene name columns from output file (requires -p and -i arguments)"
    )

    # parser.add_argument(
    #     "--disable-scientific-notation",
    #     action="store_true",
    #     help="(Optional) Disable use of scientific notation in output for very small or very large numbers"
    # )

    return parser
