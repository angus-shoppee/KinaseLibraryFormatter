
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
        "-o",
        "--output",
        type=str,
        required=True,
        help="Output file path"
    )

    return parser
