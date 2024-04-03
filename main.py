
from typing import List, Union
import os
import csv

from parser import get_arg_parser


PHOSPHO_SITE_SCORE_THRESHOLD = 0.95

DEMO_STRING_1 = "PSVEPPLS(1)QETFSDL"
# "PSVEPPLS*QETFSDL"
DEMO_STRING_2 = "QMEVLY(1)AWEFLS(0.95)FAD"
# "QMEVLYAWEFLS*FAD"
# "QMEVLY*AWEFLSFAD"


def format_peptide(
    peptide: str,
    score_threshold: float
) -> List[str]:

    peptide_with_asterisks = ""
    score_buffer = ""
    for char in peptide:
        if char.isalpha():
            if score_buffer != "":
                if float(score_buffer.replace("(", "").replace(")", "")) >= score_threshold:
                    peptide_with_asterisks += "*"
                score_buffer = ""
            peptide_with_asterisks += char
        else:
            score_buffer += char

    segments = peptide_with_asterisks.split("*")

    return [peptide_with_asterisks] if len(segments) <= 2 else [
        "".join(segments[:i+1]) + "*" + "".join(segments[i+1:])
        for i in range(len(segments) - 1)
    ]


def choose_peptide(
    peptides: List[str],
    position: int
) -> Union[None, str]:

    for peptide in peptides:

        if position == len([char for char in peptide if char.isalpha()]):
            continue  # If position is the final amino acid and no subsequent asterix is present, then skip

        aa_index = 0
        for char in peptide:
            if char.isalpha():
                if aa_index + 1 == position:
                    if peptide[aa_index + 1] == "*":  # Check next character after current index
                        return peptide
                aa_index += 1

    return None

    # raise ValueError(f"No phosphorylation site (\"*\") found at position {position} in peptides: {peptides}")


def main() -> None:

    args = get_arg_parser().parse_args()

    print(args)

    if not os.path.exists(args.dataset_path):
        raise ValueError(f"Invalid path (file does not exist): {args.dataset_path}")

    with open(args.dataset_path, "r", encoding="utf-8-sig") as f:

        dataset = csv.reader(f)

        header = next(dataset)
        sequence_index = header.index(args.sequence)
        position_index = header.index(args.position)
        fc_index = header.index(args.fold_change)
        p_index = header.index(args.p_value)

        with open(args.output, "w") as output_file:

            writer = csv.writer(output_file)

            writer.writerow(["peptide", "fc", "p"])

            # n_lines = 0
            while dataset:

                try:
                    row = next(dataset)
                except StopIteration:
                    break

                print("Sequence:", row[sequence_index])
                peptide_options = format_peptide(row[sequence_index], score_threshold=PHOSPHO_SITE_SCORE_THRESHOLD)
                print("-->", peptide_options)

                # if len([True for option in peptide_options if "*" in option]) == 0:
                #     continue  # Skip rows where no phospho-sites are above probability score threshold

                position = int(row[position_index])
                print("Position:", position)
                peptide = choose_peptide(peptide_options, position)
                print("Selected peptide:", peptide)

                print()

                if peptide is None:
                    continue  # Skip rows where phospho-site is below probability score threshold

                writer.writerow([peptide, row[fc_index], row[p_index]])

                # n_lines += 1

    # raw_peptides = [
    #     DEMO_STRING_1,
    #     DEMO_STRING_2
    # ]
    #
    # formatted_peptides = [
    #     format_peptide(peptide, score_threshold=PHOSPHO_SITE_SCORE_THRESHOLD)
    #     for peptide in raw_peptides
    # ]


if __name__ == "__main__":

    main()
