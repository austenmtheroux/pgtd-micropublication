#!/usr/bin/env python3
"""Check the PgtD sequence inputs and trim the gap-filtered codon alignment."""

import argparse
import re
from pathlib import Path


DNA_ALPHABET = set("ACGT")
PROTEIN_ALPHABET = set("ACDEFGHIKLMNPQRSTVWY")
CODONS = (a + b + c for a in "TCAG" for b in "TCAG" for c in "TCAG")
AMINO_ACIDS = (
    "FFLLSSSSYY**CC*W"
    "LLLLPPPPHHQQRRRR"
    "IIIMTTTTNNKKSSRR"
    "VVVVAAAADDEEGGGG"
)
GENETIC_CODE = dict(zip(CODONS, AMINO_ACIDS))


def read_fasta(path):
    records = {}
    identifier = None
    for raw_line in path.read_text().splitlines():
        line = raw_line.strip()
        if not line:
            continue
        if line.startswith(">"):
            header = line[1:].split()
            if not header:
                raise ValueError(f"Empty FASTA identifier in {path}")
            identifier = header[0]
            if identifier in records:
                raise ValueError(f"Duplicate FASTA identifier in {path}: {identifier}")
            records[identifier] = []
        else:
            if identifier is None:
                raise ValueError(f"Sequence encountered before a header in {path}")
            records[identifier].append(line)

    if not records:
        raise ValueError(f"No FASTA records found in {path}")
    for identifier, lines in records.items():
        records[identifier] = "".join(lines).upper()
        if not records[identifier]:
            raise ValueError(f"Empty sequence in {path}: {identifier}")
    return records


def write_fasta(records, path):
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w") as handle:
        for identifier, sequence in records.items():
            handle.write(f">{identifier}\n")
            for start in range(0, len(sequence), 80):
                handle.write(sequence[start : start + 80] + "\n")


def translate(dna, identifier):
    if len(dna) % 3:
        raise ValueError(f"CDS length is not divisible by three: {identifier}")
    invalid = set(dna) - DNA_ALPHABET
    if invalid:
        raise ValueError(f"Invalid CDS symbols for {identifier}: {sorted(invalid)}")
    protein = "".join(GENETIC_CODE[dna[i : i + 3]] for i in range(0, len(dna), 3))
    if "*" in protein:
        raise ValueError(f"Stop codon in CDS for {identifier}")
    return protein


def require_identifiers(label, observed, expected):
    observed, expected = set(observed), set(expected)
    if observed != expected:
        raise ValueError(
            f"{label} identifier mismatch; missing: {sorted(expected - observed)}; "
            f"extra: {sorted(observed - expected)}"
        )


def aligned_length(records, label):
    lengths = {len(sequence) for sequence in records.values()}
    if len(lengths) != 1:
        raise ValueError(f"Unequal sequence lengths in {label}: {sorted(lengths)}")
    return lengths.pop()


def tree_tips(path):
    # The supplied Newick uses unquoted alignment identifiers as tip labels.
    tips = re.findall(r"(?:^|[(,])\s*([^():,;\s]+)", path.read_text().strip())
    if len(tips) != len(set(tips)):
        raise ValueError("Duplicate tip identifier in the species tree")
    return tips


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--proteins", type=Path, required=True)
    parser.add_argument("--cds", type=Path, required=True)
    parser.add_argument("--protein-alignment", type=Path, required=True)
    parser.add_argument("--codon-alignment", type=Path, required=True)
    parser.add_argument("--tree", type=Path, required=True)
    parser.add_argument("--trim-after", type=int, required=True)
    parser.add_argument("--expected-sequences", type=int, required=True)
    parser.add_argument("--expected-gap-filtered-codons", type=int, required=True)
    parser.add_argument("--output", type=Path, required=True)
    parser.add_argument("--qc-output", type=Path, required=True)
    args = parser.parse_args()

    proteins = read_fasta(args.proteins)
    cds_records = read_fasta(args.cds)
    protein_alignment = read_fasta(args.protein_alignment)
    codon_alignment = read_fasta(args.codon_alignment)
    if len(proteins) != args.expected_sequences:
        raise ValueError(
            f"Expected {args.expected_sequences} proteins; found {len(proteins)}"
        )
    require_identifiers("CDS", cds_records, proteins)
    require_identifiers("protein-alignment", protein_alignment, proteins)
    require_identifiers("codon-alignment", codon_alignment, proteins)
    require_identifiers("tree-tip", tree_tips(args.tree), proteins)

    for identifier, protein in proteins.items():
        invalid = set(protein) - PROTEIN_ALPHABET
        if invalid:
            raise ValueError(
                f"Invalid protein symbols for {identifier}: {sorted(invalid)}"
            )
        if translate(cds_records[identifier], identifier) != protein:
            raise ValueError(f"CDS translation does not match protein: {identifier}")

    protein_alignment_columns = aligned_length(protein_alignment, "protein alignment")
    for identifier, aligned_protein in protein_alignment.items():
        if aligned_protein.replace("-", "") != proteins[identifier]:
            raise ValueError(
                f"MAFFT alignment does not reproduce input protein: {identifier}"
            )

    retained_columns = [
        column
        for column in range(protein_alignment_columns)
        if all(sequence[column] != "-" for sequence in protein_alignment.values())
    ]
    codon_alignment_nt = aligned_length(codon_alignment, "codon alignment")
    if codon_alignment_nt % 3:
        raise ValueError("Codon-alignment length is not divisible by three")
    gap_filtered_codons = codon_alignment_nt // 3
    if gap_filtered_codons != args.expected_gap_filtered_codons:
        raise ValueError(
            f"Expected {args.expected_gap_filtered_codons} gap-filtered codons; "
            f"found {gap_filtered_codons}"
        )
    if gap_filtered_codons != len(retained_columns):
        raise ValueError("PAL2NAL output length differs from gap-free protein columns")

    for identifier, aligned_cds in codon_alignment.items():
        expected_protein = "".join(
            protein_alignment[identifier][column] for column in retained_columns
        )
        if translate(aligned_cds, identifier) != expected_protein:
            raise ValueError(
                f"PAL2NAL output does not match protein alignment: {identifier}"
            )

    if not 0 < args.trim_after <= gap_filtered_codons:
        raise ValueError(f"--trim-after must be between 1 and {gap_filtered_codons}")
    primary_alignment = {
        identifier: sequence[: args.trim_after * 3]
        for identifier, sequence in codon_alignment.items()
    }
    write_fasta(primary_alignment, args.output)

    metrics = [
        ("protein_sequence_count", len(proteins)),
        ("cds_sequence_count", len(cds_records)),
        ("protein_cds_identifiers_match", "true"),
        ("cds_translation_match", "true"),
        ("tree_tip_identifiers_match", "true"),
        ("protein_alignment_columns", protein_alignment_columns),
        ("gap_columns_removed", protein_alignment_columns - len(retained_columns)),
        ("gap_filtered_codons", gap_filtered_codons),
        ("primary_alignment_codons", args.trim_after),
        ("terminal_codons_removed", gap_filtered_codons - args.trim_after),
        ("ambiguous_nucleotide_check", "pass"),
        ("stop_codon_check", "pass"),
    ]
    args.qc_output.parent.mkdir(parents=True, exist_ok=True)
    with args.qc_output.open("w") as handle:
        handle.write("metric\tvalue\n")
        for metric, value in metrics:
            handle.write(f"{metric}\t{value}\n")
    print(
        f"Validated {len(proteins)} sequence pairs; wrote "
        f"{args.trim_after:,} codons to {args.output}"
    )


if __name__ == "__main__":
    main()
