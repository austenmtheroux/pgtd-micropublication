#!/usr/bin/env python3
"""Extract gene-wide omega estimates and intervals from FitMG94 JSON output."""

import argparse
import csv
import json
import math
from pathlib import Path


OMEGA_LABEL = "non-synonymous/synonymous rate ratio"


def read_result(path, analysis, purpose, expected_codons):
    with path.open() as handle:
        result = json.load(handle)
    try:
        fit = result["fits"]["Standard MG94"]
        omega = fit["Rate Distributions"][OMEGA_LABEL]
        interval = fit["Confidence Intervals"][OMEGA_LABEL]
        lower, upper = interval["LB"], interval["UB"]
        sequence_count = result["input"]["number of sequences"]
        codon_count = result["input"]["number of sites"]
        analysis_version = str(result["analysis"]["version"])
    except (KeyError, TypeError) as error:
        raise SystemExit(
            f"Unexpected FitMG94 JSON structure in {path}: {error}"
        ) from error

    for label, observed, expected in (
        ("FitMG94 version", analysis_version, "0.4"),
        ("sequence count", sequence_count, 7),
        ("codon count", codon_count, expected_codons),
    ):
        if observed != expected:
            raise SystemExit(f"{path}: expected {label} {expected}; found {observed}")
    if not all(
        type(value) in (int, float) and math.isfinite(value)
        for value in (omega, lower, upper)
    ):
        raise SystemExit(f"Non-finite or non-numeric omega estimate/interval in {path}")
    if not 0 <= lower <= omega <= upper:
        raise SystemExit(f"Invalid omega confidence interval in {path}")

    return {
        "analysis": analysis,
        "purpose": purpose,
        "model": "global_MG94xREV_CF3x4",
        "fitmg94_version": analysis_version,
        "sequence_count": sequence_count,
        "codon_count": codon_count,
        "omega": f"{omega:.12g}",
        "ci_lower": f"{lower:.12g}",
        "ci_upper": f"{upper:.12g}",
    }


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--primary", type=Path, required=True)
    parser.add_argument("--complete", type=Path, required=True)
    parser.add_argument("--output", type=Path, required=True)
    args = parser.parse_args()
    rows = [
        read_result(args.primary, "gene_wide_primary", "reported", 1990),
        read_result(
            args.complete,
            "gene_wide_complete",
            "terminal-trimming sensitivity analysis",
            2027,
        ),
    ]
    args.output.parent.mkdir(parents=True, exist_ok=True)
    with args.output.open("w", newline="") as handle:
        writer = csv.DictWriter(
            handle, fieldnames=tuple(rows[0]), delimiter="\t", lineterminator="\n"
        )
        writer.writeheader()
        writer.writerows(rows)
    for row in rows:
        print(
            f"{row['analysis']}: omega = {float(row['omega']):.3f} "
            f"(95% profile-likelihood CI "
            f"{float(row['ci_lower']):.3f}-{float(row['ci_upper']):.3f}); "
            f"{row['codon_count']:,} codons"
        )


if __name__ == "__main__":
    main()
