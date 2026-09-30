#!/usr/bin/env bash
# Align seven PgtD orthologs and fit a gene-wide dN/dS model.
set -euo pipefail

ANALYSIS_DIR="$(cd -- "$(dirname -- "${BASH_SOURCE[0]}")" && pwd)"
INPUT_DIR="${INPUT_DIR:-$ANALYSIS_DIR/inputs}"
RESULTS_DIR="${RESULTS_DIR:-$ANALYSIS_DIR/results}"

PROTEIN_FASTA="$INPUT_DIR/Dictyostelium_PgtD_representative_proteins.fasta"
CDS_FASTA="$INPUT_DIR/Dictyostelium_PgtD_representative_CDS.fasta"
SPECIES_TREE="$INPUT_DIR/Dictyostelium_species_tree.nwk"
PROTEIN_ALIGNMENT="$RESULTS_DIR/protein_alignment_mafft_einsi.fasta"
GAP_FILTERED_ALIGNMENT="$RESULTS_DIR/codon_alignment_gap_filtered.fasta"
PRIMARY_ALIGNMENT="$RESULTS_DIR/codon_alignment_primary.fasta"
PRIMARY_JSON="$RESULTS_DIR/fitmg94_gene_wide_primary.json"
COMPLETE_JSON="$RESULTS_DIR/fitmg94_gene_wide_complete_sensitivity.json"
SUMMARY_TSV="$RESULTS_DIR/gene_wide_dnds_summary.tsv"

EXPECTED_HYPHY_ANALYSES_COMMIT="42a3fd041399a17b9bc9ccb64f31f2a2fb764972"
EXPECTED_FITMG94_SHA256="d03b04d93e78e514b7b9005f1e78b371b73e96e44866e108b449be8c4c13e942"

die() {
    printf 'ERROR: %s\n' "$*" >&2
    exit 1
}

sha256_file() {
    if command -v shasum >/dev/null 2>&1; then
        shasum -a 256 "$1" | awk '{print $1}'
    elif command -v sha256sum >/dev/null 2>&1; then
        sha256sum "$1" | awk '{print $1}'
    else
        die "Neither shasum nor sha256sum is available"
    fi
}

for command_name in python3 perl mafft hyphy; do
    command -v "$command_name" >/dev/null 2>&1 \
        || die "$command_name is not available in PATH"
done
for input_file in "$PROTEIN_FASTA" "$CDS_FASTA" "$SPECIES_TREE"; do
    [[ -s "$input_file" ]] || die "Missing required input: $input_file"
done

PAL2NAL_PL="${PAL2NAL_PL:-$(command -v pal2nal.pl 2>/dev/null || true)}"
[[ -n "$PAL2NAL_PL" && -s "$PAL2NAL_PL" ]] \
    || die "pal2nal.pl is not in PATH; set PAL2NAL_PL to its full path"
HYPHY_ANALYSES_DIR="${HYPHY_ANALYSES_DIR:-}"
[[ -n "$HYPHY_ANALYSES_DIR" ]] \
    || die "Set HYPHY_ANALYSES_DIR to the pinned hyphy-analyses checkout"
FITMG94_BF="$HYPHY_ANALYSES_DIR/FitMG94/FitMG94.bf"
[[ -s "$FITMG94_BF" ]] || die "FitMG94 was not found at $FITMG94_BF"

FITMG94_SHA256="$(sha256_file "$FITMG94_BF")"
[[ "$FITMG94_SHA256" == "$EXPECTED_FITMG94_SHA256" ]] \
    || die "FitMG94.bf does not match the pinned SHA-256 checksum"

HYPHY_ANALYSES_COMMIT="not_available"
HYPHY_ANALYSES_STATUS="not_a_git_checkout"
if command -v git >/dev/null 2>&1 \
    && git -C "$HYPHY_ANALYSES_DIR" rev-parse --is-inside-work-tree >/dev/null 2>&1; then
    HYPHY_ANALYSES_COMMIT="$(git -C "$HYPHY_ANALYSES_DIR" rev-parse HEAD)"
    [[ "$HYPHY_ANALYSES_COMMIT" == "$EXPECTED_HYPHY_ANALYSES_COMMIT" ]] \
        || die "hyphy-analyses is not at the pinned commit"
    [[ -z "$(git -C "$HYPHY_ANALYSES_DIR" status --porcelain --untracked-files=no)" ]] \
        || die "hyphy-analyses contains modifications to tracked files"
    HYPHY_ANALYSES_STATUS="tracked_files_clean"
fi

mkdir -p "$RESULTS_DIR"

printf 'Aligning proteins with MAFFT E-INS-i...\n'
mafft --ep 0 --genafpair --maxiterate 1000 --thread 1 "$PROTEIN_FASTA" \
    > "$PROTEIN_ALIGNMENT" 2> "$RESULTS_DIR/mafft.log"

printf 'Projecting CDS onto the alignment with PAL2NAL...\n'
perl "$PAL2NAL_PL" "$PROTEIN_ALIGNMENT" "$CDS_FASTA" -output fasta -nogap \
    > "$GAP_FILTERED_ALIGNMENT" 2> "$RESULTS_DIR/pal2nal.log"

printf 'Checking alignments and removing terminal codons 1991-2027...\n'
python3 "$ANALYSIS_DIR/prepare_alignment.py" \
    --proteins "$PROTEIN_FASTA" \
    --cds "$CDS_FASTA" \
    --protein-alignment "$PROTEIN_ALIGNMENT" \
    --codon-alignment "$GAP_FILTERED_ALIGNMENT" \
    --tree "$SPECIES_TREE" \
    --trim-after 1990 \
    --expected-sequences 7 \
    --expected-gap-filtered-codons 2027 \
    --output "$PRIMARY_ALIGNMENT" \
    --qc-output "$RESULTS_DIR/alignment_qc_summary.tsv"

run_fitmg94() {
    local alignment="$1" output_json="$2" output_log="$3"
    hyphy "$FITMG94_BF" \
        --alignment "$alignment" --tree "$SPECIES_TREE" \
        --rooted No --code Universal --type global \
        --frequencies CF3x4 --lrt No --output "$output_json" \
        > "$output_log" 2>&1
}

printf 'Fitting the primary global MG94 model...\n'
run_fitmg94 "$PRIMARY_ALIGNMENT" "$PRIMARY_JSON" \
    "$RESULTS_DIR/fitmg94_gene_wide_primary.log"
printf 'Fitting the complete-alignment sensitivity model...\n'
run_fitmg94 "$GAP_FILTERED_ALIGNMENT" "$COMPLETE_JSON" \
    "$RESULTS_DIR/fitmg94_gene_wide_complete_sensitivity.log"

python3 "$ANALYSIS_DIR/summarize_fitmg94.py" \
    --primary "$PRIMARY_JSON" --complete "$COMPLETE_JSON" --output "$SUMMARY_TSV"

PAL2NAL_SHA256="$(sha256_file "$PAL2NAL_PL")"
PAL2NAL_VERSION="$(perl -ne '
    if (/^#\s*pal2nal\.pl\s+\(v([0-9.]+)\)/) {
        print $1;
        exit;
    }
' "$PAL2NAL_PL")"

{
    printf 'software\tversion_or_identifier\n'
    printf 'Python\t%s\n' "$(python3 --version 2>&1)"
    printf 'Perl\t%s\n' "$(perl -e 'printf "%vd", $^V')"
    printf 'MAFFT\t%s\n' "$(mafft --version 2>&1 | head -n 1)"
    printf 'HyPhy\t%s\n' "$(hyphy --version 2>&1 | head -n 1)"
    printf 'PAL2NAL\t%s\n' "${PAL2NAL_VERSION:-unknown}"
    printf 'PAL2NAL_SHA256\t%s\n' "$PAL2NAL_SHA256"
    printf 'FitMG94_analysis_version\t0.4\n'
    printf 'FitMG94.bf_SHA256\t%s\n' "$FITMG94_SHA256"
    printf 'hyphy-analyses_commit\t%s\n' "$HYPHY_ANALYSES_COMMIT"
    printf 'hyphy-analyses_status\t%s\n' "$HYPHY_ANALYSES_STATUS"
} > "$RESULTS_DIR/software_versions.tsv"

printf 'Analysis complete. Summary: %s\n' "$SUMMARY_TSV"
