# PgtD glycosyltransferase-region within-species similarity search

BLASTP searches of *D. discoideum* PgtD glycosyltransferase GT-B and GT-A regions against *D. discoideum* proteome

## Inputs, methods, and results

`D_discoideum_PgtD_GT_slices.fasta` contains two subsequences extracted using one-based inclusive coordinates:

| Sequence | Residues | InterPro annotation source for window |
|----------------------|----------------------------:|----------------------|
| `PgtD_110-501` | 110–501 | CATH-Gene3D G3DSA:3.40.50.2000; candidate GT-B region |
| `PgtD_524_764` | 524–764 | CATH-Gene3D G3DSA:3.90.550.10; candidate GT-A region |

Within-species BLASTP reported that PgtD's GT-B subsequence was highly similar to the GT-B region of Gtr1 (E = 7e-42), while the GT-A subsequence matched the GT-A regions of Q54GU8 (E = 8e-10), PgtB (2e-7), and PgtC (4e-7).

## Files

| Location | Contents |
|----|----|
| `D_discoideum_PgtD_GT_slices.fasta` | GT subsequence queries |
| `blastp/` | BLASTP reports |
| `../interpro/` | InterProScan TSV files |