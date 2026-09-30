# Gene-wide dN/dS analysis of Dictyostelium PgtD

Global dN/dS estimation for the seven *Dictyostelium* PgtD homologs, with a sensitivity analysis for terminal alignment trimming. Homolog selection is documented in `../homolog_identification/`.

## Inputs

`inputs/` contains:

- `Dictyostelium_PgtD_representative_proteins.fasta`
- `Dictyostelium_PgtD_representative_CDS.fasta`
- `Dictyostelium_species_tree.nwk`

The protein and CDS files contain the same seven representatives, and terminal stop codons were removed from the CDS sequences. The species-tree topology was pruned from the six-protein AAPPRS phylogeny of Schilde et al. (2019, Fig. 2A), with tips relabeled to match the sequence identifiers. Published branch lengths were omitted. FitMG94 treats the topology as unrooted and estimates branch lengths from the codon alignment.

## Running the analysis

The recorded run used MAFFT 7.526, PAL2NAL 14.1, HyPhy 2.5.78, FitMG94 0.4, Python 3.12.2, and Perl 5.32.1. Versions and checksums are recorded in `results/software_versions.tsv`. The Python scripts require no external packages.

Make `python3`, `perl`, `mafft`, `hyphy`, and `pal2nal.pl` available in `PATH`. Clone the pinned FitMG94 source outside this directory, replacing `/path/to/hyphy-analyses` with your chosen location:

``` bash
git clone https://github.com/veg/hyphy-analyses.git /path/to/hyphy-analyses
git -C /path/to/hyphy-analyses checkout 42a3fd041399a17b9bc9ccb64f31f2a2fb764972
```

From `gene_wide_dnds/`, run:

``` bash
HYPHY_ANALYSES_DIR=/path/to/hyphy-analyses bash run_analysis.sh
```

If PAL2NAL is not in `PATH`, also set `PAL2NAL_PL=/path/to/pal2nal.pl`. The runner checks the FitMG94 checksum and (for a Git checkout) the commit and absence of tracked changes. Rerunning overwrites the corresponding files in `results/`.

## Alignment and results

Proteins were aligned with MAFFT E-INS-i (`--ep 0 --genafpair --maxiterate 1000 --thread 1`). PAL2NAL projected the CDS onto the protein alignment; `-nogap` removed codon columns containing gaps, leaving 2,027 codons.

Gap-filtered alignment codons 1991–2027 aligned a highly divergent *D. purpureum* segment with a conserved terminal segment in the other six species. These 37 ambiguous columns were excluded from all sequences for the primary analysis. The complete gap-filtered alignment was retained as a sensitivity analysis.

| Alignment                    | Codons |  dN/dS | 95% CI        |
|------------------------------|-------:|-------:|---------------|
| Terminally trimmed (primary) |  1,990 | 0.1087 | 0.1041–0.1136 |
| Complete gap-filtered        |  2,027 | 0.1077 | 0.1031–0.1124 |

Terminal trimming had little effect on the estimate. `results/` contains the protein alignment, both codon alignments, FitMG94 JSON outputs, and program logs. `summarize_fitmg94.py` writes estimates and confidence intervals to `gene_wide_dnds_summary.tsv`. `alignment_qc_summary.tsv` records validation and filtering statistics.

## References

- Katoh K, Standley DM. 2013. MAFFT multiple sequence alignment software version 7: improvements in performance and usability. *Molecular Biology and Evolution* 30:772–780. [doi:10.1093/molbev/mst010](https://doi.org/10.1093/molbev/mst010)
- Suyama M, Torrents D, Bork P. 2006. PAL2NAL: robust conversion of protein sequence alignments into the corresponding codon alignments. *Nucleic Acids Research* 34:W609–W612. [doi:10.1093/nar/gkl315](https://doi.org/10.1093/nar/gkl315)
- Kosakovsky Pond SL et al. 2020. HyPhy 2.5—a customizable platform for evolutionary hypothesis testing using phylogenies. *Molecular Biology and Evolution* 37:295–299. [doi:10.1093/molbev/msz197](https://doi.org/10.1093/molbev/msz197)
- Muse SV, Gaut BS. 1994. A likelihood approach for comparing synonymous and nonsynonymous nucleotide substitution rates, with application to the chloroplast genome. *Molecular Biology and Evolution* 11:715–724. [doi:10.1093/oxfordjournals.molbev.a040152](https://doi.org/10.1093/oxfordjournals.molbev.a040152)
- Schilde C, Lawal HM, Kin K, Shibano-Hayakawa I, Inouye K, Schaap P. 2019. A well supported multi gene phylogeny of 52 dictyostelia. *Molecular Phylogenetics and Evolution* 134:66–73. [doi:10.1016/j.ympev.2019.01.017](https://doi.org/10.1016/j.ympev.2019.01.017)
