# PgtD homolog identification

Identification of seven representative *Dictyostelium* PgtD proteins.

## Search and selection

The *D. discoideum* PgtD sequence in `D_discoideum_PgtD.fasta` (UniProt Q54XQ9; GenBank EAL67967.1) was searched against NCBI nr using BLASTP 2.17.0+, restricted to *Dictyostelium* (taxid 5782), on 2026-09-06. The E-value threshold was 1e-5.

`PgtD_BLASTp_nr_Dictyostelium.xml` is the original search output. `PgtD_BLASTp_hit_selection.tsv` records the 40 returned hits and the curation decisions. Hits with at least 90% query coverage were considered near-full-length candidates; one representative per species was selected manually (sequences in `Dictyostelium_PgtD_representative_proteins.fasta`).

The selected protein sequences were submitted to the InterProScan web server to confirm whether they have the expected PgtD architecture, and the resulting TSV files were downloaded (see `../domain_architecture/interpro/Dictyostelium_PgtD_representative_proteins`).

BLASTP searches identified one near-full-length PgtD-like protein in each sampled Dictyostelium proteome; no additional distinct near-full-length paralogs were detected among the hits. In reverse BLASTP searches against the *D. discoideum* UniProt reference proteome (ID: UP000002195), PgtD (UniProt Q54XQ9) was the highest-scoring hit for all six non-D. discoideum candidates, with 97.9–100% query coverage (see `PgtD_homologs_BLASTp_D_discoideum_proteome.xml`).