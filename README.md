## PgtD micropublication 

This repository contains the analysis scripts, source data, results, and figures supporting the micropublication on *Dictyostelium discoideum* PgtD, a predicted glycosyltransferase–REJ–GAIN protein.

The analyses outline PgtD domain architecture, identify homologs across the *Dictyostelium* genus, estimate gene-wide dN/dS among representative homologs, and reanalyze published developmental and cell-type expression data.

# Repository organization

| Location | Contents |
| --- | --- |
| `domain_architecture` | InterPro annotations, glycosyltransferase-region analyses, comparator proteins, and source files for Figure 1A
| `homolog_identification` | BLASTP searches and selection of representative PgtD homologs from seven *Dictyostelium* species
| `expression/` | Reanalysis of developmental and prespore/prestalk expression data used for Figure 1B-C
| `gene_wide_dnds/` | Protein/codon alignment workflow and HyPhy FitMG94 analysis used to estimate gene-wide dN/dS
| `figures/` | Final assembled Figure 1 export files and original Adobe Illustrator file

Individual figure panel files are located with the analyses that generated/supported them. Figure 1A files are under `domain_architecture/figures` and Figure 1B-C files under `expression/results/figures`. 