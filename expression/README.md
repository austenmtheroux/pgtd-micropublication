# pgtD expression panels B and C

Developmental and cell-type expression from Parikh et al. (2010). Input files and provenance are described in `data/raw/README.md`.

## Run

Requires R and the Source Sans 3 font installed on the OS. Install the R packages if needed:

``` r
install.packages(c("readxl", "ggplot2", "ggtext", "svglite"))
```

From this directory:

``` bash
bash run_analysis.sh
```

The scripts write extracted data to `results/tables/`, panel B and C SVGs and PNGs to `results/figures/`, and R/package versions to `results/session_info.txt`.

## Calculations

Both panels show `log2(scaled abundance + 1)` using the authors' scaled matrices directly. Panel B shows two biological replicates (grey) and their mean on the log-transformed scale (black). Panel C connects paired prestalk and prespore samples. Fold ranges are prespore/prestalk ratios calculated without a pseudocount; P values are the published whole-transcript `ebays.pval` values from Supplementary Tables 5 and 6.
