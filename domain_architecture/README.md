# Domain architecture

Annotation evidence for the PgtD domain schematic and its comparison with suREJ1 in Figure 1A.

## File locations

| Location | Contents |
| --- | --- |
| `interpro/Dictyostelium_PgtD_representative_proteins/` | InterPro TSV containing results for the seven *Dictyostelium* PgtD orthologs, including *D. discoideum* |
| `interpro/comparators/` | *D. discoideum* DDB_G0289973, PgtB, PgtC, and Gtr1 and *S. purpuratus* suREJ1 InterPro TSVs |
| `interpro/panelA_interpro_annotation_sources.tsv` | annotation sources from the InterPro analyses used to design panel A of Figure 1 |
| `gt_analysis` | *D. discoideum* PgtD GT query subsequences and within-species BLASTP results |
| `figures` | SVG and PNG of panel A |


## InterProScan

Protein sequences were submitted to the InterProScan web server (2026-08-28). Results were downloaded as TSV files and are stored under `interpro/`. InterProScan 5.78-109.0 (synchronized with InterPro 109.0) was the current public release on the analysis date.
Domain annotations were taken from the member-database matches in the TSV files. CATH-Gene3D 4.3.0 was used to annotate the two glycosyltransferase GT-B/GT-A folds. Signal peptide and transmembrane topology annotations used in panel A were taken from Phobius 1.01, and predicted intrinsically disordered regions were taken from MobiDB-lite 2.0.

## References

- Blum M, Andreeva A, Cavalcanti Florentino L, et al. 2025. InterPro: the protein sequence classification resource in 2025. *Nucleic Acids Research* 53:D444–D456. [doi:10.1093/nar/gkae1082](https://doi.org/10.1093/nar/gkae1082).
- Jones P, Binns D, Chang HY, et al. 2014. InterProScan 5: genome-scale protein function classification. *Bioinformatics* 30:1236–1240. [doi:10.1093/bioinformatics/btu031](https://doi.org/10.1093/bioinformatics/btu031).
- Käll L, Krogh A, Sonnhammer ELL. 2004. A combined transmembrane topology and signal peptide prediction method. *Journal of Molecular Biology* 338:1027–1036. [doi:10.1016/j.jmb.2004.03.016](https://doi.org/10.1016/j.jmb.2004.03.016).
- Lewis TE, Sillitoe I, Dawson N, Lam SD, Clarke T, Lee D, Orengo CA, Lees J. 2018. Gene3D: extensive prediction of globular domains in proteins. *Nucleic Acids Research* 46:D435–D439. [doi:10.1093/nar/gkx1069](https://doi.org/10.1093/nar/gkx1069).
- Necci M, Piovesan D, Dosztányi Z, Tosatto SCE. 2017. MobiDB-lite: fast and highly specific consensus prediction of intrinsic disorder in proteins. *Bioinformatics* 33:1402–1404. [doi:10.1093/bioinformatics/btx015](https://doi.org/10.1093/bioinformatics/btx015).

