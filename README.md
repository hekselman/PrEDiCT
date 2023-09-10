# Affected cell types for hundreds of Mendelian diseases revealed by analysis of human and mouse single-cell data
This GitHub repository includes codes to redo analyses applied in https://www.biorxiv.org/content/10.1101/2022.10.29.513906v1 . Each file includes functions with docstrings documenting their purpose.

 File name | Language | Methods title 
 ------------- |-------------| -----
 Cell type similarities | R | Determining matching cell types between tissues
 Human expression analysis | R | Human single-cell transcriptomics analysis
 Literature text-mining | Python | Literature text-mining of PubMed records
 Mouse clustering | R | Annotations of mouse cell clusters
 Mouse expression analysis | R | Mouse single-cell transcriptomics analysis
 Permutation tests | R | Permutation tests for similarity between affected cell types
 PrEDiCT score calculation.py | Python | PrEDiCT score calculation

Data needed for analyses are available in Data folder. Additionally, Seurat objects for the six human and mouse tissues are available in:
https://datadryad.org/stash/share/5j6T7Duzcbyx3jNVL_irARhUphpUeW0vuYGTQHofozM .

All outputs of the analysis are available in Output folder.
