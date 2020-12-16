# Initial_single_cell_analysis

This repository contains code for the initial analysis of the CeNGEN 10X 3' single-cell RNA sequencing experiments.

The Functions folder contains custom functions used in the code for analysis and visualization. Many of these came from the Trapnell/Waterston labs, but have modified slightly to work with monocle3 objects. 

Included also is a lookup table that contains cell annotations for all cells that were called by Emptydrops. Cells with greater than 20% of UMIs coming from mitochondrial protein-coding genes (after SoupX background correction) are in this list, but do not have cell or tissue annotations because they were removed before annotation. Annotation of cell and tissue types was determined by expression of many marker genes with known expression patterns in WormBase (wormbase.org) or from our research group. The list of these genes can be found in the supplemental tables of our manuscript https://www.biorxiv.org/content/10.1101/2020.12.15.422897v1.
