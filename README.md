# vep_FiVar
Script to filter VEP-annotated variants

This script takes a vep-annotated file and outputs the variants passing the input filter settings.
Filter settings: minor allele frequency (0-1), impact ("high", "low", "none")
If no MAF filter is required, use 1. The MAF is taken from the gnomAD AF grpmax - this is the maximum allele frequency across all populations in gnomAD (genomes).
The pathogenicity filter is whether a variant is predicted to be a high impact variant (according to CADD (>25), SpliceAI (>0.5), UTRAnnotator (scoring based on https://academic.oup.com/bioinformatics/article/37/8/1171/5905476 and https://www.nature.com/articles/s41467-019-10717-9) or AlphaMissense (>0.564), a low-impact variant (anything within a gene or which has a high ReMM score (>0.95) outside a gene), or whether no impact filter is applied.

A private allele frequency test is carried out to check the input cohort and filter out any variants present at a frequency over MAF+0.4. This can be deactivated by using the "ukbb" setting, although that also means the calls won't be output.
A biotype filter test is carried out to remove any transcripts clsssified as either a pseudogene or nonsense-mediated decay.

The script requires a header with the IDs in the right order and very specific VEP annotation (see the associated file for an example header output from the VEP, and the annotation sources).
The output will be vcf by default; add "tab" to command line to get tabulated file.
The "ukbb" setting means the script will only output the variants, it won't output all the individual genotypes (important for very large datasets with hundreds of thousands of participants).
The "ukbb" setting will also prevent the script from carrying out the private frequency test.

Gene-specific CADD threshold scores come from PMID: 26820543; use "set-cadd" to keep cadd threshold at 25 for all genes. 
Note that getting gene-specific CADD scores means that only the ~2950 genes with HGMD-generated CADD scores will be included in the output.
A different set of scores can be used (in the same format: tab-separated lines of 11 columns; column 9 has MSC, column 10 has "HGMD", column 11 has the Ensembl gene ID), in which case the output will include only those genes.
