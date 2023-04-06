# vepfivar
Script to filter VEP-annotated variants

This script takes a vep-annotated file and outputs the variants passing whichever filter is chosen
Filter settings:
0 - no filter, just output (default)
1 - MAF < 10%, no pathogenicity filter
2 - MAF < 0.5%, no pathogenicity filter
3 - MAF < 10%, pathogenicity filter
4 - MAF < 0.5%, pathogenicity filter
5 - no MAF filter, pathogenicity filter

The script requires a header with the IDs in the right order and very specific VEP annotation (see the associated file for an example header output from the VEP, and the annotation sources)
The output will be vcf by default; add "tab" to command line to get tabulated file
The "ukbb" setting means the script will only output the variants, it won't look at all the individual genotypes (important for very large datasets with hundreds of thousands of participants)
The "ukbb" setting will also prevent the script from carrying out the private frequency test
The "ukbb" setting will use the maximum MAF from all populations and all databases, otherwise the script will assume the Non-Finnish European setting and will check GnomAD first, then the 1000 Genomes project, then TopMed, then ESP6500
Gene-specific CADD threshold scores come from PMID: 26820543; use "set-cadd" to keep cadd threshold at 25 for all genes.
Note that getting gene-specific CADD scores means that only the ~2950 genes with HGMD-generated CADD scores will be included in the output
Pathogenicity prediction filter settings are hardwired: cadd > 25, |sutr| > 1, spliceai > 0.5
