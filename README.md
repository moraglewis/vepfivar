# vep_FiVar
Script to filter VEP-annotated variants

This script takes a vep-annotated file and outputs the variants passing whichever filter is chosen<br>
Filter settings:<br>
0 - no filter, just output (default)<br>
1 - MAF < 10%, no pathogenicity filter<br>
2 - MAF < 0.5%, no pathogenicity filter<br>
3 - MAF < 10%, pathogenicity filter<br>
4 - MAF < 0.5%, pathogenicity filter<br>
5 - no MAF filter, pathogenicity filter<br>

The script requires a header with the IDs in the right order and very specific VEP annotation (see the associated file for an example header output from the VEP, and the annotation sources)<br>
The output will be vcf by default; add "tab" to command line to get tabulated file<br>
The "ukbb" setting means the script will only output the variants, it won't output all the individual genotypes (important for very large datasets with hundreds of thousands of participants)<br>
The "ukbb" setting will also prevent the script from carrying out the private frequency test<br>
The "ukbb" setting will use the maximum MAF from all populations and all databases, otherwise the script will assume the Non-Finnish European setting and will check GnomAD first, then the 1000 Genomes project, then TopMed, then ESP6500<br>
Gene-specific CADD threshold scores come from PMID: 26820543; use "set-cadd" to keep cadd threshold at 25 for all genes<br>
Note that getting gene-specific CADD scores means that only the ~2950 genes with HGMD-generated CADD scores will be included in the output<br>
Pathogenicity prediction filter settings are hardwired: cadd > 25, |sutr| > 1, spliceai > 0.5<br>
