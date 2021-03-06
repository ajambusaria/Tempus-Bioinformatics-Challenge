# Tempus Bioinformatics Challenge

## Author
Ankit Jambusaria

## Summary
Provided a VCF file, the program outputs a table annotating each variant in the file. 
Each variant is annotated with the following pieces of information:
1. Type of variation (Substitution, Insertion, Silent, Intergenic, etc.) If there are multiple possibilities, the annotation includes the most deleterious possibility.
2. Depth of sequence coverage at the site of variation.
3. Number of reads supporting the variant.
4. Percentage of reads supporting the variant versus those supporting reference reads.
5. Allele frequency of variant from Broad Institute ExAC Project API
(API documentation is available here: http://exac.hms.harvard.edu/)
6. Additional relevant information from ExAC such as Ensembl_ID and Consequences.

## Python Dependencies
Python3.7
argparse 
os
pandas
collections
requests

## Tempus_VCF_annot Command Options

*	--input : Provide Input File For Parser
*	--output : Provide Output Folder For Annotation


### Example Run 

python Tempus_VCF_annot.py --input Challenge_data.vcf --output Output

annotated_Challenge_data.vcf in the Output folder is the annotated output file 
