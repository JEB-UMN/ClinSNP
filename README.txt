Welcome to ClinSNP! 

Mission:
This platform was designed to cross-reference an individual's SNPs with annotated entries in NCBI's ClinVar database, allowing for rapid identification of clinically-relevant alleles for particular phenotypes. 

Set-up:
After moving the ClinSNP folder to your documents, transfer the SNP data you would like to analyze into the ClinSNP folder. This is typically provided by personal genomics companies as an "xxx.SNPs.vcf" file. Open the ClinSNP.py file in a python interpreter, then run it to perform the analysis. Database download, along with the creation of output files, are performed automatically by the script.

Considerations:
-ClinSNP will likely take several minutes to run, although subsequent analyses after generating an "aligned_SNPs.csv" file should not take more than a few seconds. It also requires a 200 MB database download, which unzips into a 2 GB txt file, although both of these can be deleted after the aligned_SNPs file has been generated.
-ClinSNP should be compatible with whole genome and whole exome sequencing data.
-ClinSNP was generating using 23andMe whole genome sequencing data from 2014. As a result, some formatting parameters may be outdated. Feel free to change the file import portion of the script to suit your needs, including column names.
-ClinSNP should be compatible with SNP files from other personal genome analysis companies, as long as you adjust the script to alter the import parameters accordingly.
-If you would like to analyze SNPs from sequencing data that does not have a corresponding SNP file, you will have to generate an "xxx.SNPs.vcf file through reference genome alignment first. Other platforms exist that will perform this analysis for you. I am planning to directly incorporate this functionality in a later version of ClinSNP.
-Remember: Genomic data suggests predispositions, not concrete outcomes. A single SNP is almost never enough to cause disease.

For troubleshooting purposes or any additional questions, please contact Jacob Bridge at bridg245@umn.edu.

Thank you for using ClinSNP!