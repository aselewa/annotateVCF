### Description

This repository contains a VCF parser written in Python 3. The script parses VCF v4.1 and outputs a tab-delimited text file with at least 5 columns:

*variant_type*: contains one or two pieces of information depending on whether the variant is found in ExAC. If the variant is not in ExAC, the type is parsed from the VCF file and can be one of snp, del, ins, complex, or combinations of said types. If the variant is found in ExAC, the most severe consequence of the variant is also included in this column.  

*read_depth*: total coverage at variant position.   

*alt_count*: total number of reads supporting the alternative allele.  

*pct_variant_reads*:  percent of reads supporting alternative allele.

*exac_allele_freq*: population allele frequency of variant according to ExAC. This value is NaN if the variant is not in ExAC.

If the VCF has multiple samples, then each sample will contain its own columns of *read_depth*, *alt_count*, and *pct_variant_reads*.

### Usage

*Command line*

```{bash}
python annotateVCF.py INPUT.vcf OUTPUT_FILE_NAME.txt
```

*Jupyter Notebook*

```{python}
import VCFReader as vr

fname = 'path/to/vcf_file/'
output = 'path/to/out.txt'

vr.main(fname, output)
```