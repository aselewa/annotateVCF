import re
import sys
import json
import urllib.request     

import pandas as pd # used for outputting to csv file

ENSEMBL_SEVERITY_ORDER = [
    'transcript_ablation',
    'splice_acceptor_variant',
    'splice_donor_variant',
    'stop_gained',
    'frameshift_variant',
    'stop_lost',
    'start_lost',
    'transcript_amplification',
    'inframe_insertion',
    'inframe_deletion',
    'missense_variant',
    'protein_altering_variant',
    'splice_region_variant',
    'incomplete_terminal_codon_variant',
    'stop_retained_variant',
    'synonymous_variant',
    'coding_sequence_variant',
    'mature_miRNA_variant',
    '5_prime_UTR_variant',
    '3_prime_UTR_variant',
    'non_coding_transcript_exon_variant',
    'intron_variant',
    'NMD_transcript_variant',
    'non_coding_transcript_variant',
    'upstream_gene_variant',
    'downstream_gene_variant',
    'TFBS_ablation',
    'TFBS_amplification',
    'TF_binding_site_variant',
    'regulatory_region_ablation',
    'regulatory_region_amplification',
    'feature_elongation',
    'regulatory_region_variant',
    'feature_truncation',
    'intergenic_variant'
    ]

class exac(object):
    '''
    Object representing query results from ExAC database. 
    '''
    def __init__(self,variantIDs):

        url = 'http://exac.hms.harvard.edu/rest/bulk/variant' 

        data = json.dumps(variantIDs)
        data = data.encode('ASCII')

        exacRequest = urllib.request.Request(url, data=data, method='POST')
        exac = urllib.request.urlopen(exacRequest)
        queryResult = json.loads(exac.read().decode())
        
        self.variantIDs = variantIDs
        self.queryResult = queryResult
    
    def allele_freq(self):
        '''
        Returns allele frequency from query result for each varient. 
        If the variant is not found, it returns nan
        '''
        AF = []
        for i in range(len(self.variantIDs)):
            curr = self.queryResult[self.variantIDs[i]]['variant']
            if 'allele_freq' in curr.keys():
                AF.append([curr['allele_freq']])
            else:
                AF.append([float("nan")])
        return AF

    def consequence(self):
        '''
        Returns most severe consequence from ExAC query results.
        If the variant is not found, it returns an empty string.
        '''
        anno = []
        for i in range(len(self.variantIDs)):
            curr = self.queryResult[self.variantIDs[i]]
            if curr['consequence']:
                cons = list(curr['consequence'].keys())
                severity = []
                for c in cons:
                    if c in ENSEMBL_SEVERITY_ORDER:
                        severity.append(ENSEMBL_SEVERITY_ORDER.index(c))
                    else:
                        severity.append(float("inf"))
                idx = min([(idx, val) for (val, idx) in enumerate(severity)])[1] # index at which there is maximum severity
                anno.append(cons[idx])
            else:
                anno.append('')
        return anno
    
class VCFReader(object):
    '''
    Object representing a VCF file
    '''
    def __init__(self, fname):
        inVCF = []
        file = open(fname, 'r')
        for line in file:
            if '#CHROM' in line: 
                keys = line.strip('\n').split('\t')
                samples = re.findall('FORMAT\t(.+)', line)[0].split('\t')
                header = ['variant_type'] # header for annotated VCF file
                for s in samples:
                    header += [s + ' read_depth', s + ' alt_counts', s + ' pct_variant_reads']
                header += ['exac_allele_freq']
            if line[0] != '#':
                currLine = line.rstrip('\n').split('\t')
                inVCF.append(currLine)
        file.close()

        self.keys = keys
        self.vcf = inVCF
        self.samples = samples
        self.annotation_header = header


    def extractType(self):
        '''
        Extracts the variant type from the VCF file for each variant
        '''   
        varType = []
        for entry in self.vcf:
            currDict = dict(zip(self.keys, entry))
            currInfo = currDict["INFO"]
            currAnno = re.findall('TYPE=(.+)',currInfo)
            varType.append(currAnno)
        return varType

    def extractAD(self):
        '''
        Returns the total read depth, alternative read depth, and computes alternative allele frequency for each sample and each variant.
        '''          
        perSampleAD = []
        for entry in self.vcf:
            currDict = dict(zip(self.keys, entry))
            AD = []
            for i in range(len(self.samples)):
                sampleInfo = currDict[self.samples[i]]
                depth = int(sampleInfo.split(':')[3].split(',')[0])    # Third entry contains counts per allele
                altcount = int(sampleInfo.split(':')[3].split(',')[1])
                AF = altcount/depth
                AD += [depth, altcount, AF]
            perSampleAD.append(AD)
        return perSampleAD
    
    def extractVariantIDs(self):
        '''
        Returns a variant ID composed of 'chromosome-position-reference-alternative' which is used in performing ExAC queries
        '''  
        variantIDs = []
        for entry in self.vcf:
            currDict = dict(zip(self.keys, entry))
            varInfo = '-'.join([currDict['#CHROM'],
                                currDict['POS'],
                                currDict['REF'],
                                currDict['ALT']])
            variantIDs.append(varInfo)
        return variantIDs

def main(fname,outname):

    vcfr = VCFReader(fname = fname)

    varType = vcfr.extractType()
    allele_depths = vcfr.extractAD()
    variantIDs = vcfr.extractVariantIDs()

    exac_obj = exac(variantIDs)

    allele_freq = exac_obj.allele_freq()
    conseq = exac_obj.consequence()

    L = len(varType)
    combinedFeats = []
    for i in range(L):
        if conseq[i]:
            varType[i][0] += ' (' + conseq[i] + ')'
        combinedFeats.append(
            varType[i] + allele_depths[i] + allele_freq[i]
        )

    combinedFeats_df = pd.DataFrame(combinedFeats, columns=vcfr.annotation_header)
    combinedFeats_df.to_csv(outname, sep='\t', header=True, index=False)

if __name__ == "__main__":
    fname = sys.argv[1]
    outname = sys.argv[2]
    main(fname, outname)





    

