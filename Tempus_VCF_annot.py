import argparse
import os
import pandas as pd
from collections import *
import requests


####PARSER
parser = argparse.ArgumentParser(description="Annotate VCF Tempus Challenge")
parser.add_argument("--input", type=str, action='store')
parser.add_argument("--output", type=str, action='store')
args = parser.parse_args()

def split_VCF(input_VCF, row):
    
#   transform INFO into a dictionary
    dict_pairs=row.split(';')
    info={}
    for pair in dict_pairs:
        split = pair.split('=')
        info[split[0]] = split[1]
    return info

def variant_rank(input_VCF, variants):
#   Deleterious Variant Ranking: Complex > Insertions and Deletions > Substitutions, etc.
    for variant in variants:
        if variant == 'complex':
            return variant
        elif variant == 'del' or variant == 'ins':
            return variant
        else:
            return variant

def parse_vcf_columns(input_VCF):
    variant_type2=[]
    seq_coverage=[]
    AO_list=[]
    perc_support_variant=[]
    #INFO column
    #print(input_VCF['INFO'])
    for row in input_VCF['INFO']:
        split = split_VCF(input_VCF, row)
        variant_type=split['TYPE'].split(',')
        
        #Check if more than one, if yes, select most deleterious
        if len(variant_type) > 1:
            variant_type = variant_rank(input_VCF,variant_type)
        else:
            variant_type = variant_type[0]
            
        #get TYPE
        type_annotation = ''
        if variant_type2 == 'snp' or variant_type == 'mnp':
            type_annotation = 'Substitution'
        elif variant_type2 == 'del':
            type_annotation = 'Deletion'
        elif variant_type2 == 'ins':
            type_annotation = 'Insertion'
        elif variant_type2 =='complex':
            type_annotation = 'Complex'
        else:
            type_annotation = variant_type
        variant_type2.append(type_annotation)
        
        #sequence coverage at variant
        seq_coverage.append(split['DPB'])
        AO_param = split['AO'].split(',')
        RO_param = int(split['RO'])
        AO_sum = sum(map(int, AO_param))
        AO_list.append(AO_sum)
        
        #Percentage of reads supporting the variant vs reference reads
        #Using try except because there is a divide by 0 error
        try:
            perc_calc = float(AO_sum / RO_param)
        except ZeroDivisionError:
            perc_calc = 1
        perc_support_variant.append(perc_calc)

    returned_df = pd.DataFrame(OrderedDict({'CHROM': input_VCF['CHROM'],
                                              'POS': input_VCF['POS'],
                                              'VARIANT': input_VCF['ALT'],
                                              'REF': input_VCF['REF'],
                                              'TYPE': variant_type2,
                                              'DBP': seq_coverage,
                                              'AO': AO_list,
                                              'AO/RO': perc_support_variant}))
    return returned_df

def get_Exac_data(input_VCF):
#       Obtain ExAC allele frequency, ensembl_id, and consequence using
#       ExAC REST API. 
#       
#       A example of input varID is "14-21853913-T-C", which must consist 
#       Chromosome number, Position, Reference, and Alternative, separated
#       by "-".
#
#       This function outputs ALLELE_FREQ, ensembl_id, and CONSEQUENCE
#       corresponding to allele frequency, ensemble_id, genotypic consequences from this variantion from ExAC.

    key_values = input_VCF[['CHROM', 'POS', 'REF', 'ALT']].values
    keys = []
    for row in key_values:
        row = [str(x) for x in row]
        keys.append('-'.join(row))
    allele_frequency = []
    ensembl_id = []
    consequences = []
    # extracting desired information and fetching from ExAC database through API
    for key in keys:
        ExAC_url = ('http://exac.hms.harvard.edu/rest/variant/' + key)
        exac_return = requests.get(url=ExAC_url).json()
        if len(exac_return) != 0:
            parsed_output = parse_Exac_data(exac_return)
            allele_frequency.append(parsed_output['FREQ'])
            ensembl_id.append(parsed_output['GENES'])
            consequences.append(parsed_output['CONSEQUENCE'])
        else:
            # Append missing values.
            allele_freq.append('null')
            ensembl_id.append('null')
            consequences.append('null')
        #print(exac_return)
    
    ExAC_API_df = pd.DataFrame({'ALLELE_FREQ': allele_frequency,
                                    'GENE': ensembl_id,
                                    'CONSEQUENCE': consequences})
    return ExAC_API_df

def parse_Exac_data(exac_return):
    parsed_exac_values = {}
    try:
        variant = exac_return['variant']
    except KeyError:
        parsed_exac_values['FREQ'] = 'null'
        parsed_exac_values['GENES'] = 'null'
        return parsed_exac_values
    try:
        parsed_exac_values['FREQ'] = variant['allele_freq']
    except KeyError:
        parsed_exac_values['FREQ'] = 'null'
    try:
        parsed_exac_values['GENES'] = ';'.join(variant['ensembl_id'])
    except KeyError:
        parsed_exac_values['GENES'] = 'null'
    try:
        consequence = exac_return['consequence']
        if consequence is not None:
            parsed_exac_values['CONSEQUENCE'] = ';'.join(consequence.keys())
        else:
            parsed_exac_values['CONSEQUENCE'] = 'null'
    except KeyError:
        parsed_exac_values['CONSEQUENCE'] = 'null'
    # print(parsed_exac_values)
    return parsed_exac_values

def main():

    #print(args.input)
    #print(args.output)
    
    #read in VCF file
    #Metadata information of VCF File
    with open(args.input) as vcf_file:
        #lines starting with #
        vcf_pdata_lines = [line.rstrip('\n') for line in vcf_file if line.startswith('#')]
    vcf_pdata_lines_length = len(vcf_pdata_lines)
    #Get Headers
    vcf_pdata_header = vcf_pdata_lines[-1][1:].split('\t')   

    ####Create pandas dataframe for parsing using vcf_pdata from above        
    vcf_df = pd.read_csv(args.input,sep='\t',names=vcf_pdata_header,skiprows=vcf_pdata_lines_length)
    
    #print(vcf_pdata_linesLength)
    #print(vcf_df.head())
    print("Parsing VCF Information")
    outTable = parse_vcf_columns(vcf_df) 
    print("Obtaining Variant Information using ExAC REST API")
    ExAC_API_output = get_Exac_data(vcf_df)
    print("Writing To CSV File In " + args.output)    
    outTable = pd.concat([outTable, ExAC_API_output], axis=1)
    outFile = args.output + '/annotated_' + args.input
    outTable.to_csv(path_or_buf=outFile,na_rep='null',index=False)
    print("Done")
    
if __name__ == '__main__':
    main()