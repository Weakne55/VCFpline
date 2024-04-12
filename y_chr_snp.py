import io
import os
import pandas as pd
import sys
#import argparse

#parser = argparse.ArgumentParser()
#parser.add_argument('input_file', help='Path to the input file.')
#args = parser.parse_args()  

print(os.getcwd())

def read_vcf(path):
    with open(path, 'r') as f:
        lines = [l for l in f if not l.startswith('##')]
    return pd.read_csv(io.StringIO(''.join(lines)),
        dtype={'#CHROM': str, 'POS': int, 'ID': str, 'REF': str, 'ALT': str,'QUAL': str, 'FILTER': str, 'INFO': str},
        sep='\t').rename(columns={'#CHROM': 'CHROM'})

file = sys.argv[1]
data = read_vcf('chrY.vcf')
haplo = pd.read_excel('../../Haplogroups.xlsx',index_col=None).drop(columns=['Unnamed: 0'])
haplo[['REF', 'ALT']] = haplo['Mutation info'].str.split('->', expand=True)
haplo.drop(columns=['Mutation info'], inplace=True)
intersect_df = pd.merge(data, haplo, how ='inner', on =['POS', 'REF', 'ALT']) 
intersect_df = intersect_df.drop(columns=['ID','FILTER','INFO','FORMAT',file,'rs #'])

def remove_nan_from_list(lst):
    return [x for x in lst if x != 'nan']

def bt(el):
    return [x for x in el if len(x)>3]

intersect_df['SNPs'] = intersect_df['SNPs'].str.split(';').apply(lambda x: [item for item in x if item.strip() != 'nan']).apply(lambda x: ', '.join(x))

intersect_df['Haplogroup'] = intersect_df['Haplogroup'].str.replace('~', 'QQQ')

intersect_df.sort_values(by='Haplogroup', inplace=True, kind='mergesort')
intersect_df['Haplogroup'] = intersect_df['Haplogroup'].str.replace('QQQ', '~')

intersect_df.to_csv('temp.csv', sep='\t', encoding='utf-8')
