import os
import numpy as np
import pandas as pd
import re
from datetime import datetime
from sys import argv
#choose whether reads are already blasted
if len(argv) > 1:
    need_to_BLAST = bool(int(argv[1]))
else:
    need_to_BLAST = True

home = 'C:/Users/zokmi/space/'
ref_dir = 'ref/'
blast_tmp = 'blast_intermediates/'
directory = '2026_article_suppr/'
path = home + directory
sum_dir = 'data/'
print(f'{datetime.now():%d.%m.%Y %H:%M:%S} {path}')
# read reference file
ref = pd.read_csv(f'{home}{ref_dir}reference.txt')
print(f'{datetime.now():%d.%m.%Y %H:%M:%S}')
print(ref.head(5))

dmp = home + blast_tmp
if need_to_BLAST:
    if not os.path.exists(dmp): os.mkdir(dmp)
    if not os.path.exists(dmp + 'out/'): os.mkdir(dmp + 'out/')
    # make fasta file for BLAST database
    f = open(f'{dmp}ref_db.fasta', "w")
    for i in range(len(ref['UPTAG_seqs'])):
        i1 = ref['UPTAG_seqs'][i]
        i2 = i1 + '-refseqname'
        f.write(f'>{i2}\n{i1}\n')
    f.close()
    # BLAST command, making nucleotide database 'ref'
    os.system(f'makeblastdb -in {dmp}ref_db.fasta -dbtype nucl -out {dmp}ref')

sumblast = pd.DataFrame({
    'exp': [], 'all_count': [], 'm_count': [], 'percentage': [], 'qual': []
})
sumblastdm = pd.DataFrame({
    'exp': [], 'all_count': [], 'm_count': [], 'percentage': [], 'qual': []
})

files = list()
for x,y,z in os.walk(path):#os.listdir(path):  # list of directory and file names
    for i in z:
        file = os.path.join(x, i)
        if file.endswith('not_matched.csv') or file.endswith('not_matched_dm.csv'):
            files.append(file)
# dirs.append('')
filelist = '\n'.join(files)
print(f"{datetime.now():%d.%m.%Y %H:%M:%S}\n{filelist}")

# parsing not_matched files in folders
for j in files:
    p = os.path.abspath(f'{j}/..')

    print(f'{datetime.now():%d.%m.%Y %H:%M:%S} {j}') 
    nb = pd.read_csv(j)
    barcodes = nb['barcode']
    print(f'{datetime.now():%d.%m.%Y %H:%M:%S}', len(barcodes), 'barcodes')

    f = open(f'{dmp}tmp.fasta', "w")  # making fasta with BLAST queries
    for k in barcodes:
        f.write(f'>{k}-unknownseqname\n{k}\n')
    f.close()
    print(f'{datetime.now():%d.%m.%Y %H:%M:%S} {path}tmp.fasta')
    # search with blastn not_matched sequences in database ref. output file end with blastout3.txt. e-value 10, alignment initiating word size 6, search on the same strand, finds only the best hit. outfile contains different information about alignment in tsv format.
    os.system(
        f'blastn -query {dmp}tmp.fasta -db {dmp}ref -out {dmp}out/tmp_blastout.txt -evalue 10 -word_size 6 -strand plus -max_target_seqs 1 -max_hsps 1 -outfmt \"6 qacc qlen sacc slen length nident evalue qstart qend sstart send\" -num_threads 6')  # -num_threads 8
    # subtract -unknownseqname in outfile.
    with open(f'{dmp}out/tmp_blastout.txt', 'r') as f:
        old_data = f.read()
    new_data = 'qacc\tqlen\tsacc\tslen\tlength\tnident\tevalue\tqstart\tqend\tsstart\tsend\n' + old_data.replace(
        '-unknownseqname\t', '\t')
    with open(f'{dmp}out/tmp_blastout.txt', 'w') as f:
        f.write(new_data)
    with open(f'{dmp}out/tmp_blastout.txt', 'r') as f:
        old_data = f.read()
    new_data = old_data.replace('-refseqname\t', '\t')
    with open(f'{dmp}out/tmp_blastout.txt', 'w') as f:
        f.write(new_data)

    check = pd.read_csv(f'{dmp}out/tmp_blastout.txt', sep='\t')
    print(f'{datetime.now():%d.%m.%Y %H:%M:%S}')
    print(check)
    # filter distance between sequences
    check['dist'] = (check['slen'] + check['qlen'] - check['length'] - check['nident'])
    check = check[check['dist'] <= 2]
    print(len(check))
    print(f'{datetime.now():%d.%m.%Y %H:%M:%S} BLASTed')
    print(f'{datetime.now():%d.%m.%Y %H:%M:%S} adding reference and not_matched data')

    barcode = pd.DataFrame({'barcode': barcodes})
    full = pd.merge(barcode, check, left_on='barcode', right_on='qacc', how='left')
    full = pd.merge(full, ref, left_on='sacc', right_on='UPTAG_seqs', how='left')
    full = full[~pd.isna(full['qacc'])].reset_index()[['barcode', 'Confirmed_deletion', 'UPTAG_seqs', 'UPTAG_notes']]
    print(f'{datetime.now():%d.%m.%Y %H:%M:%S}')
    print(full)

    if j.endswith('not_matched.csv'):
        name = os.path.basename(j).split('_not_matched.csv')[0]
        f = pd.merge(full, nb, left_on='barcode', right_on='barcode', how='right')
        m = f[~pd.isna(f['Confirmed_deletion'])].reset_index()[
            ['barcode', 'Confirmed_deletion', 'UPTAG_seqs', 'UPTAG_notes', 'n', 'count', 'qual']]
        # print(str(datetime.now()))
        m = m.groupby('Confirmed_deletion').agg(
            barcode=('barcode', lambda x: '|'.join(x[~pd.isna(x)].unique())),
            n=('n', 'sum'),
            count=('count', 'sum'),
            notes=('UPTAG_notes', lambda x: '|'.join(x[~pd.isna(x)].unique())),
            original_barcode=('UPTAG_seqs', lambda x: '|'.join(x[~pd.isna(x)].unique()))).reset_index()

        m = m.sort_values(by=['count'], ascending=[False])
        print(f'{datetime.now():%d.%m.%Y %H:%M:%S} {p}/{name}_blasted.csv')
        m.to_csv(f'{p}/{name}_blasted.csv', index=False)

        sumblast.loc[len(sumblast.index)] = [
            name,
            round(nb['count'].sum()),
            round(m['count'].sum()),
            round(m['count'].sum()) / round(nb['count'].sum()),
            nb['qual'].mean().item()
        ]


    if j.endswith('not_matched_dm.csv'):
        name = os.path.basename(j).split('_not_matched_dm.csv')[0]
        f = pd.merge(full, nb, left_on='barcode', right_on='barcode', how='right')
        m = f[~pd.isna(f['Confirmed_deletion'])].reset_index()[
            ['barcode', 'Confirmed_deletion', 'UPTAG_seqs', 'UPTAG_notes', 'n', 'count', 'qual']]
        m = m.groupby('Confirmed_deletion').agg(
            barcode=('barcode', lambda x: '|'.join(x[~pd.isna(x)].unique())),
            n=('n', 'sum'),
            count=('count', 'sum'),
            notes=('UPTAG_notes', lambda x: '|'.join(x[~pd.isna(x)].unique())),
            original_barcode=('UPTAG_seqs', lambda x: '|'.join(x[~pd.isna(x)].unique()))).reset_index()

        m = m.sort_values(by=['count'], ascending=[False])
        print(f'{datetime.now():%d.%m.%Y %H:%M:%S} {p}/{name}_blasted_dm.csv')
        m.to_csv(f'{p}/{name}_blasted_dm.csv', index=False)

        sumblastdm.loc[len(sumblastdm.index)] = [
            name,
            round(nb['count'].sum()),
            round(m['count'].sum()),
            round(m['count'].sum()) / round(nb['count'].sum()),
            nb['qual'].mean().item()
        ]


sumblast.to_csv(f'{path}{sum_dir}sumblast.csv', index=False)
sumblastdm.to_csv(f'{path}{sum_dir}sumblastdm.csv', index=False)
