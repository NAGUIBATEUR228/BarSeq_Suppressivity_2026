import os
import pandas as pd
import numpy as np
import re
from datetime import datetime
from sys import argv

# choose wether the reads are already blasted
if len(argv) > 1:
    need_to_BLAST = bool(int(argv[1]))
else:
    need_to_BLAST = True

#script of blasting each file
#need_to_BLAST = True
# adapter sequences
u1 = "GATGTCCACGAGGTCTCT"
u2 = "CGTACGCTGCAGGTCGAC"

# reverse-complement string
def rc(seq):
    complement = {'A': 'T', 'C': 'G', 'G': 'C', 'T': 'A'}
    reverse_complement = ""
    for nt in seq:
        if nt in complement:
            reverse_complement = complement[nt] + reverse_complement
    return reverse_complement


# function for checking if there are sequences with two alignments with one primer
def ch(check):
    a = set(check.groupby(['qacc'])['sacc'].value_counts().values)
    a.add(1)
    a=set(a)
    print(a, max(a))
    return max(a) != 1

# make a hash table with names and sequences of adapters
ud = {'u1': u1, 'u2': u2, 'u1_rc': rc(u1), 'u2_rc': rc(u2)}

home = 'C:/Users/zokmi/space/'
ref_dir = 'ref/'
blast_tmp = 'blast_intermediates/'
directory = '2026_article_suppr/'
path = home + directory
sum_dir = 'data/'
print(f'{datetime.now():%d.%m.%Y %H:%M:%S} {path}')
# read reference file
ref = pd.read_csv(f'{home}{ref_dir}reference.txt')

# make intermediate output directories
dmp = home + blast_tmp
if need_to_BLAST:
    if not os.path.exists(dmp): os.mkdir(dmp)
    if not os.path.exists(dmp + 'out/'): os.mkdir(dmp + 'out/')
    f = open(f'{dmp}query_u.fasta', "w")
    for i, j in ud.items():
        f.write(f'>{i}\n{j}\n')
    f.close()
    os.system(f'makeblastdb -in {dmp}query_u.fasta -dbtype nucl -out {dmp}uref')

# make summary table
sumdm = pd.DataFrame({
    'exp': [], 'total_count': [], 'u_count': [], 'barcoded': [], 'u_barcoded': [], 'not_barcoded': [], 'u_nb': [],
    'matched': [], 'u_matched': [], 'not_matched': [], 'u_nm': []
})

files = list()
for x,y,z in os.walk(path):#os.listdir(path):  # list of directory and file names
    for i in z:
        file = os.path.join(x, i)
        if file.endswith('not_barcoded_raw_qual_count.csv'):
            files.append(file)
# dirs.append('')
filelist = '\n'.join(files)
print(f"{datetime.now():%d.%m.%Y %H:%M:%S}\n{filelist}")  # directories with fastq files

# parsing not_barcoded files in folders
for j in files:
    p = os.path.abspath(f'{j}/..')

    print(f'{datetime.now():%d.%m.%Y %H:%M:%S} {j}') 
    nb = pd.read_csv(j)  # not_barcoded table
    # print(len(nb[~nb['barcode'].isin(seqs)]))
    seqs = nb['seq']#pd.Series(pd.concat([seqs, , ignore_index=True).unique())

    print(f'{datetime.now():%d.%m.%Y %H:%M:%S}', len(seqs), 'reads')
    name = os.path.basename(j).split('_not_barcoded_raw_qual_count.csv')[0]
    print(f'{datetime.now():%d.%m.%Y %H:%M:%S} tmp_dm.fasta')
    with open(f'{dmp}tmp_dm.fasta', "w") as f:
        for k in seqs:
            f.write(f'>{k}-unknownseqname\n{k}\n')

    os.system(
        f'blastn -query {dmp}tmp_dm.fasta -db {dmp}uref -out {dmp}out/tmp_blastout.txt -evalue 0.001 -word_size 6 -strand plus -max_target_seqs 2 -outfmt \"6 qacc qlen sacc slen length nident evalue qstart qend sstart send\" -num_threads 6')  #

    print(f'{datetime.now():%d.%m.%Y %H:%M:%S} BLASTed')

    
    #subtract -unknownseqname in outfile and add a header.
    with open(f'{dmp}out/tmp_blastout.txt', 'r') as f:
        old_data = f.read()
    new_data = 'qacc\tqlen\tsacc\tslen\tlength\tnident\tevalue\tqstart\tqend\tsstart\tsend\n' + old_data.replace(
        '-unknownseqname\t', '\t')
    with open(f'{dmp}out/tmp_blastout.txt', 'w') as f:
        f.write(new_data)

    print(f'{datetime.now():%d.%m.%Y %H:%M:%S} filtering')
    # read blastn result
    check = pd.read_csv(f'{dmp}out/tmp_blastout.txt', sep='\t')
    # qacc - read sequence, sacc - primer name in 'ud'. remove sequences with only one primer found
    check = check[check.groupby(['qacc'])['sacc'].transform('count') >= 2]
    # split table in two parts. one of them contains sequences with more than one alignments with one primer. it is needed to choose only one.
    mask = check.groupby(['qacc', 'sacc'])['sacc'].transform('count') != 1
    check1 = check[~mask]
    check = check[mask]


    # pass some selection procedures of the best alignments
    changed = True
    if changed:
        if ch(check):
            check = check[(check.groupby(['qacc', 'sacc'])['length'].transform('max')) == check['length']]
        else:
            changed = False
    if changed:
        if ch(check):
            sacc = check['sacc'].values
            conditions = [
                sacc == 'u1',
                sacc == 'u2',
                sacc == 'u2_rc',
                sacc == 'u1_rc'
            ]

            values = [
                check['send'].values,
                -check['sstart'].values,
                check['send'].values,
                -check['sstart'].values
            ]
            check['tof'] = np.select(conditions, values, default=None)
            check = check[(check.groupby(['qacc', 'sacc'])['tof'].transform('max')) == check['tof']]
        else:
            changed = False

    if changed:
        if ch(check):
            check = check[(check.groupby(['qacc', 'sacc'])['nident'].transform('max')) == check['nident']]
        else:
            changed = False
    if changed:
        if ch(check):
            check = check[(check.groupby(['qacc', 'sacc'])['evalue'].transform('min')) == check['evalue']]
        else:
            changed = False
    if changed:
        if ch(check):
            sacc = check['sacc'].values
            conditions = [
                sacc == 'u1',
                sacc == 'u2',
                sacc == 'u2_rc',
                sacc == 'u1_rc'
            ]

            values = [
                check['qend'].values,
                -check['qstart'].values,
                check['qend'].values,
                -check['qstart'].values
            ]
            check['tof'] = np.select(conditions, values, default=None)
            check = check[(check.groupby(['qacc', 'sacc'])['tof'].transform('min')) == check['tof']]
        else:
            changed = False
    print(f'{datetime.now():%d.%m.%Y %H:%M:%S} {ch(check)}')
    # combine tables
    check = pd.concat([check, check1], ignore_index=True)
    sacc = check['sacc'].values
    conditions = [
        sacc == 'u1',
        sacc == 'u2',
        sacc == 'u2_rc',
        sacc == 'u1_rc'
    ]

    values = [
        check['qend'].values,
        check['qstart'].values - 1,
        check['qend'].values,
        check['qstart'].values - 1
    ]
    check['ind'] = np.select(conditions, values, default=None)

    # transform table to wide format of coordinates of alignments' borders
    print(f'{datetime.now():%d.%m.%Y %H:%M:%S} pivot')
    w = pd.pivot(data = check, values = 'ind', index = 'qacc', columns = 'sacc').reset_index()
    w = w.where((pd.notnull(w)), None)
    u = w.columns.values.tolist()
    # check for alignments that allows to extract barocodes
    if sum([x in u for x in ud.keys()]) == 4:
        mask = (~pd.isna(w['u1']) & ~pd.isna(w['u2']) & (w['u1'] <= w['u2'])) | (
                    ~pd.isna(w['u1_rc']) & ~pd.isna(w['u2_rc']) & (w['u2_rc'] <= w['u1_rc']))
    elif 'u1' in u and 'u2' in u:
        mask = (~pd.isna(w['u1']) & ~pd.isna(w['u2']) & (w['u1'] <= w['u2']))
    elif 'u1_rc' in u and 'u2_rc' in u:
        mask = (~pd.isna(w['u1_rc']) & ~pd.isna(w['u2_rc']) & (w['u2_rc'] <= w['u1_rc']))
        if 'u1' in u: w.drop(columns=['u1'], inplace=True)
        if 'u2' in u: w.drop(columns=['u2'], inplace=True)
    else:
        import sys
        print('data have no barcode')
        sys.exit()
    w = w[mask]

    # extract barcodes
    print(f'{datetime.now():%d.%m.%Y %H:%M:%S} extraction')
    v = np.vectorize(lambda x, y, z: x[y:z])
    vecrc = np.vectorize(rc)
    w['barcode'] = w['qacc']
    if 'u1' in u:
        w['barcode'] = v(w['qacc'].values, w['u1'].values, w['u2'].values)
    if 'u1_rc' in u:
        w['barcode'] = np.where(w['barcode'].values == w['qacc'].values,
                                vecrc(v(w['qacc'].values, w['u2_rc'].values, w['u1_rc'].values)),
                                w['barcode'].values)
    w.set_index('qacc', inplace=True)
    print(f'{datetime.now():%d.%m.%Y %H:%M:%S} filtered and barcoded')



    inds = nb['seq'][nb['seq'].isin(set(w.index))]
    to_merge = w.loc[inds, :].reset_index()
    print(f'{datetime.now():%d.%m.%Y %H:%M:%S} {j}')
    full = pd.merge(to_merge, nb, left_on='qacc', right_on='seq', how='outer')
    full = pd.merge(full, ref, left_on='barcode', right_on='UPTAG_seqs', how='left')
    # stats
    exess = round(full[pd.isna(full['qacc'])]['count'].sum())
    exessu = round(len(full[pd.isna(full['qacc'])]['seq'].unique()))

    full = full[~pd.isna(full['qacc'])]

    print(f'{datetime.now():%d.%m.%Y %H:%M:%S} aggregation')
    m = full[~pd.isna(full['Confirmed_deletion'])].copy()
    m['count'] = m.groupby('barcode')['count'].transform('sum')
    m['n'] = m.groupby('barcode')['count'].transform('count')
    m['notes'] = m['UPTAG_notes']
    m = m[['Confirmed_deletion', 'barcode', 'n', 'count', 'notes']].drop_duplicates()
    m = m.groupby('Confirmed_deletion').agg(
        barcode=('barcode', lambda x: '|'.join(x[~pd.isna(x)].unique())),
        n=('n', 'sum'),
        count=('count', 'sum'),
        notes=('notes', lambda x: '|'.join(x[~pd.isna(x)].unique()))
    ).reset_index()

    print(f'{datetime.now():%d.%m.%Y %H:%M:%S} {p}/{name}_output_dm_count.csv')
    m = m.sort_values(by=['count'], ascending=[False])
    m.to_csv(f'{p}/{name}_output_dm_count.csv', index=False)

    print(f'{datetime.now():%d.%m.%Y %H:%M:%S} nm')
    nm = full[pd.isna(full['Confirmed_deletion'])].copy()
    nm['qual'] = nm['qual'] * nm['count']
    nm['count'] = nm.groupby('barcode')['count'].transform('sum')
    nm['n'] = nm.groupby('barcode')['count'].transform('count')
    nm['qual'] = nm.groupby('barcode')['qual'].transform('sum')
    nm['qual'] = nm['qual'] / nm['count']
    nm = nm[['barcode', 'n', 'count', 'qual']].drop_duplicates()
    print(nm['qual'].mean())
    print(nb['qual'].mean())

    print(f'{datetime.now():%d.%m.%Y %H:%M:%S} {p}/{name}_not_matched_dm_count.csv')
    nm.to_csv(f'{p}/{name}_not_matched_dm.csv', index=False)

    sumdm.loc[len(sumdm.index)] = [
        name,
        round(nb['count'].sum()),
        round(len(nb['seq'].unique())),
        round(full['count'].sum()),
        round(len(full['barcode'].unique())),
        exess,
        exessu,
        round(m['count'].sum()),
        round(len(m['barcode'].unique())),
        round(nm['count'].sum()),
        round(len(nm['barcode'].unique()))
    ]


# sum information in percent
sumdm['perc_barcoded'] = sumdm['barcoded'] / sumdm['total_count']
sumdm['perc_matched'] = sumdm['matched'] / sumdm['total_count']
sumdm.to_csv(f'{path}{sum_dir}sumdm.csv', index=False)
print(f'{datetime.now():%d.%m.%Y %H:%M:%S} DONE')
