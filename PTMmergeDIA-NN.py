import pandas as pd
from Bio import SeqIO
import re
import numpy as np
import glob


#--------------------------------------------------------------------------------------------------------
#--------------------------------------------------------------------------------------------------------
# searches a FASTA-file for proteins listed in the Protein.Id table
# df = pd.read_csv('report.pr_matrix.tsv', sep='\t')
df = pd.read_excel('test.xlsx')
uni_prot_proteins = list(df['Protein.Ids'])
uni_prot_proteins = {i: '' for i in uni_prot_proteins}

fasta_file = glob.glob('*.fasta')[0]
fasta_file = SeqIO.parse(fasta_file, 'fasta')
for sequence in fasta_file:
    sequence_id = sequence.id.split('|')[1]
    if sequence_id in uni_prot_proteins.keys():
        uni_prot_proteins[sequence_id] = str(sequence.seq)

df['found_protein'] = uni_prot_proteins
df['found_protein'] = df['Protein.Ids'].apply(lambda x: uni_prot_proteins[x])

for index, row in df.iterrows():
    df['stripped_seq_pos'] = df.loc[index, 'found_protein'].find(df.loc[index, 'Stripped.Sequence'])



#---------------------------------------------------------------------------------------------------------
#---------------------------------------------------------------------------------------------------------
# finds the positions of strings starting with ( and endint with ) in a string
def return_length(string, position, protein_ids):
    # pattern for splitting based on e.g.: (UniMod:21)
    pattern = r"\([^)]*\)"
    # pattern for finding the specific UniMod (e.g. 21)
    uni_mod_pattern = r"(\d+)"
    parts = re.split(pattern, string)
    parts = [i for i in parts if i]
    if len(parts) == 1 or 0:
        return
    else:
        length = np.array([len(i) for i in parts][0:-1])
        length = length.cumsum()+position
        length = [str(i) for i in length.tolist()]

        aa_of_interest = [i[-1] for i in parts][0:-1]

        uni_mod_matches = re.findall(uni_mod_pattern, string)

        return '-'.join((protein_ids, '-'.join(length), '-'.join(aa_of_interest), '-'.join(uni_mod_matches)))


for index, row in df.iterrows():
    modified_sequence = df.loc[index, 'Modified.Sequence']
    stripped_seq_pos = df.loc[index, 'stripped_seq_pos']
    protein_ids = df.loc[index, 'Protein.Ids']

    df.loc[index, 'for_grouping'] = return_length(modified_sequence, stripped_seq_pos, protein_ids)


#---------------------------------------------------------------------------------------------------------


grouped_df = df.groupby(['for_grouping']).sum()
grouped_df = grouped_df.iloc[:, 10:-3]
df.set_index('for_grouping', inplace=True)
df = df.iloc[:,:10 ]


def filter_longest_strings(group):
    max_length = group['Modified.Sequence'].str.len().max()
    return group[group['Modified.Sequence'].str.len() == max_length]


df = df.groupby(level=0).apply(filter_longest_strings)
joined_df = df.join(grouped_df)
joined_df.to_excel('joined_df.xlsx', index=False)
