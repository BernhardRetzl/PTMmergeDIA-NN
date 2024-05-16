import pandas as pd
from Bio import SeqIO
import re
import numpy as np
import glob


#--------------------------------------------------------------------------------------------------------
#--------------------------------------------------------------------------------------------------------
# searches a FASTA-file for proteins listed in the Protein.Id table
df = pd.read_csv('report.pr_matrix.tsv', sep='\t')
# df = pd.read_excel('test.xlsx')

df = df[df['Modified.Sequence'].str.contains('Mod', na=False)]


# df = pd.read_excel('test.xlsx')
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
    # if len(parts) == 1 or 0:
    #     return
    # else:
    length = np.array([len(i) for i in parts])
    length = length.cumsum()+position
    length = [str(i) for i in length.tolist()]

    aa_of_interest = [i[-1] for i in parts]

    uni_mod_matches = re.findall(uni_mod_pattern, string)



    protein_ids = [protein_ids] * len(length)
    to_return = list(zip(protein_ids, length, aa_of_interest, uni_mod_matches))
    to_return = ['-'.join(i) for i in to_return]
    print(to_return)
    return to_return

    return '-'.join((protein_ids, '-'.join(length), '-'.join(aa_of_interest), '-'.join(uni_mod_matches)))


corrected_df = pd.DataFrame()



# Assuming return_length is a function that returns a list of items for grouping.
def optimized_processing(df):
    # Process each row and return a new DataFrame row as a dictionary in the list
    rows = []
    for index, row in df.iterrows():
        modified_sequence = row['Modified.Sequence']
        stripped_seq_pos = row['stripped_seq_pos']
        protein_ids = row['Protein.Ids']

        for_grouping = return_length(modified_sequence, stripped_seq_pos, protein_ids)
        for item in for_grouping:
            # Create a dictionary for the new row, modifying it as needed
            new_row = row.to_dict()
            new_row['for_grouping'] = item
            rows.append(new_row)

    # Create DataFrame from list of dictionaries
    corrected_df = pd.DataFrame(rows)
    return corrected_df

# Usage
corrected_df = optimized_processing(df)

grouped_df = corrected_df.groupby(['for_grouping']).sum()
grouped_df = grouped_df.iloc[:, 10:-3]
corrected_df.set_index('for_grouping', inplace=True)
corrected_df = corrected_df.iloc[:,:10 ]
corrected_df = corrected_df.groupby(['for_grouping']).first()


joined_df = corrected_df.join(grouped_df)
joined_df.to_excel('joined_df.xlsx')
