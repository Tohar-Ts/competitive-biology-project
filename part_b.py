from cmath import nan
from tkinter.messagebox import NO
from Bio import SeqIO
import numpy as np
import matplotlib.pyplot as plt
import pandas as pd
from part_a import *
import io   

UNIPORT_PATH = "data/BS168.xlsx"

def open_xlsx_file():
    return pd.read_excel(UNIPORT_PATH, sheet_name='Sheet0')

# ------------B.1--------------
def clean_id(id_arr):
    to_delete = ["_", ";"]
    to_replace = ["/"]
    clean_arr = []
    
    for id in id_arr:
        if str(id) == 'nan':
            continue
        
        for d in to_delete:
            id = str(id).replace(d, "")
            
        for r in to_replace:
            id = str(id).replace(r, " ")

        mult_id = str(id).split()
        for i in mult_id:
            clean_arr.append(i)
    
    return clean_arr


def find_missing(first, second):
    """
    Function for finding elements which
    are there in first[] but not in second[].
    """
    missing = []
    for f in first:
        count = 0
        
        for j in second:
            if (f == j):
                count +=1
                break
        
        if count == 0:
            missing.append(f)
   
    return missing


def compare_gb_to_up(cds_id, up_df):
    clean_cds_id = clean_id(cds_id)
    clean_up_df =  clean_id(up_df)

    compare_gb_to_up = find_missing(clean_cds_id, clean_up_df) 
    compare_up_to_gb = find_missing(clean_up_df, clean_cds_id)
   
    gb_total = len(clean_cds_id)
    gb_missing_len = len(compare_gb_to_up)
    gb_not_missing_len = gb_total - gb_missing_len
    
    up_total = len(clean_up_df)
    up_missing_len = len(compare_up_to_gb)
    up_not_missing_len = up_total - up_missing_len

    up_nan = [id for id in up_df if id == 'nan']
    up_nan_len = up_total - up_missing_len
    
    print("In Genebank but not in Uniport: \n", compare_gb_to_up)
    print("Missing: ", gb_missing_len, " from total:", gb_total)
    print("\nIn Uniport but not in Genebank (after removing nan): \n", compare_up_to_gb)
    print("Missing: ", up_missing_len, " from total:", up_total)
    
    plt.figure(figsize=(10, 6))

    plt.subplot(1, 2, 1)
    plt.title("GeneBank Genes")
    plt.pie([gb_missing_len, gb_not_missing_len], 
            labels = ['not in uniport', 'in uniport'], explode = [0.2, 0.1], 
            shadow = True, autopct='%1.3f%%')

    plt.subplot(1, 2, 2)
    plt.title("Uniport Genes")
    plt.pie([up_missing_len, up_nan_len, up_not_missing_len], 
            labels =  ['not in genebank', 'nan values', 'in genebank'], explode = [0.2, 0.1, 0.1],
            shadow = True, autopct='%1.3f%%')
       
    plt.suptitle("Comparing GeneBank and Uniport Proteins")
    plt.tight_layout()
    plt.show()

# ------------B.2--------------
def create_trans_table(up_df):
    up_clean_df = up_df.dropna(subset=['Transmembrane'])
    
    start_list, end_list, len_list, sequence_list, id_list = [], [], [], [], []
    
    for __, row in up_clean_df.iterrows():
        trans = row['Transmembrane']    
        sequence = row['Sequence']
        id = row['Gene ID']
        
        trans = trans.replace('TRANSMEM', '').replace(' ', '').replace('..', '-')
        trans = trans.split(';')
        trans_indexes = [t for t in trans if '\"' not in t]
        
        for i in trans_indexes:
            index = i.split('-')

            start = int(index[0])
            end = int(index[1])
            seq = sequence[start:end + 1]
            length = len(seq)
            
            start_list.append(start)
            end_list.append(end)
            sequence_list.append(seq)
            len_list.append(length)
            id_list.append(id)
    
    # creating the data frame
    df = pd.DataFrame(zip(id_list, start_list, end_list, sequence_list, len_list), 
                      columns=['id', 'start', 'end', 'sequence', 'length'])
    return df        
        
        
def plot_trans_len_stat(trans_len_arr):
    print("Transmembrane Lengths stats:")
    stat(trans_len_arr)
    
    plt.title("Transmembrane Lengths")
    plt.hist(trans_len_arr)
    plt.tight_layout()
    plt.show() 
    
    
if __name__ == "__main__":
    rec, gb_df = open_and_create_gb_dataframe()
    up_df = open_xlsx_file()

    # Q1
    # print("\n------------Question 1------------")
    # cds, __ = group_genes(gb_df)       # group to CDS and others
    # cds_id = np.asarray(cds['id'])           # converting id col to arr
    # up_df = np.asarray(up_df['Gene ID'])
    # compare_gb_to_up(cds_id, up_df)            # printing the comparing results
 
    # Q2 
    print("\n------------Question 2------------")
    trans_df = create_trans_table(up_df)
    trans_len_arr = np.asarray(trans_df['length'])        # get the arr of len for each group
    plot_trans_len_stat(trans_len_arr)
    
    # Q3
    # print("\n------------Question 3------------")
  