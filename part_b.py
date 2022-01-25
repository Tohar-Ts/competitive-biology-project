import numpy as np
import matplotlib.pyplot as plt
import pandas as pd
from part_a import *

UNIPORT_PATH = "data/BS168.xlsx"
HYDROPHOBIC_AMINO = ['A', 'F', 'L', 'I', 'V', 'M', 'P', 'W']
    
def open_xlsx_file():
    return pd.read_excel(UNIPORT_PATH, sheet_name='Sheet0')

def open_part_a_result_file():
    return pd.read_csv(RESULT_PATH)

# ------------B.1--------------
def clean_id(id_arr, remove_nan = True):
    to_delete = ["_", ";"]
    to_replace = ["/"]
    clean_arr = []
    
    for id in id_arr:
        if str(id) == 'nan':
            if remove_nan:
                continue
            else:
                clean_arr.append('nan')
        
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


def compare_gb_to_up(cds_id, up_id, clean_cds_id, clean_up_id):
    compare_gb_to_up = find_missing(clean_cds_id, clean_up_id) 
    compare_up_to_gb = find_missing(clean_up_id, clean_cds_id)
   
    gb_total_len = len(cds_id)
    gb_missing_len = len(compare_gb_to_up)
    gb_not_missing_len = gb_total_len - gb_missing_len
    
    up_total_len = len(up_id)
    up_missing_len = len(compare_up_to_gb) 
    up_nan_len = sum([1 for id in up_id if str(id) == 'nan'])
    up_not_missing_len = up_total_len - up_missing_len - up_nan_len

    print("In Genebank but not in Uniport: \n", compare_gb_to_up)
    print("Missing:", gb_missing_len, ", from total:", gb_total_len)
    print("\nIn Uniport but not in Genebank (after removing nan): \n", compare_up_to_gb)
    print("Missing:", up_missing_len, ", Nan values: ", up_nan_len ,", from total:", up_total_len)
    
    plt.figure(figsize=(10, 6))

    plt.subplot(1, 2, 1)
    plt.title("GeneBank Genes")
    plt.pie([gb_missing_len, gb_not_missing_len], 
            labels = ['not in uniport', 'in uniport'], explode = [0.5, 0.2], 
            shadow = True, autopct='%1.3f%%')

    plt.subplot(1, 2, 2)
    plt.title("Uniport Genes")
    plt.pie([up_missing_len, up_not_missing_len, up_nan_len], 
            labels =  ['not in genebank', 'in genebank' ,'nan values'], explode = [0.5, 0.1, 0.1],
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
            length = singel_gene_length(start, end)
            
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
    plot_hist("Transmembrane Lengths", trans_len_arr, "length", [0, 45], [0, 8000])
    plt.tight_layout()
    plt.show() 
    

def plot_trans_amino_stat(trans_seq_arr):
    print("Transmembrane Hydrophobic Amino stats:")
    
    # counting how many Hydrophobic amino are in every sequence
    trans_seq_percent = []
    for seq in trans_seq_arr:
        count = sum([1 for a in seq if a in HYDROPHOBIC_AMINO])
        hydro_percent = (count / len(seq)) * 100
        trans_seq_percent.append(hydro_percent)
    
    stat(trans_seq_percent)
    
    plot_hist("Transmembrane Hydrophobic Amino Percent", trans_seq_percent, "percent"
              ,[0, 100], [0, 4000])
    plt.tight_layout()
    plt.show() 
    
# ------------B.3--------------
def group_cds_by_trans(cds, clean_cds_id, clean_trans_id, clean_up_id):
    with_mask = [(id in clean_trans_id) for id in clean_cds_id]
    without_mask = [(not (id in clean_trans_id) and (id in clean_up_id))
                    for id in clean_cds_id] 
    
    cds_with_trans = cds[with_mask]
    cds_without_trans = cds[without_mask] 
    
    return cds_with_trans, cds_without_trans
    

def plot_cds_trans_gc_percent(cds, cds_with_trans, cds_without_trans):
    all_GC_percent = np.asarray(cds['GC percent'])
    with_GC_percent = np.asarray(cds_with_trans['GC percent'])
    without_GC_percent = np.asarray(cds_without_trans['GC percent'])
    
    print("All proteins gc percent stats:")
    stat(all_GC_percent)
    print("Proteins with transmembrane gc percent stats:")
    stat(with_GC_percent)
    print("Proteins without transmembrane gc percent stats:")
    stat(without_GC_percent)

    x_title = "GC percent"
    x_range = [0, 100]
    y_range = [0, 1500]
    plt.figure(figsize=(10, 6))

    plt.subplot(1, 3, 1)
    plot_hist("All Proteins", all_GC_percent, x_title, x_range, y_range)

    plt.subplot(1, 3, 2)
    plot_hist("Proteins with transmembrane", with_GC_percent, x_title, x_range, y_range)

    plt.subplot(1, 3, 3)
    plot_hist("Proteins without transmembrane", without_GC_percent, x_title, x_range, y_range)
       
    plt.suptitle("GC Percent Histograms")
    plt.tight_layout()
    plt.show()
      
       
if __name__ == "__main__":
    gb_df = open_part_a_result_file()
    up_df = open_xlsx_file()

    # Q1
    print("\n------------Question 1------------")
    cds_df, __ = group_genes(gb_df)             # group to CDS and others
    cds_id = np.asarray(cds_df['id'])           # converting id col to arr
    up_id = np.asarray(up_df['Gene ID'])
    clean_cds_id = clean_id(cds_id)
    clean_up_id =  clean_id(up_id)
    compare_gb_to_up(cds_id, up_id, clean_cds_id, clean_up_id)          # printing the comparing results
 
    # Q2 
    print("\n------------Question 2------------")
    trans_df = create_trans_table(up_df)      # creating a table of all the trans parts
    trans_len_arr = np.asarray(trans_df['length']) 
    trans_seq_arr = np.asarray(trans_df['sequence']) 
    plot_trans_len_stat(trans_len_arr)        # printing trans lengths stats
    plot_trans_amino_stat(trans_seq_arr)      # printing trans hydro amino stats
    
    # Q3
    print("\n------------Question 3------------")
    trans_id = np.asarray(trans_df['id'])
    clean_trans_id = clean_id(trans_id)
    cds_with_trans, cds_without_trans = group_cds_by_trans(cds_df, clean_cds_id, clean_trans_id, clean_up_id)
    plot_cds_trans_gc_percent(cds_df, cds_with_trans, cds_without_trans)