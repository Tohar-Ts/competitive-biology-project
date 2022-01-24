from cmath import nan
from Bio import SeqIO
import numpy as np
import matplotlib.pyplot as plt
import pandas as pd
from part_a import *

UNIPORT_PATH = "data/BS168.xlsx"

def open_xlsx_file():
    return pd.read_excel(UNIPORT_PATH, sheet_name='Sheet0')

# ------------B.1--------------
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
    clean_cds_id = [id.replace("_", "") for id in cds_id]
    clean_up_df = [str(id).replace("_", "") for id in up_df if not pd.isna(id)]

    compare_gb_to_up = find_missing(clean_cds_id, clean_up_df)
    compare_up_to_gb = find_missing(clean_up_df, clean_cds_id)
    
    total_gb = len(clean_cds_id)
    missing_gb_len = len(compare_gb_to_up)
    
    total_up = len(clean_up_df)
    missing_up_len = len(compare_up_to_gb)
    
    print("In Genebank but not in Uniport: \n", compare_gb_to_up)
    print("Missing: ", missing_gb_len, " from total:", total_gb)
    print("\nIn Uniport but not in Genebank: \n", compare_up_to_gb)
    print("Total: ", missing_up_len, " from total:", total_up)
    
    labels = ['not in', 'in']
    plt.figure(figsize=(10, 6))

    plt.subplot(1, 2, 1)
    plt.title("GeneBank Genes")
    plt.pie([missing_gb_len, total_gb-missing_gb_len], labels = labels)

    plt.subplot(1, 2, 2)
    plt.title("Uniport Genes")
    plt.pie([missing_up_len, total_up-missing_up_len], labels = labels)
       
    plt.suptitle("Comparing GeneBank and Uniport Proteins")
    plt.tight_layout()
    plt.show()
    
    
if __name__ == "__main__":
    rec, gb_df = open_and_create_gb_dataframe()
    up_df = open_xlsx_file()

    # Q1
    print("\n------------Question 1------------")
    cds, other_gene = group_genes(gb_df)       # group to CDS and others
    
    cds_id = np.asarray(gb_df['id'])
    up_df = np.asarray(up_df['Gene ID'])
    
    compare_gb_to_up(cds_id, up_df)
 
    # Q2 
    print("\n------------Question 2------------")

    # Q3
    print("\n------------Question 3------------")
  