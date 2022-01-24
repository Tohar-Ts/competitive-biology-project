from Bio import SeqIO
import numpy as np
import matplotlib.pyplot as plt
import pandas as pd

UNIPORT_PATH = "data/BS168.gb"


def open_and_create_gb_dataframe():
    with open(UNIPORT_PATH, "r") as file:
        gen = SeqIO.parse(file, "genbank")
        rec = next(gen)  # content of 1st record
        
    id, start, end, type, strand = [], [], [], [], []
    # only for cds
    table, translation, codon_start = [], [], []
    
    feats = rec.features[1:]
    for feat in feats:
        if feat.type == 'gene':
            continue
        
        id.append(feat.qualifiers["locus_tag"][0])
        start.append(feat.location.start.position)
        end.append(feat.location.end.position)
        type.append(feat.type)
        strand.append(feat.location.strand)
        
        if feat.type == 'CDS':
            table.append(feat.qualifiers["transl_table"][0])
            translation.append(feat.qualifiers["translation"][0])
            codon_start.append(feat.qualifiers["codon_start"][0])
        else:
            table.append(None)
            translation.append(None)
            codon_start.append(None)
    
    # creating the data frame
    df = pd.DataFrame(zip(id, start, end, strand, type, table, translation, codon_start), 
                      columns=['id', 'start', 'end', 'strand', 'type', 'table', 'translation', 'codon_start'])
    return rec, df    

# ------------A.1--------------
def create_type_dict(df):
    return df.type.value_counts().to_dict()

# ------------A.2.a------------
def singel_gene_length(row):
    return np.abs(row['end'] - row['start']) + 1


def add_len_col(df):
    new_col = df.apply(lambda row: singel_gene_length(row), axis=1)
    df['len'] = new_col
    return df

# ------------A.2.b------------
def group_genes(df):
    cds = df[df['type'] == 'CDS']
    other_gene = df[df['type'] != 'CDS']
    return cds, other_gene

# ------------A.2.c------------
def stat(arr):
    min = np.min(arr)
    max = np.max(arr)
    avg = np.average(arr)
    print("average: {:.2f},".format(avg), "minimum: {:.2f},".format(min), "maximum: {:.2f}".format(max))

# ------------A.2.d------------
def plot_len_stat(cds_len, other_gene_len):
    stat(cds_len)
    stat(other_gene_len)

    all_genes = np.concatenate((cds_len,other_gene_len))
    plt.figure(figsize=(10, 6))

    plt.subplot(1, 3, 1)
    plt.title("All Genes")
    plt.hist(all_genes)

    plt.subplot(1, 3, 2)
    plt.title("CDS Genes")
    plt.hist(cds_len)

    plt.subplot(1, 3, 3)
    plt.title("Others")
    plt.hist(other_gene_len)
       
    plt.suptitle("Lengths Histograms")
    plt.tight_layout()
    plt.show()
    
# ------------A.3.a------------
def GC_percentage(sequence):
    G_count = sequence.count('G')
    C_count = sequence.count('C')
    return (G_count + C_count) / len(sequence)


def source_GC(rec):
    genome = rec.seq.upper()
    GC_percent = GC_percentage(genome)
    print("Geonome GC percent: {:.2f}%".format(GC_percent))
    return genome

# ------------A.3.b------------
def GC_percent_for_row(row, genome):
    if row['type'] != 'CDS':
        return None
    return GC_percentage(genome[row['start']:row['end'] + 1])


def add_GC_percent_col(df, genome):
    new_col = df.apply(lambda row: GC_percent_for_row(row, genome), axis=1)
    df['GC percent'] = new_col
    return df

# ------------A.3.c------------
def avg_GC_percent_col(df):
    print("Proteins average GC percent: {:.2f}%".format(df["GC percent"].mean()))

# ------------A.3.d------------
def plot_GC_stat(GC_percent_arr):
    plt.title("GC Percent Histogram")
    plt.hist(GC_percent_arr)
    plt.tight_layout()
    plt.show()

# ------------A.3.d------------
def extreme_GC_percents_genes(df, n = 5):
    top = df.nlargest(n, 'GC percent')
    bottom = df.nsmallest(n, 'GC percent')
    print("Top ", n, " GC percents genes details: ")
    print(top)
    print("Bottom ", n, " GC percents genes details: ")
    print(bottom)


if __name__ == "__main__":
    rec, gb_df = open_and_create_gb_dataframe()
    
    # Q1
    print("\n------------Question 1------------")
    t_dict = create_type_dict(gb_df)
    print("Result:", t_dict)
    
    # Q2 
    print("\n------------Question 2------------")
    gb_df = add_len_col(gb_df)                    # calculate the len
    cds, other_gene = group_genes(gb_df)       # group to CDS and others
    cds_arr = np.asarray(cds['len'])        # get the arr of len for each group
    other_gene_arr = np.asarray(other_gene['len'])
    plot_len_stat(cds_arr, other_gene_arr)  # print the stat and plot the histogram

    # Q3
    print("\n------------Question 3------------")
    genome = source_GC(rec)                 # print GC percent of the genome 
    gb_df = add_GC_percent_col(gb_df, genome)     # create GC percent column
    avg_GC_percent_col(gb_df)                  # print average GC percent of proteins
    GC_percent_arr = np.asarray(gb_df['GC percent'])
    plot_GC_stat(GC_percent_arr)            # plot GC histogram
    extreme_GC_percents_genes(gb_df)
    gb_df.to_csv('data/part_a.csv')                 # save final results to csv file
