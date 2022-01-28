from part_a import *

UNIPORT_PATH = os.path.join(DATA_PATH, "BS168.xlsx")
HYDROPHOBIC_AMINO = ['A', 'F', 'L', 'I', 'V', 'M', 'P', 'W']


def open_xlsx_file():
    return pd.read_excel(UNIPORT_PATH, sheet_name='Sheet0')


def open_part_a_result_file():
    return pd.read_csv(RESULT_PATH)

# ------------B.1--------------


def clean_id(id_arr, remove_nan=True):
    to_delete = ["_", ";"]
    to_replace = ["/"]
    clean_arr = []

    for name in id_arr:
        if str(name) == 'nan':
            if remove_nan:
                continue
            else:
                clean_arr.append('nan')

        for d in to_delete:
            name = str(name).replace(d, "")

        for r in to_replace:
            name = str(name).replace(r, " ")

        mult_id = str(name).split()
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
            if f == j:
                count += 1
                break

        if count == 0:
            missing.append(f)

    return missing


def show_compare_gb_to_up(cds_id_col, up_id_col):
    cds_id_arr = np.asarray(cds_id_col)
    up_id_arr = np.asarray(up_id_col)

    clean_cds_id = clean_id(cds_id_arr)
    clean_up_id = clean_id(up_id_arr)

    compare_gb_to_up = find_missing(clean_cds_id, clean_up_id)
    compare_up_to_gb = find_missing(clean_up_id, clean_cds_id)

    gb_total_len = len(cds_id_arr)
    gb_missing_len = len(compare_gb_to_up)
    gb_not_missing_len = gb_total_len - gb_missing_len

    up_total_len = len(up_id_arr)
    up_missing_len = len(compare_up_to_gb)
    up_nan_len = sum([1 for up_id in up_id_arr if str(up_id) == 'nan'])
    up_not_missing_len = up_total_len - up_missing_len - up_nan_len

    print(f"\nIn Genbank but not in Uniport: \n{compare_gb_to_up}")
    print(f"\nMissing:{gb_missing_len}, from total: {gb_total_len}")
    print(f"\nIn Uniport but not in Genbank (after removing nan): \n{compare_up_to_gb}")
    print(f"\nMissing: {up_missing_len}, Nan values: {up_nan_len}, from total: {up_total_len}")

    plt.figure(figsize=(10, 6))

    plt.subplot(1, 2, 1)
    plt.title("GenBank Genes")
    plt.pie([gb_missing_len, gb_not_missing_len],
            labels=['not in uniport', 'in uniport'], explode=[0.5, 0.2],
            shadow=True, autopct='%1.3f%%')

    plt.subplot(1, 2, 2)
    plt.title("Uniport Genes")
    plt.pie([up_missing_len, up_not_missing_len, up_nan_len],
            labels=['not in genbank', 'in genbank', 'nan values'], explode=[0.5, 0.1, 0.1],
            shadow=True, autopct='%1.3f%%')

    plt.suptitle("Comparing Genbank and Uniport Proteins")
    plt.tight_layout()
    plt.show()

# ------------B.2--------------


def create_trans_table(df):
    up_clean_df = df.dropna(subset=['Transmembrane'])

    start_list, end_list, len_list, sequence_list, id_list = [], [], [], [], []

    for __, row in up_clean_df.iterrows():
        trans = row['Transmembrane']
        sequence = row['Sequence']
        gene_id = row['Gene ID']

        trans = trans.replace('TRANSMEM', '').replace(' ', '').replace('..', '-')
        trans = trans.split(';')
        trans_indexes = [t for t in trans if '\"' not in t]

        for index_lst in trans_indexes:
            index = index_lst.split('-')

            start = int(index[0])
            end = int(index[1])
            sub_seq = sub_sequence(sequence, start, end)
            length = single_gene_length(start, end)

            start_list.append(start)
            end_list.append(end)
            sequence_list.append(sub_seq)
            len_list.append(length)
            id_list.append(gene_id)

    # creating the data frame
    df = pd.DataFrame(zip(id_list, start_list, end_list, sequence_list, len_list),
                      columns=['id', 'start', 'end', 'sequence', 'length'])
    return df


def plot_trans_len_stat(trans_len):
    trans_len_arr = np.asarray(trans_len)

    print("\nTransmembrane Lengths stats:")
    len_max = show_stat(trans_len_arr)

    plot_hist("Transmembrane Lengths",
              trans_len_arr, "length", len_max, 8000)
    plt.tight_layout()
    plt.show()


def calculate_hydro_percent(seq_arr):
    # counting how many Hydrophobic amino are in every sequence
    seq_percent = []
    for trans_seq in seq_arr:
        count = sum([1 for a in trans_seq if a in HYDROPHOBIC_AMINO])
        hydro_percent = (count / len(trans_seq)) * 100
        seq_percent.append(hydro_percent)

    return seq_percent


def plot_trans_amino_stat(trans_seq_arr):
    print("\nTransmembrane Hydrophobic Amino stats:")

    trans_seq_percent = calculate_hydro_percent(trans_seq_arr)
    show_stat(trans_seq_percent)

    plot_hist("Transmembrane Hydrophobic Amino Percent",
              trans_seq_percent, "percent", 100, 4000)
    plt.tight_layout()
    plt.show()

# ------------B.3--------------


def group_cds_by_trans(cds, trans):
    cds_id_arr = np.asarray(cds['id'])
    trans_id_arr = np.asarray(trans['id'])

    clean_cds_id = clean_id(cds_id_arr)
    clean_trans_id = clean_id(trans_id_arr)

    mask = np.array([(cds_id in clean_trans_id) for cds_id in clean_cds_id])

    cds_with_trans = cds[mask]
    cds_without_trans = cds[~mask]
    return cds_with_trans, cds_without_trans


def plot_cds_trans_gc_percent(cds, cds_with_trans, cds_without_trans):
    all_gc_percent = np.asarray(cds['GC percent'])
    with_gc_percent = np.asarray(cds_with_trans['GC percent'])
    without_gc_percent = np.asarray(cds_without_trans['GC percent'])

    print("\nAll proteins gc percent stats:")
    show_stat(all_gc_percent)
    print("\nProteins with transmembrane gc percent stats:")
    show_stat(with_gc_percent)
    print("\nProteins without transmembrane gc percent stats:")
    show_stat(without_gc_percent)

    x_title = "GC percent"
    x_max = 100
    y_max = 1500
    plt.figure(figsize=(15, 6))

    plt.subplot(1, 4, 1)
    plot_hist("All Proteins", all_gc_percent, x_title, x_max, y_max)

    plt.subplot(1, 4, 2)
    plot_hist("Proteins with transmembrane",
              with_gc_percent, x_title, x_max, y_max)

    plt.subplot(1, 4, 3)
    plot_hist("Other proteins",
              without_gc_percent, x_title, x_max, y_max)

    plt.subplot(1, 4, 4)
    plot_hist("B and C",
              [with_gc_percent, without_gc_percent], x_title, x_max, y_max, labels=['with trans', 'others'])

    plt.suptitle("GC Percent Histograms")
    plt.tight_layout()
    plt.show()


if __name__ == "__main__":
    gb_df = open_part_a_result_file()
    up_df = open_xlsx_file()

    print("\n------------Question 1------------")
    cds_df, __ = group_genes(gb_df)                             # group to CDS and others
    show_compare_gb_to_up(cds_df['id'], up_df['Gene ID'])       # printing the comparing results

    print("\n------------Question 2------------")
    trans_df = create_trans_table(up_df)                        # creating a table of all the trans parts
    plot_trans_len_stat(trans_df['length'])                     # printing trans lengths stats
    plot_trans_amino_stat(trans_df['sequence'])                 # printing trans hydro amino stats

    print("\n------------Question 3------------")
    b_group, c_group = group_cds_by_trans(cds_df, trans_df)
    plot_cds_trans_gc_percent(cds_df, b_group, c_group)
