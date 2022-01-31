from part_a import PART_A_RESULT_PATH
from generic_functions import *

UNIPORT_PATH = os.path.join(DATA_PATH, "BS168.xlsx")


class GenbankWithUniport:

    def __init__(self, genbank_path=PART_A_RESULT_PATH, uniport_path=UNIPORT_PATH):
        self.gb_df = open_csv_file(genbank_path)
        self.up_df = open_xlsx_file(uniport_path, 'Sheet0')
        self.trans_df = None
        gb_cds_df, __ = group_genes(self.gb_df)
        self.gb_cds_df = gb_cds_df
        self.gc_cds_with_trans = None
        self.gc_cds_without_trans = None

    # ------------B.1--------------

    def show_compare_gb_to_up(self):
        cds_id_arr = np.asarray(self.gb_cds_df['id'])
        up_id_arr = np.asarray(self.up_df['Gene ID'])

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

    def create_trans_table(self):
        up_clean_df = self.up_df.dropna(subset=['Transmembrane'])

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
        self.trans_df = pd.DataFrame(zip(id_list, start_list, end_list, sequence_list, len_list),
                                     columns=['id', 'start', 'end', 'sequence', 'length'])

    def plot_trans_len_stat(self):
        trans_len_arr = np.asarray(self.trans_df['length'])

        print("\nTransmembrane Lengths stats:")
        len_max = show_stat(trans_len_arr)

        plot_hist("Transmembrane Lengths",
                  trans_len_arr, "length", len_max, 8000)
        plt.tight_layout()
        plt.show()

    def plot_trans_amino_stat(self):
        print("\nTransmembrane Hydrophobic Amino stats:")

        trans_seq_hydro_percent = calculate_hydro_percent(self.trans_df['sequence'])
        show_stat(trans_seq_hydro_percent)

        plot_hist("Transmembrane Hydrophobic Amino Percent",
                  trans_seq_hydro_percent, "percent", 100, 4000)
        plt.tight_layout()
        plt.show()

    # ------------B.3--------------

    def set_cds_with_and_without_trans(self):
        cds_id_arr = np.asarray(self.gb_cds_df['id'])
        trans_id_arr = np.asarray(self.trans_df['id'])

        clean_cds_id = clean_id(cds_id_arr)
        clean_trans_id = clean_id(trans_id_arr)

        mask = np.array([(cds_id in clean_trans_id) for cds_id in clean_cds_id])

        cds_with_trans = self.gb_cds_df[mask]
        cds_without_trans = self.gb_cds_df[~mask]

        self.gc_cds_with_trans = cds_with_trans
        self.gc_cds_without_trans = cds_without_trans

    def plot_cds_trans_gc_percent(self):
        all_gc_percent = np.asarray(self.gb_cds_df['GC percent'])
        with_gc_percent = np.asarray(self.gc_cds_with_trans['GC percent'])
        without_gc_percent = np.asarray(self.gc_cds_without_trans['GC percent'])

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
        plot_hist("All Genbank Proteins", all_gc_percent, x_title, x_max, y_max)

        plt.subplot(1, 4, 2)
        plot_hist("Genbank proteins with transmembrane",
                  with_gc_percent, x_title, x_max, y_max)

        plt.subplot(1, 4, 3)
        plot_hist("Other Genbank proteins",
                  without_gc_percent, x_title, x_max, y_max)

        plt.subplot(1, 4, 4)
        plot_hist("B and C",
                  [with_gc_percent, without_gc_percent], x_title, x_max, y_max, labels=['with trans', 'others'])

        plt.suptitle("GC Percent Histograms")
        plt.tight_layout()
        plt.show()


if __name__ == "__main__":
    gb_and_up = GenbankWithUniport(PART_A_RESULT_PATH, UNIPORT_PATH)

    print("\n------------Question 1------------")
    gb_and_up.show_compare_gb_to_up()       # printing the comparing results

    print("\n------------Question 2------------")
    gb_and_up.create_trans_table()      # creating a table of all the trans parts
    gb_and_up.plot_trans_len_stat()       # printing trans lengths stats
    gb_and_up.plot_trans_amino_stat()       # printing trans hydro amino stats
    print("\n", gb_and_up.trans_df)
    print("\n------------Question 3------------")
    gb_and_up.set_cds_with_and_without_trans()       # group cds by trans
    gb_and_up.plot_cds_trans_gc_percent()       # plot a, b, c, b+c subplots
