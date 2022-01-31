from generic_functions import *


UNIPORT_PATH = os.path.join(DATA_PATH, 'BS168.gb')
PART_A_RESULT_PATH = os.path.join(DATA_PATH, "part_a.csv")
EXCEPTION_PATH = os.path.join(DATA_PATH, "gene_exceptions.csv")


class GenBank:
    def __init__(self, file_path, id_header="locus_tag",  has_trans_table=True):
        genes, seq, df = open_and_create_gb_dataframe(file_path, id_header,  has_trans_table)
        self.df = df
        self.seq = seq
        self.genes = genes
        self.cds = None
        self.not_cds = None

    # ------------A.1--------------

    def show_type_dict(self):
        gene_dict = {'gene': len(self.genes)}
        types_dict = self.df.type.value_counts().to_dict()
        total_dict = dict(gene_dict)
        total_dict.update(types_dict)
        print(f"\nTypes dictionary:\n{total_dict}")

    # ------------A.2.a------------

    def add_len_col(self):
        new_col = self.df.apply(lambda row: single_gene_length(row['start'], row['end']), axis=1)
        self.df['len'] = new_col

    # ------------A.2.b------------

    def set_cds_and_not_cds(self):
        cds, not_cds = group_genes(self.df)
        self.cds = cds
        self.not_cds = not_cds

    # ------------A.2.d------------

    def plot_len_stat(self):
        all_genes_len = np.asarray(self.df['len'])
        cds_len = np.asarray(self.cds['len'])
        others_len = np.asarray(self.not_cds['len'])

        print("\nProteins lengths stats:")
        cds_len_max = show_stat(cds_len)

        print("\nNon proteins lengths stats:")
        other_len_max = show_stat(others_len)

        x_max = max([cds_len_max, other_len_max])
        y_max = 4500

        x_label = "length"
        plt.figure(figsize=(10, 6))

        plt.subplot(1, 3, 1)
        plot_hist("All Genes", all_genes_len, x_label, x_max, y_max)

        plt.subplot(1, 3, 2)
        plot_hist("Proteins", cds_len, x_label, x_max, y_max)

        plt.subplot(1, 3, 3)
        plot_hist("Non Proteins", others_len, x_label, x_max, y_max)

        plt.suptitle("Lengths Histograms")
        plt.tight_layout()
        plt.show()

    # ------------A.3.a------------

    def show_source_gc(self):
        gc_percent = gc_percentage(self.seq)
        print(f"\nGenome GC percent: {gc_percent:.2f}%")

    # ------------A.3.b------------

    def add_sub_sequence_col(self):
        new_col = self.df.apply(lambda row: sub_sequence(self.seq, row['start'], row['end']), axis=1)
        self.df['sub sequence'] = new_col

    def add_gc_percent_col(self):
        new_col = self.df.apply(lambda row: gc_percentage(row['sub sequence']), axis=1)
        self.df['GC percent'] = new_col
        self.set_cds_and_not_cds()

    # ------------A.3.c------------

    def plot_cds_gc_percent_stats(self):
        cds_gc_percent = self.cds["GC percent"]

        average = cds_gc_percent.mean()
        print(f"\nProteins average GC percent: {average:.2f}%")

        plot_hist("Proteins GC Percent Histogram", cds_gc_percent, "GC percent", 100, 1600)
        plt.show()

    # ------------A.3.d------------

    def show_extreme_gc_percents_genes(self, n=5):
        top = self.df.nlargest(n, 'GC percent')
        bottom = self.df.nsmallest(n, 'GC percent')

        print("\nExtreme GC percents genes:")
        print(f"\nTop {n} GC percents genes details: ")
        print(top)
        print(f"\nBottom {n} GC percents genes details: ")
        print(bottom)

    # ------------A.4------------

    def report_conflicts(self, results_path):
        self.df['check'] = self.df.apply(lambda row: check_translation(row["type"], row["translation"],
                                                                       row["sub sequence"], row["table"],
                                                                       row["strand"], row["codon_start"]), axis=1)
        check_df = self.df.dropna(subset=['check'])
        conf_df = check_df[check_df['check'] != 'OK']

        print(f"\nGenes with conflict in translations:\n{conf_df}")
        conf_df.to_csv(results_path)  # save exceptions results to csv file


if __name__ == "__main__":
    genbank = GenBank(UNIPORT_PATH)

    print("\n------------Question 1------------")
    genbank.show_type_dict()

    print("\n------------Question 2------------")
    genbank.add_len_col()        # calculate the len
    genbank.set_cds_and_not_cds()       # group to CDS and others
    genbank.plot_len_stat()     # print the stat and plot the histogram

    print("\n------------Question 3------------")
    genbank.show_source_gc()        # print GC percent of the genome
    genbank.add_sub_sequence_col()
    genbank.add_gc_percent_col()        # create GC percent column
    genbank.plot_cds_gc_percent_stats()       # print average GC percent of proteins and plot GC histogram
    genbank.show_extreme_gc_percents_genes()

    print("\n------------Question 4------------")
    genbank.report_conflicts(EXCEPTION_PATH)
