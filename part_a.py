import os
from Bio import SeqIO
import numpy as np
import matplotlib.pyplot as plt
import pandas as pd
from Bio.Seq import Seq
from Bio.Data.CodonTable import TranslationError

DATA_PATH = os.path.join(os.path.dirname(os.path.abspath(__file__)), "data")
UNIPORT_PATH = os.path.join(DATA_PATH, 'BS168.gb')
RESULT_PATH = os.path.join(DATA_PATH, "part_a.csv")
EXCEPTION_PATH = os.path.join(DATA_PATH, "gene_exceptions.csv")

# ------------Global Functions--------------


def open_and_create_gb_dataframe(path):
    # opening the file
    with open(path, "r") as file:
        gen = SeqIO.parse(file, "genbank")
        source = next(gen)  # content of 1st record

    gene_names = []
    locus_id, start, end, feat_types, strand = [], [], [], [], []
    # only for cds
    tables, translations, codon_starts = [], [], []

    feats = source.features[1:]
    for feat in feats:
        locus_name = feat.qualifiers["locus_tag"][0]
        feat_type = feat.type

        if feat_type == 'gene':
            gene_names.append(locus_name)
            continue

        if feat_type == 'CDS':
            table = feat.qualifiers["transl_table"][0]
            translation = feat.qualifiers["translation"][0]
            codon_start = feat.qualifiers["codon_start"][0]
        else:
            table = None
            translation = None
            codon_start = None

        locus_id.append(locus_name)
        feat_types.append(feat_type)
        start.append(feat.location.start.position)
        end.append(feat.location.end.position)
        strand.append(feat.location.strand)
        tables.append(table)
        translations.append(translation)
        codon_starts.append(codon_start)

    # creating the data frame
    df = pd.DataFrame(zip(locus_id, start, end, strand, feat_types, tables, translations, codon_starts),
                      columns=['id', 'start', 'end', 'strand', 'type', 'table', 'translation', 'codon_start'])

    # getting the full genome
    sequence = source.seq.upper()

    return gene_names, sequence, df


def single_gene_length(start, end):
    return np.abs(end - start)


def show_stat(arr):
    if type(arr) != np.ndarray:
        arr = np.asarray(arr)

    minimum = np.min(arr)
    maximum = np.max(arr)
    avg = np.average(arr)
    med = np.median(arr)

    print(f"average: {avg:.2f} minimum: {minimum:.2f} maximum: {maximum:.2f} median: {med:.2f}")
    return maximum


def plot_hist(title, from_arr, x_label, x_max, y_max, labels=None):
    plt.style.use('seaborn-deep')
    plt.title(title)
    if labels:
        plt.hist(from_arr, label=labels)
        plt.legend(loc='upper right')
    else:
        plt.hist(from_arr)
    plt.xlabel(x_label)
    plt.ylabel("number count")
    plt.xlim([0, x_max])
    plt.ylim([0, y_max])

# ------------GenBank Object--------------


class GenBank:
    def __init__(self, file_path):
        genes, seq, df = open_and_create_gb_dataframe(file_path)
        self.df = gb_df
        self.seq = seq
        self.genes = genes

    # ------------A.1--------------

    def show_type_dict(self):
        gene_dict = {'gene': self.genes}
        types_dict = self.df.type.value_counts().to_dict()
        total_dict = dict(gene_dict)
        total_dict.update(types_dict)
        print(f"\nTypes dictionary:\n{total_dict}")

    # ------------A.2.a------------

    def add_len_col(self):
        new_col = self.df.apply(lambda row: single_gene_length(row['start'], row['end']), axis=1)
        self.df['len'] = new_col

    # ------------A.2.b------------

    def group_genes(self):
        cds = self.df[self.df['type'] == 'CDS']
        others = self.df[self.df['type'] != 'CDS']
        return cds, others


# ------------A.2.d------------




def plot_len_stat(cds_len, others_len):
    print("\nProteins lengths stats:")
    cds_len_max = show_stat(cds_len)

    print("\nNon proteins lengths stats:")
    other_len_max = show_stat(others_len)

    cds_len_arr = np.asarray(cds_len)
    others_len_arr = np.asarray(others_len)

    x_max = max([cds_len_max, other_len_max])
    y_max = 4500

    all_genes = np.concatenate((cds_len_arr, others_len_arr))

    x_label = "length"
    plt.figure(figsize=(10, 6))

    plt.subplot(1, 3, 1)
    plot_hist("All Genes", all_genes, x_label, x_max, y_max)

    plt.subplot(1, 3, 2)
    plot_hist("Proteins", cds_len_arr, x_label, x_max, y_max)

    plt.subplot(1, 3, 3)
    plot_hist("Non Proteins", others_len_arr, x_label, x_max, y_max)

    plt.suptitle("Lengths Histograms")
    plt.tight_layout()
    plt.show()

# ------------A.3.a------------


def gc_percentage(sequence):
    g_count = sequence.count('G')
    c_count = sequence.count('C')
    return ((g_count + c_count) / len(sequence)) * 100


def show_source_gc(sequence):
    gc_percent = gc_percentage(sequence)
    print(f"\nGenome GC percent: {gc_percent:.2f}%")

# ------------A.3.b------------


def sub_sequence(sequence, start, end):
    return sequence[start:end]


def add_sub_sequence_col(df, sequence):
    new_col = df.apply(lambda row: sub_sequence(sequence, row['start'], row['end']), axis=1)
    df['sub sequence'] = new_col
    return df


def gc_percent_for_row(row):
    if row['type'] != 'CDS':
        return None
    return gc_percentage(row['sub sequence'])


def add_gc_percent_col(df):
    new_col = df.apply(lambda row: gc_percent_for_row(row), axis=1)
    df['GC percent'] = new_col
    return df

# ------------A.3.c------------


def mean(lst):
    return lst.mean()


def show_avg_gc_percent_col(df):
    average = mean(df["GC percent"])
    print(f"\nProteins average GC percent: {average:.2f}%")

# ------------A.3.d------------


def plot_gc_stat(gc_percent_arr):
    plot_hist("GC Percent Histogram", gc_percent_arr, "GC percent", 100, 1600)
    plt.show()

# ------------A.3.d------------


def show_extreme_gc_percents_genes(df, n=5):
    top = df.nlargest(n, 'GC percent')
    bottom = df.nsmallest(n, 'GC percent')
    print("\nExtreme GC percents genes:")
    print(f"\nTop {n} GC percents genes details: ")
    print(top)
    print(f"\nBottom {n} GC percents genes details: ")
    print(bottom)

# ------------A.4------------


def check_translation(row):
    row_type = row["type"]

    if row_type != "CDS":
        return None

    org_trans = row["translation"]
    sequence = str(row["sub sequence"])
    table = row["table"]
    strand = row["strand"]
    codon_start = int(row["codon_start"])

    org_seq = Seq(sequence[(codon_start-1):])
    if strand != 1:
        org_seq = org_seq.reverse_complement()

    try:
        trans = str(org_seq.translate(table=table, cds=True))
        check = (trans == org_trans)
        if check:
            return "OK"
        return "Translation doesn't match"

    except TranslationError as err:
        return str(err)


def find_conflicts(df):
    new_col = df.apply(lambda row: check_translation(row), axis=1)
    df['check'] = new_col

    check_df = df.dropna(subset='check')
    conf_df = check_df[check_df['check'] != 'OK']

    print(f"\nGenes with conflict in translations:\n{conf_df}")
    return conf_df


if __name__ == "__main__":
    total_genes, seq, gb_df = open_and_create_gb_dataframe(UNIPORT_PATH)

    print("\n------------Question 1------------")
    show_type_dict(gb_df, total_genes)

    print("\n------------Question 2------------")
    gb_df = add_len_col(gb_df)                          # calculate the len
    cds_df, other_gene = group_genes(gb_df)             # group to CDS and others
    plot_len_stat(cds_df['len'], other_gene['len'])     # print the stat and plot the histogram

    print("\n------------Question 3------------")
    show_source_gc(seq)                                 # print GC percent of the genome
    gb_df = add_sub_sequence_col(gb_df, seq)
    gb_df = add_gc_percent_col(gb_df)                   # create GC percent column
    show_avg_gc_percent_col(gb_df)                      # print average GC percent of proteins
    plot_gc_stat(gb_df['GC percent'])                   # plot GC histogram
    show_extreme_gc_percents_genes(gb_df)
    gb_df.to_csv(RESULT_PATH)                           # save final results to csv file

    print("\n------------Question 4------------")
    conflicts_df = find_conflicts(gb_df)
    conflicts_df.to_csv(EXCEPTION_PATH)                 # save exceptions results to csv file
