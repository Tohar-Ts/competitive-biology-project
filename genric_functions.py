import pandas as pd
import numpy as np
from Bio import SeqIO
import matplotlib.pyplot as plt
import os
from Bio.Data.CodonTable import TranslationError
from Bio.Seq import Seq

pd.set_option('display.expand_frame_repr', False)
DATA_PATH = os.path.join(os.path.dirname(os.path.abspath(__file__)), "data")
HYDROPHOBIC_AMINO = ['A', 'F', 'L', 'I', 'V', 'M', 'P', 'W']


def open_xlsx_file(path, sheet_name):
    return pd.read_excel(path, sheet_name=sheet_name)


def open_csv_file(path):
    return pd.read_csv(path)


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


def sub_sequence(sequence, start, end):
    return sequence[start:end]


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


def group_genes(df):
    cds = df[df['type'] == 'CDS']
    not_cds = df[df['type'] != 'CDS']
    return cds, not_cds


def gc_percentage(sequence):
    if sequence is None:
        return None
    g_count = sequence.count('G')
    c_count = sequence.count('C')
    return ((g_count + c_count) / len(sequence)) * 100


def check_translation(gene_type, org_trans, sequence, trans_table, strand, codon_start):
    if gene_type != "CDS":
        return None

    codon_start = int(codon_start)
    org_seq = Seq(sequence[(codon_start - 1):])
    if strand != 1:
        org_seq = org_seq.reverse_complement()

    try:
        trans = str(org_seq.translate(table=trans_table, cds=True))
        check = (trans == org_trans)
        if check:
            return "OK"
        return "Translation doesn't match"

    except TranslationError as err:
        return str(err)


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

        multi_id = str(name).split()
        for i in multi_id:
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


def calculate_hydro_percent(seq_arr):
    # counting how many Hydrophobic amino are in every sequence
    seq_percent = []
    for seq in seq_arr:
        count = sum([1 for a in seq if a in HYDROPHOBIC_AMINO])
        hydro_percent = (count / len(seq)) * 100
        seq_percent.append(hydro_percent)

    return seq_percent
