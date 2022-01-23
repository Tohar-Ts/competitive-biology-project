import os
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
from Bio.SeqFeature import SeqFeature, FeatureLocation
from Bio import SeqIO
from pytest import skip
import numpy as np
from statistics import mean
import matplotlib.pyplot as plt

file_path = "../BS168.gb"


def open_file(file_path):
    # read file:
    with open(file_path, "r") as file:
        gen = SeqIO.parse(file, "genbank")
        record_gb = next(gen)  # content of 1st record
    return record_gb


def create_type_dict(rec):
    type_dic = {}
    feats = rec.features[1:]  # skip "Source" feture
    for feat in feats:
        feat_type = feat.type
        if feat_type not in type_dic:
            type_dic[feat_type] = 1
        else:
            type_dic[feat_type] += 1
    return type_dic


def group_genes(rec):
    feats = rec.features[1:]  # skip "Source" feture

    cds, other_gene = [], []
    for feat in feats:
        if feat.type == 'gene':
            continue
        if feat.type == 'CDS':
            cds.append(feat)
        else:
            other_gene.append(feat)

    return cds, other_gene


def singel_gene_length(gene):
    start_s = gene.location.start.position
    end_s = gene.location.end.position
    length = np.abs(end_s - start_s) + 1
    return length


def Q1(rec):
    print("Question 1:")
    t_dict = create_type_dict(rec)
    print("Result:\n", t_dict, "\n")


def stat(lst):
    arr = np.asarray(lst)
    min = np.min(arr)
    max = np.max(arr)
    avg = np.average(arr)

    return avg, min, max


def create_histogram(lst1, name1, lst2, name2, lst3, name3):
    plt.title("histogram graphs")

    plt.subplot(1, 3, 1)
    plt.title(name1)
    plt.hist(lst1)

    plt.subplot(1, 3, 2)
    plt.title(name2)
    plt.hist(lst2)

    plt.subplot(1, 3, 3)
    plt.title(name3)
    plt.hist(lst3)

    plt.show()


def Q2(rec):
    print("Question 2:")
    cds_arr, other_gene_arr = group_genes(rec)

    cds_length = [singel_gene_length(feat) for feat in cds_arr]
    other_gene_length = [singel_gene_length(feat) for feat in other_gene_arr]

    avg_cds, min_cds, max_cds = stat(cds_length)
    avg_others, min_others, max_others = stat(other_gene_length)

    create_histogram(cds_length, "cds", other_gene_length, "Other gene length", cds_length + other_gene_length,
                     "All the genes")

    print("Result A:\n")

    return cds_arr, other_gene_arr


def GC_percentage(gene):
    G_count = gene.count('G')
    C_count = gene.count('C')
    GC_count = G_count + C_count

    return GC_count / len(gene)


def CDS_GC_avg_precent(genome, cds_arr):
    CDS_GC_percent = []
    for cds in cds_arr:
        start_s = cds.location.start.position
        end_s = cds.location.end.position
        gene = genome[start_s:end_s + 1]
        CDS_GC_percent.append(GC_percentage(gene))
    arr = np.asarray(CDS_GC_percent)
    CDS_GC_avg = np.average(arr)
    print("CDS_GC_avg", CDS_GC_avg)

    return arr, CDS_GC_avg


def Q3(rec, cds_arr):
    print("Question 3:")
    genome = rec.seq.upper()
    GC_percent = GC_percentage(genome)
    print("3a:\nGC percent: ", "{:.2f}".format(GC_percent), "%")

    CDS_GC_percent_arr, CDS_GC_avg = CDS_GC_avg_precent(genome, cds_arr)

    plt.title("histogram CDS GC percent arr")
    plt.hist(CDS_GC_percent_arr)
    plt.show()

def Q4(rec):
    print("Question 4:")


if __name__ == "__main__":
    rec = open_file(file_path)
    # Q1(rec)
    cds_arr, other_gene_arr = Q2(rec)
    Q3(rec, cds_arr)
    # Q4(rec)
