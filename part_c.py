import numpy as np
import matplotlib.pyplot as plt
import pandas as pd
from part_a import *
from part_b import *

UNIPORT_PATH_2020 = os.path.join(DATA_PATH, "covid19 - 2020.gb")
UNIPORT_PATH_2022 = os.path.join(DATA_PATH, "covid19 - 2022.gb")

TRANS_TABLE_1 = {
    'TTT': 'F', 'TCT': 'S', 'TAT': 'Y', 'TGT': 'C',
    'TTC': 'F', 'TCC': 'S', 'TAC': 'Y', 'TGC': 'C',
    'TTA': 'L', 'TCA': 'S', 'TAA': '*', 'TGA': '*',
    'TTG': 'L', 'TCG': 'S', 'TAG': '*', 'TGG': 'W',
    'CTT': 'L', 'CCT': 'P', 'CAT': 'H', 'CGT': 'R',
    'CTC': 'L', 'CCC': 'P', 'CAC': 'H', 'CGC': 'R',
    'CTA': 'L', 'CCA': 'P', 'CAA': 'Q', 'CGA': 'R',
    'CTG': 'L', 'CCG': 'P', 'CAG': 'Q', 'CGG': 'R',
    'ATT': 'I', 'ACT': 'T', 'AAT': 'N', 'AGT': 'S',
    'ATC': 'I', 'ACC': 'T', 'AAC': 'N', 'AGC': 'S',
    'ATA': 'I', 'ACA': 'T', 'AAA': 'K', 'AGA': 'R',
    'ATG': 'M', 'ACG': 'T', 'AAG': 'K', 'AGG': 'R',
    'GTT': 'V', 'GCT': 'A', 'GAT': 'D', 'GGT': 'G',
    'GTC': 'V', 'GCC': 'A', 'GAC': 'D', 'GGC': 'G',
    'GTA': 'V', 'GCA': 'A', 'GAA': 'E', 'GGA': 'G',
    'GTG': 'V', 'GCG': 'A', 'GAG': 'E', 'GGG': 'G'}


def count_synonymous():
    DNA = ['A', 'C', 'T', 'G']

    syn_dict = {}
    for key, value in TRANS_TABLE_1.items():
        synon = 0
        for position in range(3):
            d_count = 0
            for d in DNA:
                if key[position] == d:  # not change
                    continue
                else:
                    temp = list(key)
                    temp[position] = d
                    temp = "".join(temp)

                    if value == TRANS_TABLE_1[temp]:
                        d_count += 1
            if d_count == len(DNA) - 1:
                synon += 1
        syn_dict[key] = synon

    return syn_dict


def create_dataframe(path):
    # opening the file
    with open(path, "r") as file:
        gen = SeqIO.parse(file, "genbank")
        rec = next(gen)  # content of 1st record

    id, start, end, feat_type, strand = [], [], [], [], []
    # only for cds
    translation, codon_start = [], []

    feats = rec.features[1:]
    for feat in feats:
        if feat.type == 'gene':
            continue
        if "gene" in feat.qualifiers.keys():
            name = feat.qualifiers["gene"][0]
        else:
            name = None
        id.append(name)
        start.append(feat.location.start.position)
        end.append(feat.location.end.position)
        feat_type.append(feat.type)
        strand.append(feat.location.strand)

        if feat.type == 'CDS':
            translation.append(feat.qualifiers["translation"][0])
            codon_start.append(feat.qualifiers["codon_start"][0])
        else:
            translation.append(None)
            codon_start.append(None)

    # creating the data frame
    df = pd.DataFrame(zip(id, start, end, strand, feat_type, translation, codon_start),
                      columns=['id', 'start', 'end', 'strand', 'type', 'translation', 'codon_start'])
    return rec, df


if __name__ == "__main__":
    # Q1
    print("\n------------Question 1------------")
    dict = count_synonymous()
    print("synonymous count dict:\n", dict)

    # Q2
    print("\n------------Question 2------------")
    __, df20 = create_dataframe(UNIPORT_PATH_2020)
    __, df22 = create_dataframe(UNIPORT_PATH_2022)

    df20_id = np.asarray(df20['id'])
    df22_id = np.asarray(df22['id'])

    id_only_20 = find_missing(df20_id, df22_id)
    id_only_22 = find_missing(df22_id, df20_id)
    print("id_only_20: ", id_only_20)
    print("id_only_22: ", id_only_22)
