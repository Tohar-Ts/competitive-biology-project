import numpy as np
import matplotlib.pyplot as plt
import pandas as pd
from Bio.Seq import Seq
from Bio import Align
from Bio.Align import substitution_matrices
from Bio.codonalign.codonseq import CodonSeq, cal_dn_ds
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


def split_align(alignments):
    position = int(len(alignments)/3)
    align1 = alignments[:position]
    align2 = alignments[position*2:]
    # align3 = alignments[position:position*2]
    # print(alignments)
    # print("align1: ", align1)
    # # print("seq3: ", seq3)
    # print("align2: ", align2)
    return align1, align2


def pro_align_to_rna(seq, align):
    # print("seq : ", seq)
    # print("align : ", align)
    rna_align = []
    for ind, c in enumerate(align):
        if c == '-':
            rna_align = np.concatenate((rna_align, ['-']*3))
            continue
        else:
            # print("seq[ind*3:ind*3 + 3]", seq[ind*3:ind*3 + 3])
            if ind*3 + 3 <= len(seq):
                rna_align = np.concatenate((rna_align, seq[ind*3:ind*3 + 3]))
                print("seq[ind*3:ind*3 + 3]", seq[ind*3:ind*3 + 3])
    print("".join(rna_align))
    return "".join(rna_align)


def seq_codon(seq1, translation1, seq2, translation2):

    aligner = Align.PairwiseAligner()
    aligner.substitution_matrix = substitution_matrices.load("BLOSUM62")
    # print("translation1", translation1)
    # print("translation2", translation2)
    alignments = aligner.align(translation1, translation2)
    alignments = list(alignments)
    print("alignments:\n", alignments[0])
    align1, align2 = split_align(str(alignments[0]))
    rna_align1 = pro_align_to_rna(seq1, align1)
    rna_align2 = pro_align_to_rna(seq2, align2)
    print("lens:", len(rna_align1), len(rna_align2))
    seq1 = CodonSeq(rna_align1)
    seq2 = CodonSeq(rna_align2)

    dN, dS = cal_dn_ds(seq1, seq2)
    dN_dS_ratio = float(dN/dS)
    print("dN:%0.3f " % dN)
    print("dS:%0.3f " % dS)
    print("dN/dS:%0.3f " % dN_dS_ratio)
    # return dN, dS, dN_dS_ratio
    return 0, 0, 0


def Comparison_genes(seq1, df20, seq2, df22, names_genes):

    for gene in names_genes:
        print('gene', gene)
        # print(df20['id'])
        loc1 = df20[df20['id'] == gene]
        loc2 = df22[df22['id'] == gene]
        # print("loc1:", loc1)
        # print(list(loc1['translation'])[0], list(loc2['translation'])[0])
        start1 = list(loc1['start'])[0]
        end1 = list(loc1['end'])[0]
        start2 = list(loc2['start'])[0]
        end2 = list(loc2['end'])[0]
        dN, dS, dN_dS_ratio = seq_codon(seq1[start1: end1 + 1],
                                        list(loc1['translation'])[0], seq2[start2: end2 + 1], list(loc2['translation'])[0])


if __name__ == "__main__":
    # Q1
    print("\n------------Question 1------------")
    dict = count_synonymous()
    print("synonymous count dict:\n", dict)

    # Q2
    print("\n------------Question 2------------")
    rec1, df20 = create_dataframe(UNIPORT_PATH_2020)
    rec2, df22 = create_dataframe(UNIPORT_PATH_2022)
    # print(df20)
    up_clean_df20 = df20.dropna(subset=['translation'])
    up_clean_df22 = df22.dropna(subset=['translation'])
    df20_id = np.asarray(up_clean_df20['id'])
    df22_id = np.asarray(up_clean_df22['id'])

    id_only_20 = find_missing(df20_id, df22_id)
    id_only_22 = find_missing(df22_id, df20_id)
    print("id_only_20: ", id_only_20)
    print("id_only_22: ", id_only_22)

    # print(np.asarray(df20.keys()))
    # print(np.asarray(df22.keys()))
    print(up_clean_df20)
    five_common_genes = ['ORF3a', 'ORF6', 'ORF7a', 'ORF8', 'ORF10']
    Comparison_genes(rec1.seq.upper(), up_clean_df20,
                     rec2.seq.upper(), up_clean_df22, five_common_genes)
