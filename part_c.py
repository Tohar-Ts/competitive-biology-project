import numpy as np
import matplotlib.pyplot as plt
import pandas as pd
from Bio.Seq import Seq
from Bio import Align
from Bio.Align import substitution_matrices
from Bio.codonalign.codonseq import CodonSeq, cal_dn_ds
from part_a import *
from part_b import *
from Bio.Data import CodonTable
from Bio.Align import MultipleSeqAlignment
from Bio.SeqRecord import SeqRecord
from Bio import pairwise2


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

# ------------C.1--------------
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

# ------------C.2--------------
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


def get_seq_and_trans(df, seq, gene_name):
    gene_row = df[df['id'] == gene_name]
    
    start = list(gene_row['start'])[0]
    end = list(gene_row['end'])[0]
    
    gene_seq = seq[start: end + 1]
    gene_trans = list(gene_row['translation'])[0]
    
    return gene_seq, gene_trans


def align(trans1, trans2):
    alignment = pairwise2.align.globalxx(trans1, trans2)
    align1 = alignment[0][0]
    align2 = alignment[0][1]
    return align1, align2      


def pro_align_to_rna(seq, align):
    rna_align = ""
    rna_index = 0
    for c in align:
        if c == '-':
            continue
            
        else:
            rna_align += str(seq[rna_index:rna_index + 3])
            rna_index += 3

    return rna_align


def dn_ds_stats(rna_align1, rna_align2):
    codon_seq1 = CodonSeq(rna_align1)
    codon_seq2 = CodonSeq(rna_align2)
    
    dN, dS = cal_dn_ds(codon_seq1, codon_seq2,
                       codon_table=CodonTable.generic_by_id[1])
    
    print("dN:%0.3f " % dN)
    print("dS:%0.3f " % dS)


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


def genes_stats(seq1, df20, seq2, df22, names_genes):

    for gene in names_genes:
        print(f'**********gene {gene} **********')
        
        gene_seq1, trans1 = get_seq_and_trans(df20, seq1, gene)
        gene_seq2, trans2 = get_seq_and_trans(df22, seq2, gene)
        
        align1, align2 = align(trans1, trans2)

        rna_align1 = pro_align_to_rna(gene_seq1, align1)
        rna_align2 = pro_align_to_rna(gene_seq2, align2)

        dn_ds_stats(rna_align1, rna_align2)
    

if __name__ == "__main__":
    # Q1
    print("\n------------Question 1------------")
    dict = count_synonymous()
    print("synonymous count dict:\n", dict)

    # # Q2
    print("\n------------Question 2a------------")
    rec1, df20 = create_dataframe(UNIPORT_PATH_2020)
    rec2, df22 = create_dataframe(UNIPORT_PATH_2022)

    up_clean_df20 = df20.dropna(subset=['translation'])
    up_clean_df22 = df22.dropna(subset=['translation'])
    df20_id = np.asarray(up_clean_df20['id'])
    df22_id = np.asarray(up_clean_df22['id'])

    id_only_20 = find_missing(df20_id, df22_id)
    id_only_22 = find_missing(df22_id, df20_id)
    print("id_only_20: ", id_only_20)
    print("id_only_22: ", id_only_22)
    
    print("\n------------Question 2b------------")
    five_common_genes = ['ORF3a', 'M', 'ORF7a', 'ORF3a', 'ORF10']
    
    genes_stats(rec1.seq.upper(), up_clean_df20,
                     rec2.seq.upper(), up_clean_df22, five_common_genes)
