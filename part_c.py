from Bio.codonalign.codonseq import CodonSeq, cal_dn_ds
from generic_functions import *
from Bio.Data import CodonTable
from part_a import GenBank
from collections import Counter
import pprint
from Bio import Align

UNIPORT_PATH_2020 = os.path.join(DATA_PATH, "covid19 - 2020.gb")
UNIPORT_PATH_2022 = os.path.join(DATA_PATH, "covid19 - 2022.gb")


# ------------C.1--------------

def replace_char(string, at, to):
    temp = list(string)
    temp[at] = to
    return "".join(temp)


def show_count_synonymous():
    syn_dict = {}
    for key, value in TRANS_TABLE_1.items():
        possible_mutations = {}

        for position in range(3):
            for dna_letter in DNA:
                if key[position] == dna_letter:  # no change
                    continue

                else:
                    mutation = replace_char(key, position, dna_letter)
                    amino = TRANS_TABLE_1[mutation]
                    if amino != '*':        # stop codon
                        possible_mutations[mutation] = TRANS_TABLE_1[mutation]

        synonymous_count = Counter(possible_mutations.values())[value]
        syn_dict[key] = round((3 * synonymous_count) / len(possible_mutations), 2)

    print("\nsynonymous count dict:\n")
    pprint.pprint(syn_dict)

# ------------C.2.a--------------


def show_common_names(names1, names2):
    id_only_1 = find_missing(names1, names2)
    id_only_2 = find_missing(names2, names1)

    common_1_names = subtract_lists(names1, id_only_1)
    common_2_names = subtract_lists(names2, id_only_2)
    total_common_names = common(common_1_names, common_2_names)

    print(f"Gene names only in 2020: {id_only_1}, Common: {common_1_names}")
    print(f"Gene names only in 2022: {id_only_2}, Common: {common_2_names}")
    print(f"Total names in 2020: {len(names1)},"
          f" Total names in 2022: {len(names2)},"
          f" Total common names: {len(total_common_names)}")

    return total_common_names

# ------------C.2.b--------------


def remove_extra_from_cds_sequence(sequence):
    seq_len = len(sequence)
    seq_extra = seq_len % 3
    new_seq = sequence[:seq_len-seq_extra]
    return new_seq


def get_seq_and_trans(df, gene_name):
    gene_rows = df[df['id'] == gene_name]
    
    gene_seq = list(gene_rows['sub sequence'])[0]
    gene_trans = list(gene_rows['translation'])[0]

    gene_seq = remove_extra_from_cds_sequence(gene_seq)

    return gene_seq, gene_trans


def align(trans1, trans2):
    aligner = Align.PairwiseAligner()
    alignments = aligner.align(trans1, trans2)

    alignments = sorted(alignments)
    alignment_str = str(alignments[len(alignments) - 1])

    split = alignment_str.split('\n')
    align1 = split[0]
    align2 = split[2]

    return align1, align2


def pro_align_to_rna(sequence, alignment):
    rna_align = ""
    rna_index = 0
    for c in alignment:
        if c == '-':
            rna_align += '---'
            
        else:
            if (rna_index + 3) <= len(sequence):
                rna_align += str(sequence[rna_index:rna_index + 3])
                rna_index += 3

    rna_align = remove_extra_from_cds_sequence(rna_align)
    return rna_align


def selection_type(dn, ds):
    if dn == ds:
        return "Natural"
    elif dn > ds:
        return "Positive"
    else:
        return "Negative"


def show_dn_ds_stats(rna_align1, rna_align2):
    codon_seq1 = CodonSeq(rna_align1)
    codon_seq2 = CodonSeq(rna_align2)
    try:
        dn, ds = cal_dn_ds(codon_seq1, codon_seq2, codon_table=CodonTable.generic_by_id[1])
        print(f"dN:{dn}")
        print(f"dS:{ds}")
        print(f"{selection_type(dn, ds)} selection.")
    except KeyError:
        print('key error')


def calculate_dn_ds(gene_seq1, trans1, gene_seq2, trans2):
    align1, align2 = align(trans1, trans2)
    assert (len(align1) == len(align2))

    rna_align1 = pro_align_to_rna(gene_seq1, align1)
    rna_align2 = pro_align_to_rna(gene_seq2, align2)
    assert (len(rna_align1) == len(rna_align2))

    show_dn_ds_stats(rna_align1, rna_align2)


def genes_dn_ds(genes_names, df1, df2):
    for gene in genes_names:
        print(f'********** gene {gene} **********')

        gene_seq1, trans1 = get_seq_and_trans(df1, gene)
        gene_seq2, trans2 = get_seq_and_trans(df2, gene)

        calculate_dn_ds(gene_seq1, trans1, gene_seq2, trans2)


if __name__ == "__main__":
    genbank20 = GenBank(UNIPORT_PATH_2020, 'gene', False)
    genbank22 = GenBank(UNIPORT_PATH_2022, 'gene', False)
    genbank20.add_sub_sequence_col()
    genbank22.add_sub_sequence_col()

    print("\n------------Question 1------------")
    show_count_synonymous()

    print("\n------------Question 2a------------")
    common_names = show_common_names(genbank20.genes, genbank22.genes)

    print("\n------------Question 2b------------")
    chosen_genes = ['S', 'N', 'M', 'ORF6', 'ORF7a']
    genes_dn_ds(chosen_genes, genbank20.df, genbank22.df)
    mask = [i in chosen_genes for i in genbank20.df['id']]
    print(genbank20.df[mask])
