from Bio.codonalign.codonseq import CodonSeq, cal_dn_ds
from part_b import *
from Bio.Data import CodonTable
from Bio import pairwise2
from Bio.Align import substitution_matrices

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


def show_count_synonymous():
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

    print(f"\nsynonymous count dict:\n{syn_dict}")

# ------------C.2.a--------------


def create_dataframe(path):
    # opening the file
    with open(path, "r") as file:
        gen = SeqIO.parse(file, "genbank")
        source = next(gen)  # content of 1st record

    id_name, start, end, feat_type, strand = [], [], [], [], []
    # only for cds
    translations, codon_starts = [], []

    feats = source.features[1:]
    for feat in feats:
        if feat.type == 'gene':
            continue

        if "gene" in feat.qualifiers.keys():
            name = feat.qualifiers["gene"][0]
        else:
            name = None

        if feat.type == 'CDS':
            trans = feat.qualifiers["translation"][0]
            codon_start = feat.qualifiers["codon_start"][0]
        else:
            trans = None
            codon_start = None

        id_name.append(name)
        start.append(feat.location.start.position)
        end.append(feat.location.end.position)
        feat_type.append(feat.type)
        strand.append(feat.location.strand)
        translations.append(trans)
        codon_starts.append(codon_start)

    # creating the data frame
    df = pd.DataFrame(zip(id_name, start, end, strand, feat_type, translations, codon_starts),
                      columns=['id', 'start', 'end', 'strand', 'type', 'translation', 'codon_start'])

    # getting the full genome
    sequence = source.seq.upper()

    return sequence, df


def show_df_common_names(df1, df2):
    up_clean_df20 = df1.dropna(subset=['translation'])
    up_clean_df22 = df2.dropna(subset=['translation'])

    df20_id = np.asarray(up_clean_df20['id'])
    df22_id = np.asarray(up_clean_df22['id'])

    id_only_20 = find_missing(df20_id, df22_id)
    id_only_22 = find_missing(df22_id, df20_id)

    print("id_only_20: ", id_only_20)
    print("id_only_22: ", id_only_22)

# ------------C.2.b--------------


def get_seq_and_trans(df, sequence, gene_name):
    gene_row = df[df['id'] == gene_name]
    
    start = list(gene_row['start'])[0]
    end = list(gene_row['end'])[0]
    
    gene_seq = sub_sequence(sequence, start, end)
    gene_trans = list(gene_row['translation'])[0]
    
    return gene_seq, gene_trans


def align(trans1, trans2):
    matrix = substitution_matrices.load("BLOSUM62")
    alignment = pairwise2.align.globaldx(trans1, trans2, matrix)
    align1 = alignment[0][0]
    align2 = alignment[0][1]
    print(align1, align2)

    return align1, align2      


def pro_align_to_rna(sequence, alignment):
    rna_align = ""
    rna_index = 0
    for c in alignment:
        if c == '-':
            rna_align += '---'
            
        else:
            rna_align += str(sequence[rna_index:rna_index + 3])
            rna_index += 3

    return rna_align


def dn_ds_stats(rna_align1, rna_align2):
    codon_seq1 = CodonSeq(rna_align1)
    codon_seq2 = CodonSeq(rna_align2)
    try:
        dn, ds = cal_dn_ds(codon_seq1, codon_seq2, codon_table=CodonTable.generic_by_id[1])
        print(f"dN:{dn}")
        print(f"dS:{ds}")
    except KeyError:
        print('key error')


def show_genes_dnds_stats(seq1, df1, seq2, df2):
    clean_df1 = df1.dropna(subset=['id'])
    genes_names = np.asarray(clean_df1['id'])

    for gene in genes_names:
        print(f'\n**********gene {gene} **********\n')
        
        gene_seq1, trans1 = get_seq_and_trans(df1, seq1, gene)
        gene_seq2, trans2 = get_seq_and_trans(df2, seq2, gene)
        
        align1, align2 = align(trans1, trans2)

        rna_align1 = pro_align_to_rna(gene_seq1, align1)
        rna_align2 = pro_align_to_rna(gene_seq2, align2)

        # dn_ds_stats(rna_align1, rna_align2)
    

if __name__ == "__main__":
    print("\n------------Question 1------------")
    show_count_synonymous()

    print("\n------------Question 2a------------")
    genome1, df20 = create_dataframe(UNIPORT_PATH_2020)
    genome2, df22 = create_dataframe(UNIPORT_PATH_2022)
    show_df_common_names(df20, df22)

    print("\n------------Question 2b------------")    
    show_genes_dnds_stats(genome1, df20, genome2, df22)
 