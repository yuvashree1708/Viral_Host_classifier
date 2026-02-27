import csv
from Bio import SeqIO

def cpg_oe(seq):
    seq = seq.upper().replace("T", "U")
    L = len(seq)
    if L < 2:
        return 0.0
    f_C = seq.count("C") / L
    f_G = seq.count("G") / L
    f_CpG = seq.count("CG") / (L - 1)
    if f_C * f_G == 0:
        return 0.0
    return f_CpG / (f_C * f_G)

def upa_oe(seq):
    # match CpG function: T -> U, then compute UpA O/E
    seq = seq.upper().replace("T", "U")
    L = len(seq)
    if L < 2:
        return 0.0
    f_U = seq.count("U") / L
    f_A = seq.count("A") / L
    f_UpA = seq.count("UA") / (L - 1)
    if f_U * f_A == 0:
        return 0.0
    return f_UpA / (f_U * f_A)

def gc_content(seq):
    seq = seq.upper()
    L = len(seq)
    if L == 0:
        return 0.0
    # parentheses are important: (G + C) / L
    return (seq.count("G") + seq.count("C")) / L

def extract_features(fasta_file, csv_output):
    rows = []
    for record in SeqIO.parse(fasta_file, "fasta"):
        seq_id = record.id
        sequence = str(record.seq).strip()
        if len(sequence) == 0:
            continue
        cpg_value = cpg_oe(sequence)
        upa_value = upa_oe(sequence)
        gc_value = gc_content(sequence)
        rows.append([seq_id, cpg_value, upa_value, gc_value])

    with open(csv_output, "w", newline="") as f:
        writer = csv.writer(f)
        writer.writerow(["Sequence_ID", "CpG_OE", "UpA_OE", "GC_content"])
        writer.writerows(rows)
    print(f"Wrote {len(rows)} rows to {csv_output}")


extract_features("/content/human_sequences.fasta", "human_seq_features.csv")
extract_features("/content/insect_sequences.fasta", "insect_seq_features.csv")
