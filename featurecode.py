import numpy as np
from Bio import SeqIO
def cpg_oe(seq):
    seq = seq.upper().replace("T", "U")
    f_C = seq.count("C") / len(seq)
    f_G = seq.count("G") / len(seq)
    f_CpG = seq.count("CG") / (len(seq)-1)
    return f_CpG / (f_C*f_G) if (f_C*f_G) != 0 else 0
def upa_oe(seq):
    seq = seq.upper().replace("T", "U")
    f_U = seq.count("U") / len(seq)
    f_A = seq.count("A") / len(seq)
    f_UpA = seq.count("UA") / (len(seq)-1)
    return f_UpA / (f_U*f_A) if (f_U*f_A) != 0 else 0
def gc_content(seq):
    seq = seq.upper()
    return (seq.count("G")+seq.count("C"))/len(seq)
def extract_features(fasta_path):
    X = []
    for rec in SeqIO.parse(fasta_path, "fasta"):
         seq = str(rec.seq)
         features = [cpg_oe(seq), upa_oe(seq), gc_content(seq)]
         X.append(features)
    return np.array(X)
X_insect = extract_features("insect_sequences.fasta")
X_human = extract_features("human_sequences.fasta")
print("Insect feature matrix:", X_insect.shape)
print("Human feature matrix:", X_human.shape)







