from Bio import SeqIO 
import numpy as np
def cpg_oe(seq):
    seq = seq.upper().replace("T", "U")  # RNA-style
    f_C = seq.count("C") / len(seq)
    f_G = seq.count("G") / len(seq)
    f_CpG = seq.count("CG") / (len(seq)-1)  # divide by possible dinucleotide positions
    return f_CpG / (f_C * f_G) if (f_C * f_G) != 0 else 0
def upa_oe(seq):
    seq = seq.upper().replace("T", "U")
    f_U = seq.count("U") / len(seq)
    f_A = seq.count("A") / len(seq)
    f_UpA = seq.count("UA") / (len(seq)-1)
    return f_UpA / (f_U * f_A) if (f_U * f_A) != 0 else 0
def gc_content(seq):
    seq = seq.upper()
    return (seq.count("G") + seq.count("C")) / len(seq)
def load_and_extract_features(fasta_path, label):
    X=[]
    Y=[]
    for rec in SeqIO.parse(fasta_path,"fasta"):
        seq=str(rec.seq)
        features = [cpg_oe(seq),upa_oe(seq),gc_content(seq)]
        X.append(features)
        Y.append(label)
    return X, Y
X_insect, y_insect = load_and_extract_features("human_sequences.fasta", 0)
X_human, y_human = load_and_extract_features("insect_sequences.fasta", 1)
X = np.array(X_insect + X_human)
Y = np.array(y_insect + y_human)
print("Feature matrix shape:", X.shape) 
print("Labels shape:", Y.shape)






