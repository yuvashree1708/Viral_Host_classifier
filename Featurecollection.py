import csv
from Bio import SeqIO
def cpg_oe(seq):
    seq=seq.upper().replace("T","U")
    f_C=seq.count("C")/len(seq)
    f_G=seq.count("G")/len(seq)
    f_CpG=seq.count("CG")/len(seq)-1
    if f_C*
