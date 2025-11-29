import csv
from Bio import SeqIO
def cpg_oe(seq):
  seq=seq.upper().replace("T","U")
  f_C=seq.count("C")/len(seq)
  f_G=seq.count("G")/len(seq)
  f_CpG=seq.count("CG")/(len(seq)-1)
  if f_C*f_G==0:
    return 0
  return f_CpG/(f_C*f_G)
def upa_oe(seq):
  f_U=seq.count("U")/len(seq)
  f_A=seq.count("A")/len(seq)
  f_UpA=seq.count("UA")/(len(seq)-1)
  if f_U*f_A==0:
    return 0
  return f_UpA/(f_U*f_A)
def gc_content(seq):
  seq = seq.upper()
  return(seq.count("G")+seq.count("C")/len(seq))
fasta_file="/content/human_sequences.fasta"
fasta_file="/content/insect_sequences.fasta"
csv_output="human_seq_features.csv"
csv_output="insect_seq_features.csv"
def extract_features(fasta_file,csv_output):
  rows=[]
  for record in SeqIO.parse(fasta_file,"fasta"):
    seq_id=record.id
    sequence=str(record.seq)
    cpg_value=cpg_oe(sequence)
    upa_value=upa_oe(sequence)
    gc_value=gc_content(sequence)
    rows.append([seq_id,cpg_value,upa_value,gc_value])
  with open(csv_output, "w", newline="") as f:
    writer = csv.writer(f)
    writer.writerow(["Sequence_ID", "CpG_OE", "UpA_OE", "GC_content"])
    writer.writerows(rows)
extract_features("/content/human_sequences.fasta","human_seq_features.csv")
extract_features("/content/insect_sequences.fasta","insect_seq_features.csv")
