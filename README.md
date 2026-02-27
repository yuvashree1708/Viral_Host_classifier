🧠 Background

Viruses experience different selective pressures depending on the host they infect. One well-known mechanism of host antiviral defense is nucleic acid–based immunity, which can act on specific dinucleotide patterns.
In vertebrates:

CpG dinucleotides are often suppressed

Proteins like ZAP preferentially target CpG-rich viral RNA

In many RNA viruses:

CpG and UpA dinucleotides are underrepresented

This suppression may reflect host adaptation

In a previous project analyzing mutational frequencies of viral proteins across time, we observed unusual stabilization of certain amino acid positions during outbreak years.
This raised a broader question:

Could viral genomes be adapting at the nucleotide level in response to host-specific immune pressures?
🎯 Objective

To compare CpG observed/expected (O/E) ratios, UpA O/E ratios, and GC content of:

CHIKV sequences isolated from human hosts

CHIKV sequences isolated from insect hosts

We aimed to test whether:

There is measurable host-dependent dinucleotide adaptation.

Outbreak-associated viral populations show nucleotide-level immune adaptation signals.
⚙️ Methodology
1️⃣ Input Data

FASTA sequences of CHIKV

Separated by host:

human_sequences.fasta

insect_sequences.fasta

2️⃣ Features Computed

For each sequence:
<img width="772" height="345" alt="image" src="https://github.com/user-attachments/assets/d08f55f0-521f-4cb0-8ccf-381b44d98241" />
Where:

L = sequence length

Frequencies normalized per sequence

3️⃣ Implementation

The pipeline:

Parses FASTA using Biopython

Computes CpG O/E, UpA O/E, and GC content

Exports feature tables to CSV

Output files:

human_seq_features.csv

insect_seq_features.csv


