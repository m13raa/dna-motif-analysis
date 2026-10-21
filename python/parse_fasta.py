from Bio import SeqIO

#Path to the FASTA file
FASTA = "data/processed/human_sequences.fasta"

#Initialize counters for the number of sequences and total bases
number_sequences = 0
total_bases = 0

#Read the fasta file
for record in SeqIO.parse(FASTA, "fasta"):
    number_sequences += 1
    total_bases += len(record.seq)

#Results 
print("Sequences:", number_sequences)
print("Total_bases:", total_bases)

if number_sequences >0:
    print("Average length:", total_bases / number_sequences)

    