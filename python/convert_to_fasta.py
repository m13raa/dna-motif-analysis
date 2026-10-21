import pandas as pd

from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
from Bio import SeqIO

INPUT = "data/processed/validated_sequences.csv"
OUTPUT = "data/processed/human_sequences.fasta"

#Load the processed sequences

df = pd.read_csv(INPUT)

records = []

for index, row in df.iterrows():

    if not row["valid_dna"]:
        continue

    record = SeqRecord(
        Seq(row["sequence"]),
        id=f"seq_{index + 1:05d}",
        description=f"class={row['class']}"
    )

    records.append(record)

SeqIO.write(records, OUTPUT, "fasta")

print(f"Wrote {len(records)} sequences")
print(f"Saved to {OUTPUT}")