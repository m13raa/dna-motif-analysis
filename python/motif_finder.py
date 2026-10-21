import re
from turtle import position

import pandas as pd
from pathlib import Path

input = "data/processed/validated_sequences.csv"
summary_output = "results/motif_summary.csv"
positions_output = "results/motif_positions.csv"

#Define the DNA motifs

motifs = {
    "start_codon_ATG": "ATG",
    "stop_codon_TAA": "TAA",
    "stop_codon_TAG": "TAG",
    "stop_codon_TGA": "TGA",
    "TATA_box": "TATA"
}

#Motif search function 

def find_overlapping_motifs(sequence, motif):
    pattern = re.compile(
        f"(?={re.escape(motif)})"
    )
    positions = [
        match.start()
        for match in pattern.finditer(sequence)
    ]
    return positions

#Load the validated sequences

df = pd.read_csv(input)

#Make sure to exclude any sequences that are not valid DNA sequences

df = df[df["valid_dna"] == True].copy()

summary_results = []
position_results = []

#Analyze the sequences

for index, row in df.iterrows():
    sequence = row["sequence"]
    sequence_id = f"seq_{index + 1:05d}"
    gene_class = row["class"]
    sequence_length = len(sequence)
    for motif_name, motif in motifs.items():
        positions = find_overlapping_motifs(
            sequence,
            motif
        )
        count = len(positions)
        density = (
            count / sequence_length * 1000
            if sequence_length >0
            else 0
        )
        #Store the summary results
        summary_results.append({
            "sequence_id": sequence_id,
            "class": gene_class,
            "sequence_length": sequence_length,
            "gc_content": row["gc_content"],
            "motif_name": motif_name,
            "motif": motif,
            "count": count,
            "density_per_kb": density
        })
        #Individual motif positions
        for position in positions:
            position_results.append({
                "sequence_id": sequence_id,
                "class": gene_class,
                "motif_name": motif_name,
                "motif": motif,
                "start_0based": position,
                "start_1based": position + 1,
                "end_1based": position + len(motif),
                "relative position": (
                    position/sequence_length
                )
            })

#Export the results to CSV files

summary_df = pd.DataFrame(summary_results)
positions_df = pd.DataFrame(
    position_results,
    columns=[
        "sequence_id",
        "class",
        "motif_name",
        "motif",
        "start_0based",
        "start_1based",
        "end_1based",
        "relative position"
    ]
)

Path("results").mkdir(exist_ok=True)

summary_df.to_csv(
    summary_output,
    index=False
)

positions_df.to_csv(
    positions_output,
    index=False
)

#Print the results

print("Motif analysis completed. Results saved to CSV files.")

print("Sequences analyzed:", len(df))

print("Total motif occurrences:", len(positions_df))

print("\nMotif frequencies:")

print(
    summary_df.groupby("motif_name")["count"]
    .sum()
    .sort_values(ascending=False)
)

