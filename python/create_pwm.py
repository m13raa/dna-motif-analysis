#Create a position weight matrix

from pathlib import Path

import numpy as np
import pandas as pd

#Define the example aligned sequences

sequences = [
    "ATGA",
    "ATGT",
    "ACGA",
    "ATGA",
    "ATGG"
]

nucleotides = ["A", "C", "G", "T"]

motif_length = len(sequences[0])

#Check to see if all sequences have equal length

assert all(
    len(seq) == motif_length
    for seq in sequences
), "All sequences must have equal length."

#Calculate the nucleotide counts

count_matrix = pd.DataFrame(
    0,
    index = nucleotides,
    columns=range(1, motif_length +1)
)

for seq in sequences:
    for position, nucleotide in enumerate(seq, start=1):
        count_matrix.loc[nucleotide, position] += 1

print("\nPosition Count Matrix:")
print(count_matrix)

#Conver the counts into probabilities 

#Create pseudocount the prevent zero probabilities
pseudocount = 1

probability_matrix = (
    count_matrix + pseudocount
) / (
    len(sequences) + pseudocount * len(nucleotides)
)

print("\nPosition Probability Matrix:")
print(probability_matrix)

#Calculate the log-dds scores

#Assume equal background nucleotide probablities
background_probability = 0.25

pwm= np.log2(
    probability_matrix / background_probability
)

print("\nPosition Weght Matrix:")
print(pwm.round(3))

#Save results

Path("results").mkdir(exist_ok=True)

count_matrix.to_csv(
    "results/pwm_counts.csv"
)

probability_matrix.to_csv(
    "results/pwm_probabilities.csv"
)

pwm.to_csv(
    "results/pwm_log_odds.csv"
)

print("\nPWM files saved successfully!")