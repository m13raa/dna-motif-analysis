from pathlib import Path
import re
import time
import pandas as pd

#Load the validated DNA sequence data from the csv file

df = pd.read_csv(
    "data/processed/validated_sequences.csv"
)

df = df[df["valid_dna"] == True].copy()

sequences = df["sequence"].tolist()

motif = "ATG"

#Implement a sliding-window motif search algorithm to find the motif in the DNA sequences

def sliding_window_count(sequence, motif):
    count = 0
    motif_length = len(motif)
    for i in range(
        len(sequence) - motif_length + 1
    ):
        if sequence[i:i + motif_length] == motif:
            count += 1
    return count

#Implement regular expression

def regex_count(sequence,motif):
    pattern = re.compile(
        rf"(?={re.escape(motif)})"
    )

    return len(pattern.findall(sequence))

#Confirm that both the methods agree

for sequence in sequences:
    sliding_count = sliding_window_count(
        sequence,
        motif
    )
    regex_result = regex_count(
        sequence,
        motif
    )
    assert sliding_count == regex_result, (
        "Motif counts do not match!"
    )
print("Both algorithms produce identical counts.")

#Benchmark the performance of both algorithmss

def benchmark(function, sequences, motif, repeats=5):
    durations = []
    for _ in range(repeats):
        start = time.perf_counter()
        for sequence in sequences:
            function(sequence, motif)
        end = time.perf_counter()
        durations.append(end - start)
    return durations

sliding_times = benchmark(
    sliding_window_count,
    sequences,
    motif
)

regex_times = benchmark(
    regex_count,
    sequences,
    motif
)

#Save the results

results = pd.DataFrame({
    "method": (
        ["Sliding Window"] * len(sliding_times)
        + ["Regular Expression"] * len(regex_times)
    ),
    "time_seconds": sliding_times + regex_times
})

Path("results").mkdir(exist_ok=True)

results.to_csv(
    "results/motif_benchmark.csv",
    index=False
)

print("\nBenchmark Results:")
print(
    results.groupby("method")["time_seconds"]
    .agg(["mean", "std", "min", "max"])
)