import pandas as pd

INPUT = "data/raw/human.txt"
OUTPUT = "data/processed/validated_sequences.csv"

#Load the dataset
df = pd.read_table(INPUT)

#Standardize sequences
print("dataset loaded")
print("Number of sequences:", len(df))
df["sequence"] = df["sequence"].astype(str).str.upper()

#Define the valid bases

valid_bases = set("ATGC")

def is_valid_dna(sequence):
    """
    Return True if a sequence contains
    only A, T, G, and C (valid bases)"""

    return set(sequence).issubset(valid_bases)

#Calculate the GC content

def gc_content(sequence):
    """
    Calculate percentage GC content.
    """
    if len(sequence) == 0:
        return 0
    gc_count = (
        sequence.count("G")
        + sequence.count("C")
    )
    return (gc_count / len(sequence)) * 100

#Generate QC metrics

df["length"] = df["sequence"].str.len()

df["valid_dna"] = df["sequence"].apply(is_valid_dna)

df["gc_content"] = df["sequence"].apply(gc_content)

#Print QC summary

print("\nDNA validation:")

print(
    "Valid sequences:",
    df["valid_dna"].sum()
    )

print(
    "Invalid sequences:",
    (~df["valid_dna"]).sum()
)

print("\nSequence length summmary:")

print(
    df["length"].describe()
)

print("\nGC content summary:")

print(
    df["gc_content"].describe()
)

#Save the processed dataset

df.to_csv(
    OUTPUT,
    index=False
)

print(
    "\nProcessed dataset saved to:",
    OUTPUT
)