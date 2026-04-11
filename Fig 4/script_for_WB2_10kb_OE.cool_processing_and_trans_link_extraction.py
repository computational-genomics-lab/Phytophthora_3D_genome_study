import cooler
import numpy as np
import pandas as pd

# -----------------------------
# Parameters
# -----------------------------
cool_file = "WB2_10kb_OE.cool"
resolution = 10000
threshold = 1   # minimum OE value

# -----------------------------
# Load cooler file
# -----------------------------
c = cooler.Cooler(cool_file)

chromosomes = c.chromnames

results = []

print("Chromosomes detected:")
print(chromosomes)

# -----------------------------
# Scan chromosome pairs
# -----------------------------
for i, chr1 in enumerate(chromosomes):

    for j, chr2 in enumerate(chromosomes):

        # Only inter-chromosomal
        if i >= j:
            continue

        print(f"Processing {chr1} vs {chr2}")

        # Fetch interaction block
        block = c.matrix(balance=False).fetch(chr1, chr2)

        # Find bins above threshold
        rows, cols = np.where(block >= threshold)

        print(f"Interactions found: {len(rows)}")

        # Convert bin indices to genomic coordinates
        for r, c2 in zip(rows, cols):

            start1 = r * resolution
            end1 = start1 + resolution

            start2 = c2 * resolution
            end2 = start2 + resolution

            oe_value = float(block[r, c2])

            # log2(O/E)
            log2_oe = np.log2(oe_value)

            # -------- ADDED FILTER --------
            if log2_oe < 5:
                continue
            # --------------------------------

            results.append([
                chr1, start1, end1,
                chr2, start2, end2,
                oe_value,
                log2_oe
            ])

# -----------------------------
# Create dataframe
# -----------------------------
df = pd.DataFrame(
    results,
    columns=[
        "Chromosome1",
        "Start1",
        "End1",
        "Chromosome2",
        "Start2",
        "End2",
        "OE_value",
        "log2OE"
    ]
)

print(df.head())

# -----------------------------
# Save results
# -----------------------------
output_file = "interchromosomal_links.tsv"

df.to_csv(
    output_file,
    sep="\t",
    index=False
)

print(f"\nSaved to: {output_file}")
print(f"Total interactions detected: {len(df)}")
