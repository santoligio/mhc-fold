#!/usr/bin/env python3

import pandas as pd
import matplotlib.pyplot as plt
import numpy as np

# ------------------------------------------------------------
# Input file
# ------------------------------------------------------------
INPUT_FILE = "pdb_mhc_annotations_filtered.csv"

# ------------------------------------------------------------
# Load data
# ------------------------------------------------------------
df = pd.read_csv(INPUT_FILE)

if "target_length" not in df.columns:
    raise ValueError("Column 'target_length' not found in the CSV file.")

# Convert to numeric (invalid values like "-" become NaN)
values = pd.to_numeric(df["target_length"], errors="coerce").dropna()

# ------------------------------------------------------------
# Plot normalized histogram (%)
# ------------------------------------------------------------
plt.figure()

# Weights so total = 100%
weights = np.ones(len(values)) / len(values) * 100

plt.hist(values, bins=30, weights=weights)

plt.xlabel("Target Sequence Length")
plt.ylabel("Percentage (%)")

# Set x-axis ticks every 10 units
plt.xticks(np.arange(130, 186, 5))
plt.xlim(130, 185)
plt.tight_layout()

plt.savefig("pdb_histogram_target_length", dpi=300)
plt.close()
