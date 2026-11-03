
#!/usr/bin/env python3

import argparse
import pandas as pd

parser = argparse.ArgumentParser()
parser.add_argument("--input", required=True)
parser.add_argument("--output", default="sample_variant_results_longer.tsv")
args = parser.parse_args()

# Read TSV, preserving "." and empty values
data = pd.read_csv(
    args.input,
    sep="\t",
    keep_default_na=False
)

# Convert binary annotation columns to booleans
boolean_columns = [
    "IsPseudo",
    "IsGenic",
    "IsTransition",
    "IsSynonymous"
]

for col in boolean_columns:
    data[col] = data[col].replace({
        1: 'TRUE',
        0: 'FALSE',
        9: 'NA'
    })

# First 19 columns are identifiers
# TODO: Parameterize ID vs Sample columns
id_cols = data.columns[:19].tolist()

# Remaining columns represent individual samples
value_cols = data.columns[19:].tolist()

# Convert wide to long
cdata = data.melt(
    id_vars=id_cols,
    value_vars=value_cols,
    var_name="Sample",
    value_name="Percentage"
)

# Remove missing variant observations
cdata = cdata[cdata["Percentage"] != "."]

# Move Sample and Percentage to the beginning
cdata = cdata[
    ["Sample", "Percentage"] +
    [col for col in cdata.columns if col not in ["Sample", "Percentage"]]
]

# Write to TSV
cdata.to_csv(args.output, sep = "\t", index=False)
