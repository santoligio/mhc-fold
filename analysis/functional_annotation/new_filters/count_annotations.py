import pandas as pd

def count_gene_name_entries(csv_file):
    # Read CSV file
    df = pd.read_csv(csv_file)
    
    # Check if column exists
    if "gene_name" not in df.columns:
        raise ValueError("Column 'gene_name' not found in the CSV file.")
    
    # Remove missing values
    df = df.dropna(subset=["gene_name"])
    
    # Count occurrences
    counts = df["gene_name"].value_counts().sort_index()
    
    # Print in requested format
    for gene, count in counts.items():
        print(f"{gene}: {count} entries")
    
    # Save to CSV
    output_df = counts.reset_index()
    output_df.columns = ["gene_name", "entries"]
    output_df.to_csv("annotations_number.csv", index=False)


if __name__ == "__main__":
    csv_path = "pdb_mhc_annotations_filtered.csv"
    count_gene_name_entries(csv_path)
