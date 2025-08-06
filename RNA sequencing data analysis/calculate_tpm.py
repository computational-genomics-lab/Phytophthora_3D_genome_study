import pandas as pd
import re
import os

# --- Configuration ---
# !!! IMPORTANT: Update these filenames if they differ slightly !!!
GTF_FILE = 'Phytophthora_capsici_WB2_name_parse_gene_name_updated_fixed_Hic.gtf'
COUNTS_FILE = 'featureCounts_marge_0_16_48h_filt10reads_myc_early_late_infection_filt100reads.csv'
OUTPUT_TPM_FILE = 'TPM_normalized_counts_1.csv'

print(f"Starting TPM calculation...")
print(f"GTF file specified: {GTF_FILE}")
print(f"Counts file specified: {COUNTS_FILE}")

# --- Step 1: Parse GTF file to calculate gene lengths ---
print("\nStep 1: Calculating gene lengths from GTF...")
try:
    # Read the GTF file
    # GTF files typically have columns separated by tabs, and comments start with '#'
    gtf_df = pd.read_csv(
        GTF_FILE,
        sep='\t',
        header=None,
        comment='#'
    )
    # Assign column names for easier access
    gtf_df.columns = [
        'seqname', 'source', 'feature', 'start', 'end',
        'score', 'strand', 'frame', 'attribute'
    ]
    print(f"Successfully loaded GTF file with {len(gtf_df)} entries.")
except FileNotFoundError:
    print(f"Error: GTF file '{GTF_FILE}' not found.")
    print("Please ensure the file is in the same directory as the script and the name is exactly correct.")
    exit()
except Exception as e:
    print(f"An unexpected error occurred while reading the GTF file: {e}")
    exit()


# Function to extract gene_id from the 'attribute' column
# This regex looks for 'gene_id "your_gene_id";'
def extract_gene_id(attributes):
    match = re.search(r'gene_id "([^"]+)"', attributes)
    return match.group(1) if match else None

gtf_df['gene_id'] = gtf_df['attribute'].apply(extract_gene_id)

# Filter for 'exon' features and calculate length for each exon
# Gene length for TPM is usually the sum of all its exon lengths.
exon_df = gtf_df[gtf_df['feature'] == 'exon'].copy()
exon_df['length'] = exon_df['end'] - exon_df['start'] + 1 # Length is end - start + 1 (inclusive)

# Group by gene_id and sum exon lengths to get total gene length per gene
gene_lengths = exon_df.groupby('gene_id')['length'].sum().reset_index()
gene_lengths.rename(columns={'length': 'GeneLength'}, inplace=True)

print(f"Calculated lengths for {len(gene_lengths)} unique genes from GTF.")
# print("\nPreview of Gene Lengths:")
# print(gene_lengths.head())

# --- Step 2: Load featureCounts data ---
print("\nStep 2: Loading featureCounts data...")
try:
    # Read the CSV file. The first column usually contains gene IDs, so set it as index.
    counts_df = pd.read_csv(COUNTS_FILE, index_col=0)
    print(f"Successfully loaded counts for {counts_df.shape[0]} genes and {counts_df.shape[1]} samples.")
except FileNotFoundError:
    print(f"Error: Counts file '{COUNTS_FILE}' not found.")
    print("Please ensure the file is in the same directory as the script and the name is exactly correct.")
    exit()
except Exception as e:
    print(f"An unexpected error occurred while reading the counts file: {e}")
    exit()

# print("\nPreview of Raw Counts:")
# print(counts_df.head())

# --- Step 3: Merge gene lengths with counts data ---
print("\nStep 3: Merging gene lengths with counts data...")
# Merge the counts DataFrame with the gene lengths DataFrame based on gene_id
# 'left_index=True' means use the index of counts_df (gene IDs) for merging
# 'right_on='gene_id'' means use the 'gene_id' column from gene_lengths for merging
# 'how='left'' ensures all genes from the counts_df are kept
merged_df = pd.merge(counts_df, gene_lengths, left_index=True, right_on='gene_id', how='left')
merged_df.set_index('gene_id', inplace=True) # Set gene_id back as the index

# Handle genes that might be in counts but not in GTF (i.e., no length information)
# These genes will have NaN in the 'GeneLength' column after the merge.
if merged_df['GeneLength'].isnull().any():
    missing_lengths_count = merged_df['GeneLength'].isnull().sum()
    print(f"Warning: {missing_lengths_count} genes in the counts file do not have corresponding length information in the GTF.")
    print("These genes will be excluded from TPM calculation as their lengths are unknown.")
    merged_df.dropna(subset=['GeneLength'], inplace=True) # Remove rows with NaN in 'GeneLength'
    print(f"Proceeding with {merged_df.shape[0]} genes after filtering for missing lengths.")

# Identify sample columns (all columns except the newly added 'GeneLength')
# We take the original columns from counts_df to ensure we only get sample names.
sample_columns = counts_df.columns.tolist()

# Convert gene lengths from base pairs to kilobases for RPK calculation
merged_df['GeneLength_kb'] = merged_df['GeneLength'] / 1000

# --- Step 4: Calculate TPM ---
print("\nStep 4: Calculating TPM values...")
tpm_df = pd.DataFrame(index=merged_df.index) # Create a new DataFrame for TPM results

for col in sample_columns:
    # Calculate RPK (Reads Per Kilobase) for each gene in the current sample
    # merged_df[col] are the raw counts for the current sample
    # merged_df['GeneLength_kb'] are the gene lengths (in kilobases)
    rpk = merged_df[col] / merged_df['GeneLength_kb']

    # Calculate the scaling factor for the current sample: sum of all RPKs for that sample
    scaling_factor = rpk.sum()

    if scaling_factor == 0:
        # If the sum of RPKs for a sample is zero (e.g., no reads, very short genes),
        # TPM for that sample will also be zero to avoid division by zero.
        tpm_df[f'TPM_{col}'] = 0
        print(f"  Warning: Total RPK for sample '{col}' is zero. TPMs for this sample will be 0.")
    else:
        # Calculate TPM: (RPK / Scaling Factor) * 1,000,000
        tpm_df[f'TPM_{col}'] = (rpk / scaling_factor) * 1e6

print("TPM calculation complete for all samples.")

# --- Step 5: Save the TPM results ---
print(f"\nStep 5: Saving TPM data to '{OUTPUT_TPM_FILE}'...")
tpm_df.to_csv(OUTPUT_TPM_FILE)

print("\n--- Calculation Finished ---")
print(f"TPM data successfully saved to: {OUTPUT_TPM_FILE}")
print("\nFirst 10 rows of calculated TPM data:")
print(tpm_df.head(10)) # Display more rows for a better preview
print("\nSummary of TPM values (e.g., sum per sample - should be close to 1,000,000):")
print(tpm_df.sum()) # Sum of TPMs for each sample should ideally be 1,000,000
