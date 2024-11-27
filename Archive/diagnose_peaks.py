import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns
import os

def diagnose_peak_data(peak_comparison_file, output_dir):
    """Analyze and visualize peak comparison data."""
    
    # Read the data
    df = pd.read_csv(peak_comparison_file, sep='\t')
    
    # Basic statistics
    stats = {
        'Total peaks': len(df),
        'Unique chromosomes': df['chr'].nunique(),
        'Missing BG values': df['bg_mean'].isna().sum(),
        'Missing BM values': df['bm_mean'].isna().sum(),
        'Zero BG values': (df['bg_mean'] == 0).sum(),
        'Zero BM values': (df['bm_mean'] == 0).sum()
    }
    
    # Create diagnostic plots
    plt.figure(figsize=(15, 10))
    
    # 1. Distribution of peak intensities
    plt.subplot(2, 2, 1)
    sns.histplot(data=df, x='bm_mean', bins=50)
    plt.title('Distribution of BM Peak Intensities')
    plt.xlabel('Peak Intensity')
    
    # 2. Peak lengths
    plt.subplot(2, 2, 2)
    df['peak_length'] = df['end'] - df['start']
    sns.histplot(data=df, x='peak_length', bins=50)
    plt.title('Distribution of Peak Lengths')
    plt.xlabel('Peak Length (bp)')
    
    # 3. Peaks per chromosome
    plt.subplot(2, 2, 3)
    df['chr'].value_counts().plot(kind='bar')
    plt.title('Peaks per Chromosome')
    plt.xticks(rotation=45)
    
    # Save results
    plt.tight_layout()
    plt.savefig(os.path.join(output_dir, 'peak_diagnostics.png'))
    
    # Save statistics
    with open(os.path.join(output_dir, 'peak_statistics.txt'), 'w') as f:
        for key, value in stats.items():
            f.write(f"{key}: {value}\n")
    
    return stats

def main():
    import argparse
    parser = argparse.ArgumentParser()
    parser.add_argument('--input', required=True, help='Peak comparison file')
    parser.add_argument('--output-dir', required=True, help='Output directory')
    args = parser.parse_args()
    
    os.makedirs(args.output_dir, exist_ok=True)
    diagnose_peak_data(args.input, args.output_dir)

if __name__ == '__main__':
    main() 