import argparse
import subprocess
import os

def merge_peaks(input_peaks, output_file):
    # Concatenate all peak files and pipe to sort, then to bedtools merge
    cat_cmd = f"cat {' '.join(input_peaks)} | sort -k1,1 -k2,2n | bedtools merge > {output_file}"
    subprocess.run(cat_cmd, shell=True, check=True)

def main():
    parser = argparse.ArgumentParser(description='Merge peak files for a group of samples')
    parser.add_argument('--input-peaks', nargs='+', required=True,
                        help='List of peak files to merge')
    parser.add_argument('--output', required=True,
                        help='Output merged bed file')
    
    args = parser.parse_args()
    
    # Create output directory if it doesn't exist
    os.makedirs(os.path.dirname(args.output), exist_ok=True)
    
    merge_peaks(args.input_peaks, args.output)

if __name__ == '__main__':
    main() 