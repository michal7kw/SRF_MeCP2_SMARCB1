import argparse
import subprocess
import os

def find_common_cpg_peaks(bg_peaks, bm_peaks, cpg_islands, output_file):
    # Create the pipeline command
    cmd = f"""
    bedtools intersect -a {bg_peaks} -b {bm_peaks} -wa | \
    bedtools intersect -a - -b {cpg_islands} -wa | \
    sort -k1,1 -k2,2n > {output_file}
    """
    subprocess.run(cmd, shell=True, check=True)

def main():
    parser = argparse.ArgumentParser(description='Find common peaks in CpG islands')
    parser.add_argument('--bg-peaks', required=True,
                        help='BG merged peaks file')
    parser.add_argument('--bm-peaks', required=True,
                        help='BM merged peaks file')
    parser.add_argument('--cpg-islands', required=True,
                        help='CpG islands bed file')
    parser.add_argument('--output', required=True,
                        help='Output bed file')
    
    args = parser.parse_args()
    
    # Create output directory if it doesn't exist
    os.makedirs(os.path.dirname(args.output), exist_ok=True)
    
    find_common_cpg_peaks(args.bg_peaks, args.bm_peaks, 
                         args.cpg_islands, args.output)

if __name__ == '__main__':
    main() 