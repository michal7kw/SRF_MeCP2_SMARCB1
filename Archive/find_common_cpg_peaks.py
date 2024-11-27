import argparse
import subprocess
import os

def find_common_cpg_peaks(bg_peaks, bm_peaks, cpg_islands, output_file, threads=1):
    # Create temporary files
    bg_cpg = bg_peaks + ".cpg.tmp"
    bm_cpg = bm_peaks + ".cpg.tmp"
    
    try:
        # First find peaks that overlap with CpG islands for each sample
        print("Finding peaks overlapping with CpG islands...")
        subprocess.run(f"bedtools intersect -a {bg_peaks} -b {cpg_islands} -u > {bg_cpg}",
                     shell=True, check=True)
        subprocess.run(f"bedtools intersect -a {bm_peaks} -b {cpg_islands} -u > {bm_cpg}",
                     shell=True, check=True)
        
        # Count peaks in CpG islands for each sample
        bg_count = int(subprocess.check_output(f"wc -l {bg_cpg}", shell=True).decode().split()[0])
        bm_count = int(subprocess.check_output(f"wc -l {bm_cpg}", shell=True).decode().split()[0])
        print(f"Found {bg_count} BG peaks and {bm_count} BM peaks in CpG islands")
        
        # Find peaks that overlap between samples with minimum overlap of 1bp
        print("Finding overlapping peaks between samples...")
        cmd = f"""
        bedtools intersect -a {bg_cpg} -b {bm_cpg} -u | \
        sort -k1,1 -k2,2n > {output_file}
        """
        subprocess.run(cmd, shell=True, check=True)
        
        # Count final peaks
        final_count = int(subprocess.check_output(f"wc -l {output_file}", shell=True).decode().split()[0])
        print(f"Found {final_count} common peaks in CpG islands")
        
    finally:
        # Cleanup temporary files
        for f in [bg_cpg, bm_cpg]:
            if os.path.exists(f):
                os.remove(f)

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
    parser.add_argument('--threads', type=int, default=1,
                       help='Number of threads to use')
    
    args = parser.parse_args()
    
    # Create output directory if it doesn't exist
    os.makedirs(os.path.dirname(args.output), exist_ok=True)
    
    find_common_cpg_peaks(args.bg_peaks, args.bm_peaks, 
                         args.cpg_islands, args.output, threads=args.threads)

if __name__ == '__main__':
    main() 