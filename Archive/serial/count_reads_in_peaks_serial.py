import argparse
import subprocess
import os

def count_reads(peaks_file, bam_file, output_file):
    # Check if input files exist
    if not os.path.exists(peaks_file):
        raise FileNotFoundError(f"Peaks file not found: {peaks_file}")
    if not os.path.exists(bam_file):
        raise FileNotFoundError(f"BAM file not found: {bam_file}")
    
    # First, let's check the chromosome names in both files
    print(f"Checking chromosome names in {peaks_file}...")
    try:
        peaks_chroms = subprocess.check_output(f"cut -f1 {peaks_file} | sort | uniq", 
                                             shell=True).decode().strip().split('\n')
        print(f"Peaks chromosomes: {peaks_chroms}")
    except subprocess.CalledProcessError as e:
        print(f"Error checking peaks file: {str(e)}")
        raise

    print(f"Checking chromosome names in {bam_file}...")
    try:
        bam_chroms = subprocess.check_output(f"samtools view -H {bam_file} | grep '^@SQ' | cut -f2 | cut -d':' -f2", 
                                           shell=True).decode().strip().split('\n')
        print(f"BAM chromosomes: {bam_chroms}")
    except subprocess.CalledProcessError as e:
        print(f"Error checking BAM file: {str(e)}")
        raise

    # Run bedtools coverage with detailed error output
    cmd = f"bedtools coverage -sorted -split -a {peaks_file} -b {bam_file} -counts > {output_file}"
    print(f"Running command: {cmd}")
    
    try:
        result = subprocess.run(cmd, shell=True, check=False, 
                              stderr=subprocess.PIPE, stdout=subprocess.PIPE)
        if result.returncode == 134:
            raise MemoryError("bedtools coverage ran out of memory. Try sorting input files or using less memory-intensive options")
        elif result.returncode != 0:
            print(f"Command failed with return code: {result.returncode}")
            print(f"stderr: {result.stderr.decode()}")
            print(f"stdout: {result.stdout.decode()}")
            raise subprocess.CalledProcessError(result.returncode, cmd, 
                                             stderr=result.stderr)
    except Exception as e:
        print(f"Error running bedtools coverage: {str(e)}")
        raise

def main():
    parser = argparse.ArgumentParser(description='Count reads in peaks')
    parser.add_argument('--peaks', required=True,
                        help='Peaks bed file')
    parser.add_argument('--bam', required=True,
                        help='BAM file')
    parser.add_argument('--output', required=True,
                        help='Output counts file')
    
    args = parser.parse_args()
    
    # Create output directory if it doesn't exist
    os.makedirs(os.path.dirname(args.output), exist_ok=True)
    
    try:
        count_reads(args.peaks, args.bam, args.output)
    except Exception as e:
        print(f"Error processing {args.bam}: {str(e)}")
        raise

if __name__ == '__main__':
    main() 