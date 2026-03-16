#!/usr/bin/env python3
"""
Generate test FASTQ data for pipeline testing.
Creates small simulated FASTQ files with known barcodes and UMIs.
"""

import gzip
import random
import argparse
from pathlib import Path
import string

def parse_args():
    parser = argparse.ArgumentParser(
        description="Generate test FASTQ data for pipeline testing"
    )
    parser.add_argument('--output-dir', default='test/data',
                       help='Output directory for test files')
    parser.add_argument('--num-reads', type=int, default=10000,
                       help='Number of reads per file (default: 10000)')
    parser.add_argument('--num-plates', type=int, default=2,
                       help='Number of plates to generate (default: 2)')
    parser.add_argument('--num-wells', type=int, default=4,
                       help='Number of wells to simulate (default: 4)')
    parser.add_argument('--read-length', type=int, default=100,
                       help='Read length (default: 100)')
    parser.add_argument('--umi-length', type=int, default=12,
                       help='UMI length (default: 12)')
    parser.add_argument('--barcode-length', type=int, default=8,
                       help='Barcode length (default: 8)')
    parser.add_argument('--adapter-rate', type=float, default=0.05,
                       help='Proportion of adapter dimers (default: 0.05)')
    parser.add_argument('--duplicate-rate', type=float, default=0.3,
                       help='Proportion of duplicates (default: 0.3)')
    return parser.parse_args()

def generate_random_sequence(length):
    """Generate random DNA sequence."""
    return ''.join(random.choices('ACGT', k=length))

def generate_quality_string(length, min_qual=20, max_qual=40):
    """Generate random quality string."""
    quals = [chr(random.randint(min_qual, max_qual) + 33) for _ in range(length)]
    return ''.join(quals)

def generate_barcodes(num_wells, length=8):
    """Generate unique barcodes for wells."""
    barcodes = []
    well_names = []
    
    for i in range(num_wells):
        row = chr(65 + (i // 12))  # A, B, C, ...
        col = (i % 12) + 1
        well_name = f"{row}{col:02d}"
        well_names.append(well_name)
        
        # Generate unique barcode
        barcode = generate_random_sequence(length)
        while barcode in barcodes:
            barcode = generate_random_sequence(length)
        barcodes.append(barcode)
    
    return dict(zip(well_names, barcodes))

def generate_test_fastq(args, plate_id, barcodes):
    """Generate test FASTQ files for a plate."""
    output_dir = Path(args.output_dir)
    output_dir.mkdir(parents=True, exist_ok=True)
    
    # File paths
    r1_file = output_dir / f"{plate_id}_R1.fastq.gz"
    r2_file = output_dir / f"{plate_id}_R2.fastq.gz"
    
    print(f"Generating {plate_id}...")
    
    # Generate adapter sequence
    adapter = "CAGTCAGCGT"
    
    # Prepare sequence pools for duplicates
    sequence_pool = []
    umi_pool = []
    
    # Pre-generate some sequences for duplicates
    for _ in range(int(args.num_reads * 0.3)):  # 30% of reads
        seq = generate_random_sequence(args.read_length)
        umi = generate_random_sequence(args.umi_length)
        sequence_pool.append(seq)
        umi_pool.append(umi)
    
    # Distribute reads among wells
    well_names = list(barcodes.keys())
    reads_per_well = args.num_reads // len(well_names)
    
    with gzip.open(r1_file, 'wt') as r1, gzip.open(r2_file, 'wt') as r2:
        read_id = 0
        
        for well_name in well_names:
            barcode = barcodes[well_name]
            
            for _ in range(reads_per_well):
                read_id += 1
                
                # Decide if this is a duplicate
                is_duplicate = random.random() < args.duplicate_rate
                is_dimer = random.random() < args.adapter_rate
                
                # Generate or reuse sequences
                if is_duplicate and sequence_pool:
                    insert_seq = random.choice(sequence_pool)
                    umi = random.choice(umi_pool)
                else:
                    insert_seq = generate_random_sequence(args.read_length)
                    umi = generate_random_sequence(args.umi_length)
                
                if is_dimer:
                    # Adapter dimer - starts with adapter
                    insert_seq = adapter + generate_random_sequence(args.read_length - len(adapter))
                
                # R1: barcode(8) + linker(14) + insert
                r1_seq = barcode + "GTCTTGTCTTCTAT" + insert_seq
                r1_qual = generate_quality_string(len(r1_seq))
                
                # R2: UMI(12) + barcode(8) + linker(14) + insert  
                r2_seq = umi + barcode + "AGATCGGAAGAGCA" + insert_seq
                if is_dimer:
                    r2_seq = umi + barcode + "AGATCGGAAGAGCA" + adapter + generate_random_sequence(args.read_length - len(adapter))
                r2_qual = generate_quality_string(len(r2_seq))
                
                # Write to files
                header = f"@SEQ{read_id:08d} {well_name}\n"
                
                r1.write(header)
                r1.write(r1_seq + "\n")
                r1.write("+\n")
                r1.write(r1_qual + "\n")
                
                r2.write(header)
                r2.write(r2_seq + "\n")
                r2.write("+\n")
                r2.write(r2_qual + "\n")
    
    print(f"  Created {r1_file}")
    print(f"  Created {r2_file}")
    
    return r1_file, r2_file

def create_test_config(args, barcodes, plate_files):
    """Create test configuration files."""
    output_dir = Path(args.output_dir)
    
    # Create barcode file
    barcode_file = output_dir / "test_barcodes.tsv"
    with open(barcode_file, 'w') as f:
        f.write("well_id\tbarcode\n")
        for well, bc in barcodes.items():
            f.write(f"{well}\t{bc}\n")
    print(f"Created {barcode_file}")
    
    # Create sample sheet
    sample_file = output_dir / "test_samples.tsv"
    with open(sample_file, 'w') as f:
        f.write("plate_id\tr1_fastq\tr2_fastq\tsample_group\ttreatment\tnotes\n")
        for plate_id, (r1, r2) in plate_files.items():
            f.write(f"{plate_id}\t{r1}\t{r2}\tTestSample\tcontrol\tTest data\n")
    print(f"Created {sample_file}")
    
    # Create test config
    config_file = output_dir / "test_config.yaml"
    config_content = f"""# Test configuration for single-cell pipeline

project_name: "Test_Run"
output_dir: "test_results"

samples_file: "{sample_file}"
barcodes_file: "{barcode_file}"

genome:
  fasta: "resources/genome/test_genome.fa"
  index_prefix: "resources/genome/test_genome"
  assembly: "test"
  chromosomes: ["chr1", "chr2", "chr3"]

preprocessing:
  umi_length: {args.umi_length}
  barcode1_length: {args.barcode_length}
  barcode2_length: {args.barcode_length}
  barcode1_offset: 0
  barcode2_offset: {args.umi_length}
  r1_trim_5prime: 22
  r2_trim_5prime: 34
  dedup_method: "umi_insert"
  adapter_sequence: "CAGTCAGCGT"
  case_insensitive: false

alignment:
  bowtie2_params: "--very-sensitive"
  threads: 4
  memory: "8G"

aneufinder:
  method: ["edivisive"]
  binsize: 1000000
  num_cpu: 4
  chromosomes_to_analyze: ["chr1", "chr2", "chr3"]

resources:
  default:
    threads: 1
    mem_mb: 2000
    time: "00:30:00"
  demux:
    threads: 2
    mem_mb: 4000
    time: "01:00:00"
  alignment:
    threads: 4
    mem_mb: 8000
    time: "02:00:00"

batch_size: 100000
"""
    
    with open(config_file, 'w') as f:
        f.write(config_content)
    print(f"Created {config_file}")
    
    return config_file

def create_test_genome(output_dir):
    """Create a small test reference genome."""
    genome_dir = output_dir / ".." / ".." / "resources" / "genome"
    genome_dir.mkdir(parents=True, exist_ok=True)
    
    genome_file = genome_dir / "test_genome.fa"
    
    if not genome_file.exists():
        print("Creating test reference genome...")
        with open(genome_file, 'w') as f:
            # Create 3 small chromosomes
            for i in range(1, 4):
                f.write(f">chr{i}\n")
                # 10kb chromosomes for testing
                seq = generate_random_sequence(10000)
                # Write in 80bp lines
                for j in range(0, len(seq), 80):
                    f.write(seq[j:j+80] + "\n")
        print(f"  Created {genome_file}")
    
    return genome_file

def main():
    args = parse_args()
    
    print("=" * 60)
    print("GENERATING TEST DATA")
    print("=" * 60)
    
    # Generate barcodes
    print(f"\nGenerating {args.num_wells} well barcodes...")
    barcodes = generate_barcodes(args.num_wells, args.barcode_length)
    
    # Generate FASTQ files for each plate
    plate_files = {}
    for i in range(1, args.num_plates + 1):
        plate_id = f"test_plate{i}"
        r1, r2 = generate_test_fastq(args, plate_id, barcodes)
        plate_files[plate_id] = (r1, r2)
    
    # Create configuration files
    print("\nCreating configuration files...")
    config_file = create_test_config(args, barcodes, plate_files)
    
    # Create test genome
    genome_file = create_test_genome(Path(args.output_dir))
    
    print("\n" + "=" * 60)
    print("TEST DATA GENERATION COMPLETE")
    print("=" * 60)
    print("\nTo run the test pipeline:")
    print(f"  1. Build genome index:")
    print(f"     bowtie2-build {genome_file} resources/genome/test_genome")
    print(f"  2. Run pipeline:")
    print(f"     snakemake --cores 4 --use-conda --configfile {config_file}")
    print("\nNote: Test data is small and may not produce meaningful")
    print("      biological results, but will test pipeline functionality.")

if __name__ == "__main__":
    main()
