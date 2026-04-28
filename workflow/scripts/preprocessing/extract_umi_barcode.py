#!/usr/bin/env python3
"""
Demultiplex FASTQ files by barcode and extract UMI information.
Aligned for single-cell sequencing pipeline.
"""

import gzip
import os
import sys
import argparse
import logging
from pathlib import Path
from typing import Dict, Tuple, Iterator
import json

def setup_logging(log_file: str = None):
    """Set up logging configuration."""
    level = logging.INFO
    format_str = '%(asctime)s - %(name)s - %(levelname)s - %(message)s'
    
    if log_file:
        logging.basicConfig(level=level, format=format_str,
                          handlers=[
                              logging.FileHandler(log_file),
                              logging.StreamHandler(sys.stdout)
                          ])
    else:
        logging.basicConfig(level=level, format=format_str)
    
    return logging.getLogger(__name__)

def parse_args():
    parser = argparse.ArgumentParser(
        description="Demultiplex FASTQ by barcode and extract UMI."
    )
    parser.add_argument('--r1', required=True, help="Input R1 FASTQ (gzipped)")
    parser.add_argument('--r2', required=True, help="Input R2 FASTQ (gzipped)")
    parser.add_argument('--barcodes', required=True, 
                       help="Barcode file (TSV: well_id<TAB>barcode)")
    parser.add_argument('--outdir', required=True, help="Output directory")
    parser.add_argument('--config', help="Config JSON file with parameters")
    parser.add_argument('--batch-size', type=int, default=500000,
                       help="Batch size for writing (default: 500000)")
    parser.add_argument('--log', help="Log file path")
    parser.add_argument('--stats', help="Output statistics JSON file")
    
    # Processing parameters with defaults
    parser.add_argument('--umi-length', type=int, default=12,
                       help="UMI length (default: 12)")
    parser.add_argument('--bc1-length', type=int, default=8,
                       help="Barcode 1 length from R1 (default: 8)")
    parser.add_argument('--bc2-length', type=int, default=8,
                       help="Barcode 2 length from R2 (default: 8)")
    parser.add_argument('--bc1-offset', type=int, default=0,
                       help="Barcode 1 offset in R1 (default: 0)")
    parser.add_argument('--bc2-offset', type=int, default=12,
                       help="Barcode 2 offset in R2 (default: 12)")
    parser.add_argument('--r1-trim', type=int, default=22,
                       help="Bases to trim from R1 5' end (default: 22)")
    parser.add_argument('--r2-trim', type=int, default=34,
                       help="Bases to trim from R2 5' end (default: 34)")
    
    return parser.parse_args()

def load_config(config_file: str) -> Dict:
    """Load configuration from JSON file."""
    if config_file and Path(config_file).exists():
        with open(config_file, 'r') as f:
            return json.load(f)
    return {}

def load_barcodes(barcode_file: str, logger) -> Dict[str, str]:
    """Load barcode mappings from TSV file."""
    barcode_map = {}
    
    with open(barcode_file, 'r') as f:
        # Skip header if present
        header = f.readline().strip()
        if not header.startswith('#') and '\t' in header:
            cols = header.split('\t')
            if cols[0].lower() != 'well_id':
                # Not a header, process as data
                parts = header.split('\t')
                if len(parts) >= 2:
                    barcode_map[parts[1].strip()] = parts[0].strip()
        
        # Process remaining lines
        for line in f:
            if line.strip() and not line.startswith('#'):
                parts = line.strip().split('\t')
                if len(parts) >= 2:
                    well_id, barcode = parts[0], parts[1]
                    barcode_map[barcode.strip()] = well_id.strip()
    
    logger.info(f"Loaded {len(barcode_map)} barcodes")
    return barcode_map

def fastq_reader(handle) -> Iterator[Tuple[str, str, str, str]]:
    """Read FASTQ records as tuples of 4 lines."""
    while True:
        header = handle.readline()
        if not header:
            break
        seq = handle.readline()
        plus = handle.readline()
        qual = handle.readline()
        yield header, seq, plus, qual

def write_batch(well_id: str, batch: list, handles: Dict, stats: Dict):
    """Write a batch of records to output files."""
    if not batch or well_id not in handles:
        return
    
    r1_handle, r2_handle = handles[well_id]
    for r1_record, r2_record in batch:
        r1_handle.write(r1_record)
        r2_handle.write(r2_record)
    
    stats[well_id]['written'] += len(batch)

def process_reads(args, logger):
    """Main processing function."""
    # Create output directory
    outdir = Path(args.outdir)
    outdir.mkdir(parents=True, exist_ok=True)
    
    # Load barcodes
    barcode_map = load_barcodes(args.barcodes, logger)
    
    # Initialize statistics
    stats = {
        'total_reads': 0,
        'assigned_reads': 0,
        'unassigned_reads': 0,
        'wells': {}
    }
    
    # Prepare output handles and buffers
    handles = {}
    batches = {}
    batch_counts = {}
    
    # Create unassigned outputs
    unassigned_r1 = gzip.open(outdir / "unassigned_R1.fastq.gz", 'wt')
    unassigned_r2 = gzip.open(outdir / "unassigned_R2.fastq.gz", 'wt')
    handles['unassigned'] = (unassigned_r1, unassigned_r2)
    batches['unassigned'] = []
    batch_counts['unassigned'] = 0
    stats['wells']['unassigned'] = {'total': 0, 'written': 0}
    
    # Process reads
    logger.info(f"Processing reads from {args.r1} and {args.r2}")
    
    with gzip.open(args.r1, 'rt') as r1_file, gzip.open(args.r2, 'rt') as r2_file:
        for r1_rec, r2_rec in zip(fastq_reader(r1_file), fastq_reader(r2_file)):
            stats['total_reads'] += 1
            
            if stats['total_reads'] % 100000 == 0:
                logger.info(f"Processed {stats['total_reads']:,} reads...")
            
            # Extract barcodes and UMI
            r1_header, r1_seq, r1_plus, r1_qual = r1_rec
            r2_header, r2_seq, r2_plus, r2_qual = r2_rec
            
            # Extract barcode and UMI sequences
            bc1 = r1_seq[args.bc1_offset:args.bc1_offset + args.bc1_length]
            umi = r2_seq[:args.umi_length]
            bc2 = r2_seq[args.bc2_offset:args.bc2_offset + args.bc2_length]
            
            # Trim sequences
            r1_seq_trimmed = r1_seq[args.r1_trim:].rstrip()
            r1_qual_trimmed = r1_qual[args.r1_trim:].rstrip()
            r2_seq_trimmed = r2_seq[args.r2_trim:].rstrip()
            r2_qual_trimmed = r2_qual[args.r2_trim:].rstrip()
            
            # Put UMI at the end of the read name for umi-tools dedup
            read_id = r1_header.split()[0]
            new_header = f"{read_id}_{umi} BC1:{bc1} BC2:{bc2}\n"
            
            # Determine well assignment
            if bc1 in barcode_map:
                well_id = barcode_map[bc1]
                stats['assigned_reads'] += 1
            else:
                well_id = 'unassigned'
                stats['unassigned_reads'] += 1
            
            # Initialize well if needed
            if well_id not in handles:
                r1_out = gzip.open(outdir / f"{well_id}_R1.fastq.gz", 'wt')
                r2_out = gzip.open(outdir / f"{well_id}_R2.fastq.gz", 'wt')
                handles[well_id] = (r1_out, r2_out)
                batches[well_id] = []
                batch_counts[well_id] = 0
                stats['wells'][well_id] = {'total': 0, 'written': 0}
            
            # Update statistics
            stats['wells'][well_id]['total'] += 1
            
            # Format records
            r1_record = f"{new_header}{r1_seq_trimmed}\n{r1_plus}{r1_qual_trimmed}\n"
            r2_record = f"{new_header}{r2_seq_trimmed}\n{r2_plus}{r2_qual_trimmed}\n"
            
            # Add to batch
            batches[well_id].append((r1_record, r2_record))
            batch_counts[well_id] += 1
            
            # Write batch if threshold reached
            if batch_counts[well_id] >= args.batch_size:
                write_batch(well_id, batches[well_id], handles, stats['wells'])
                batches[well_id] = []
                batch_counts[well_id] = 0
    
    # Write remaining batches
    logger.info("Writing remaining batches...")
    for well_id in batches:
        if batches[well_id]:
            write_batch(well_id, batches[well_id], handles, stats['wells'])
    
    # Close all handles
    logger.info("Closing output files...")
    for r1_handle, r2_handle in handles.values():
        r1_handle.close()
        r2_handle.close()
    
    # Log summary statistics
    logger.info("=" * 60)
    logger.info("DEMULTIPLEXING SUMMARY")
    logger.info("=" * 60)
    logger.info(f"Total reads processed: {stats['total_reads']:,}")
    logger.info(f"Assigned reads: {stats['assigned_reads']:,} ({stats['assigned_reads']/stats['total_reads']*100:.1f}%)")
    logger.info(f"Unassigned reads: {stats['unassigned_reads']:,} ({stats['unassigned_reads']/stats['total_reads']*100:.1f}%)")
    logger.info(f"Number of wells: {len(stats['wells']) - 1}")  # Excluding unassigned
    
    # Save statistics
    if args.stats:
        with open(args.stats, 'w') as f:
            json.dump(stats, f, indent=2)
        logger.info(f"Statistics saved to {args.stats}")
    
    return stats

def main():
    args = parse_args()
    
    # Set up logging
    logger = setup_logging(args.log)
    logger.info("Starting demultiplexing process")
    logger.info(f"Parameters: {vars(args)}")
    
    try:
        # Load config if provided
        if args.config:
            config = load_config(args.config)
            # Override args with config values if not specified on command line
            for key, value in config.items():
                if hasattr(args, key) and getattr(args, key) is None:
                    setattr(args, key, value)
        
        # Process reads
        stats = process_reads(args, logger)
        
        logger.info("Demultiplexing completed successfully")
        
    except Exception as e:
        logger.error(f"Error during processing: {str(e)}", exc_info=True)
        sys.exit(1)

if __name__ == "__main__":
    main()
