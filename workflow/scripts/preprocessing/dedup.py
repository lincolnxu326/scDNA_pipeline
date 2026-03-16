#!/usr/bin/env python3
"""
Deduplicate demultiplexed FASTQ files by UMI and/or insert sequence.
Aligned for single-cell sequencing pipeline.
"""

import gzip
import os
import sys
import argparse
import logging
from pathlib import Path
from typing import Dict, Set, Tuple
import json
from collections import defaultdict

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
        description="Deduplicate demultiplexed FASTQ files by UMI and/or insert."
    )
    parser.add_argument('--indir', required=True, 
                       help="Input directory of demultiplexed FASTQs")
    parser.add_argument('--outdir', required=True, 
                       help="Output directory for deduplicated FASTQs")
    parser.add_argument('--method', default='umi_insert',
                       choices=['umi_only', 'insert_only', 'umi_insert'],
                       help="Deduplication method (default: umi_insert)")
    parser.add_argument('--suffix', default='.fastq.gz',
                       help="FASTQ file suffix (default: .fastq.gz)")
    parser.add_argument('--log', help="Log file path")
    parser.add_argument('--stats', help="Output statistics JSON file")
    parser.add_argument('--wells', nargs='*', 
                       help="Specific wells to process (default: all)")
    parser.add_argument('--min-length', type=int, default=20,
                       help="Minimum read length after trimming (default: 20)")
    
    return parser.parse_args()

def extract_umi_from_header(header: str) -> str:
    """Extract UMI from the annotated header."""
    if '|UMI:' in header:
        umi = header.split('|UMI:')[1].split('|')[0].strip()
        return umi
    return 'noUMI'

def get_dedup_key(r1_lines: list, r2_lines: list, method: str) -> Tuple:
    """
    Generate deduplication key based on method.
    
    Args:
        r1_lines: Four lines from R1 FASTQ record
        r2_lines: Four lines from R2 FASTQ record
        method: Deduplication method
    
    Returns:
        Tuple key for deduplication
    """
    umi = extract_umi_from_header(r1_lines[0])
    
    if method == 'umi_only':
        return (umi,)
    elif method == 'insert_only':
        # Use concatenated insert sequences from both reads
        insert = r1_lines[1].strip() + r2_lines[1].strip()
        return (insert,)
    else:  # umi_insert
        # Use R1 insert for simplicity, could be modified
        insert = r1_lines[1].strip()
        return (umi, insert)

def process_well(well_id: str, indir: Path, outdir: Path, 
                method: str, min_length: int, 
                suffix: str, logger) -> Dict:
    """Process a single well for deduplication."""
    
    r1_in = indir / f"{well_id}_R1{suffix}"
    r2_in = indir / f"{well_id}_R2{suffix}"
    r1_out = outdir / f"{well_id}_R1.dedup{suffix}"
    r2_out = outdir / f"{well_id}_R2.dedup{suffix}"
    
    if not r1_in.exists() or not r2_in.exists():
        logger.warning(f"Input files not found for {well_id}, skipping")
        return None
    
    logger.info(f"Processing {well_id}...")
    
    seen = set()
    stats = {
        'well_id': well_id,
        'total_reads': 0,
        'unique_reads': 0,
        'duplicate_reads': 0,
        'filtered_short': 0,
        'umi_diversity': defaultdict(int)
    }
    
    try:
        with gzip.open(r1_in, 'rt') as f1, gzip.open(r2_in, 'rt') as f2, \
             gzip.open(r1_out, 'wt') as o1, gzip.open(r2_out, 'wt') as o2:
            
            while True:
                # Read four lines from each file
                r1_lines = [f1.readline() for _ in range(4)]
                r2_lines = [f2.readline() for _ in range(4)]
                
                # Check for end of file
                if not r1_lines[0] or not r2_lines[0]:
                    break
                
                stats['total_reads'] += 1
                
                # Check minimum length
                if len(r1_lines[1].strip()) < min_length or \
                   len(r2_lines[1].strip()) < min_length:
                    stats['filtered_short'] += 1
                    continue
                
                # Generate deduplication key
                key = get_dedup_key(r1_lines, r2_lines, method)
                
                # Track UMI diversity
                umi = extract_umi_from_header(r1_lines[0])
                stats['umi_diversity'][umi] += 1
                
                # Check for duplicates
                if key not in seen:
                    seen.add(key)
                    stats['unique_reads'] += 1
                    
                    # Write to output
                    for line in r1_lines:
                        o1.write(line)
                    for line in r2_lines:
                        o2.write(line)
                else:
                    stats['duplicate_reads'] += 1
        
        # Calculate additional statistics
        stats['dedup_rate'] = (stats['duplicate_reads'] / stats['total_reads'] * 100) \
                              if stats['total_reads'] > 0 else 0
        stats['unique_umis'] = len(stats['umi_diversity'])
        
        # Convert defaultdict to regular dict for JSON serialization
        stats['umi_diversity'] = dict(stats['umi_diversity'])
        
        logger.info(f"  {well_id}: {stats['unique_reads']:,}/{stats['total_reads']:,} unique reads "
                   f"({100 - stats['dedup_rate']:.1f}% retained)")
        
        return stats
        
    except Exception as e:
        logger.error(f"Error processing {well_id}: {str(e)}")
        return None

def main():
    args = parse_args()
    
    # Set up logging
    logger = setup_logging(args.log)
    logger.info("Starting deduplication process")
    logger.info(f"Parameters: {vars(args)}")
    
    # Set up directories
    indir = Path(args.indir)
    outdir = Path(args.outdir)
    outdir.mkdir(parents=True, exist_ok=True)
    
    # Determine which wells to process
    if args.wells:
        wells = args.wells
    else:
        # Find all R1 files in input directory
        r1_files = list(indir.glob(f"*_R1{args.suffix}"))
        wells = [f.name.replace(f"_R1{args.suffix}", "") for f in r1_files]
        wells = sorted(set(wells))
    
    logger.info(f"Found {len(wells)} wells to process")
    
    # Process each well
    all_stats = {
        'method': args.method,
        'total_wells': len(wells),
        'wells': {},
        'summary': {
            'total_reads': 0,
            'total_unique': 0,
            'total_duplicates': 0,
            'total_filtered_short': 0
        }
    }
    
    for well_id in wells:
        well_stats = process_well(
            well_id, indir, outdir, 
            args.method, args.min_length, 
            args.suffix, logger
        )
        
        if well_stats:
            all_stats['wells'][well_id] = well_stats
            all_stats['summary']['total_reads'] += well_stats['total_reads']
            all_stats['summary']['total_unique'] += well_stats['unique_reads']
            all_stats['summary']['total_duplicates'] += well_stats['duplicate_reads']
            all_stats['summary']['total_filtered_short'] += well_stats['filtered_short']
    
    # Calculate overall statistics
    if all_stats['summary']['total_reads'] > 0:
        all_stats['summary']['overall_dedup_rate'] = \
            all_stats['summary']['total_duplicates'] / all_stats['summary']['total_reads'] * 100
    else:
        all_stats['summary']['overall_dedup_rate'] = 0
    
    # Log summary
    logger.info("=" * 60)
    logger.info("DEDUPLICATION SUMMARY")
    logger.info("=" * 60)
    logger.info(f"Method: {args.method}")
    logger.info(f"Wells processed: {len(all_stats['wells'])}/{len(wells)}")
    logger.info(f"Total reads: {all_stats['summary']['total_reads']:,}")
    logger.info(f"Unique reads: {all_stats['summary']['total_unique']:,}")
    logger.info(f"Duplicate reads: {all_stats['summary']['total_duplicates']:,}")
    logger.info(f"Overall duplication rate: {all_stats['summary']['overall_dedup_rate']:.1f}%")
    logger.info(f"Reads filtered (too short): {all_stats['summary']['total_filtered_short']:,}")
    
    # Save statistics
    if args.stats:
        # Create summary TSV for easy viewing
        summary_tsv = outdir / "dedup_summary.tsv"
        with open(summary_tsv, 'w') as f:
            f.write("well_id\ttotal_reads\tunique_reads\tduplicate_reads\tdedup_rate\n")
            for well_id, stats in all_stats['wells'].items():
                f.write(f"{well_id}\t{stats['total_reads']}\t{stats['unique_reads']}\t"
                       f"{stats['duplicate_reads']}\t{stats['dedup_rate']:.2f}\n")
        
        # Save detailed JSON statistics
        with open(args.stats, 'w') as f:
            # Remove UMI diversity for cleaner JSON (can be large)
            stats_to_save = all_stats.copy()
            for well_stats in stats_to_save['wells'].values():
                if 'umi_diversity' in well_stats:
                    del well_stats['umi_diversity']
            json.dump(stats_to_save, f, indent=2)
        
        logger.info(f"Statistics saved to {args.stats}")
        logger.info(f"Summary table saved to {summary_tsv}")
    
    logger.info("Deduplication completed successfully")

if __name__ == "__main__":
    main()
