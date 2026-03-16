#!/usr/bin/env python3
"""
Filter out adapter dimer reads from deduplicated FASTQ files.
Aligned for single-cell sequencing pipeline.
"""

import os
import sys
import gzip
import argparse
import logging
from pathlib import Path
from typing import Tuple, Dict
import json
import io

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
        description="Filter adapter dimer reads from deduplicated FASTQs."
    )
    parser.add_argument('--indir', required=True,
                       help="Input directory with deduplicated FASTQs")
    parser.add_argument('--outdir', required=True,
                       help="Output directory for filtered FASTQs")
    parser.add_argument('--adapter', default='CAGTCAGCGT',
                       help="Adapter sequence to filter (default: CAGTCAGCGT)")
    parser.add_argument('--case-insensitive', action='store_true',
                       help="Case-insensitive adapter matching")
    parser.add_argument('--both-reads', action='store_true', default=True,
                       help="Require adapter in both R1 and R2 (default: True)")
    parser.add_argument('--suffix', default='.dedup.fastq.gz',
                       help="Input file suffix (default: .dedup.fastq.gz)")
    parser.add_argument('--out-suffix', default='.filtered.fastq.gz',
                       help="Output file suffix (default: .filtered.fastq.gz)")
    parser.add_argument('--log', help="Log file path")
    parser.add_argument('--stats', help="Output statistics JSON file")
    parser.add_argument('--wells', nargs='*',
                       help="Specific wells to process (default: all)")
    
    return parser.parse_args()

def open_gz_text(path: str, mode: str = 'rt'):
    """Open a gzipped file with optimized buffering."""
    if 'w' in mode:
        return gzip.open(path, mode)
    # Use TextIOWrapper for better read performance
    return io.TextIOWrapper(
        gzip.open(path, mode.replace('t', 'b')),
        encoding='utf-8',
        newline=''
    )

def iter_fastq_pairs(r1_handle, r2_handle):
    """Yield paired FASTQ records."""
    while True:
        r1_header = r1_handle.readline()
        if not r1_header:
            break
        r1_seq = r1_handle.readline()
        r1_plus = r1_handle.readline()
        r1_qual = r1_handle.readline()
        
        r2_header = r2_handle.readline()
        if not r2_header:
            break
        r2_seq = r2_handle.readline()
        r2_plus = r2_handle.readline()
        r2_qual = r2_handle.readline()
        
        yield (r1_header, r1_seq, r1_plus, r1_qual,
               r2_header, r2_seq, r2_plus, r2_qual)

def is_adapter_dimer(r1_seq: str, r2_seq: str, adapter: str,
                    case_insensitive: bool, both_reads: bool) -> bool:
    """
    Check if reads contain adapter dimers.
    
    Args:
        r1_seq: R1 sequence
        r2_seq: R2 sequence
        adapter: Adapter sequence to check
        case_insensitive: Whether to ignore case
        both_reads: Whether to require adapter in both reads
    
    Returns:
        True if reads are adapter dimers
    """
    if case_insensitive:
        adapter_upper = adapter.upper()
        r1_start = r1_seq[:len(adapter)].upper()
        r2_start = r2_seq[:len(adapter)].upper()
        
        if both_reads:
            return r1_start == adapter_upper and r2_start == adapter_upper
        else:
            return r1_start == adapter_upper or r2_start == adapter_upper
    else:
        if both_reads:
            return r1_seq.startswith(adapter) and r2_seq.startswith(adapter)
        else:
            return r1_seq.startswith(adapter) or r2_seq.startswith(adapter)

def process_well(well_id: str, indir: Path, outdir: Path,
                adapter: str, case_insensitive: bool,
                both_reads: bool, suffix: str, 
                out_suffix: str, logger) -> Dict:
    """Process a single well for adapter filtering."""
    
    r1_in = indir / f"{well_id}_R1{suffix}"
    r2_in = indir / f"{well_id}_R2{suffix}"
    r1_out = outdir / f"{well_id}_R1{out_suffix}"
    r2_out = outdir / f"{well_id}_R2{out_suffix}"
    
    if not r1_in.exists() or not r2_in.exists():
        logger.warning(f"Input files not found for {well_id}, skipping")
        return None
    
    logger.info(f"Processing {well_id}...")
    
    stats = {
        'well_id': well_id,
        'total_pairs': 0,
        'kept_pairs': 0,
        'removed_dimers': 0,
        'dimer_rate': 0.0
    }
    
    try:
        with open_gz_text(r1_in, 'rt') as f1, \
             open_gz_text(r2_in, 'rt') as f2, \
             open_gz_text(r1_out, 'wt') as o1, \
             open_gz_text(r2_out, 'wt') as o2:
            
            for record in iter_fastq_pairs(f1, f2):
                r1_h, r1_s, r1_p, r1_q, r2_h, r2_s, r2_p, r2_q = record
                stats['total_pairs'] += 1
                
                if is_adapter_dimer(r1_s, r2_s, adapter, 
                                   case_insensitive, both_reads):
                    stats['removed_dimers'] += 1
                    continue
                
                # Keep the pair
                stats['kept_pairs'] += 1
                o1.write(r1_h)
                o1.write(r1_s)
                o1.write(r1_p)
                o1.write(r1_q)
                
                o2.write(r2_h)
                o2.write(r2_s)
                o2.write(r2_p)
                o2.write(r2_q)
        
        # Calculate dimer rate
        if stats['total_pairs'] > 0:
            stats['dimer_rate'] = (stats['removed_dimers'] / stats['total_pairs']) * 100
        
        logger.info(f"  {well_id}: Kept {stats['kept_pairs']:,}/{stats['total_pairs']:,} pairs "
                   f"(removed {stats['removed_dimers']:,} dimers, {stats['dimer_rate']:.1f}%)")
        
        return stats
        
    except Exception as e:
        logger.error(f"Error processing {well_id}: {str(e)}")
        return None

def main():
    args = parse_args()
    
    # Set up logging
    logger = setup_logging(args.log)
    logger.info("Starting adapter dimer filtering")
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
    logger.info(f"Adapter sequence: {args.adapter}")
    logger.info(f"Matching mode: {'case-insensitive' if args.case_insensitive else 'exact'}")
    logger.info(f"Filter mode: {'both reads' if args.both_reads else 'either read'}")
    
    # Process each well
    all_stats = {
        'adapter_sequence': args.adapter,
        'case_insensitive': args.case_insensitive,
        'both_reads_required': args.both_reads,
        'total_wells': len(wells),
        'wells': {},
        'summary': {
            'total_pairs': 0,
            'total_kept': 0,
            'total_removed': 0,
            'overall_dimer_rate': 0.0
        }
    }
    
    for well_id in wells:
        well_stats = process_well(
            well_id, indir, outdir,
            args.adapter, args.case_insensitive,
            args.both_reads, args.suffix,
            args.out_suffix, logger
        )
        
        if well_stats:
            all_stats['wells'][well_id] = well_stats
            all_stats['summary']['total_pairs'] += well_stats['total_pairs']
            all_stats['summary']['total_kept'] += well_stats['kept_pairs']
            all_stats['summary']['total_removed'] += well_stats['removed_dimers']
    
    # Calculate overall statistics
    if all_stats['summary']['total_pairs'] > 0:
        all_stats['summary']['overall_dimer_rate'] = \
            (all_stats['summary']['total_removed'] / all_stats['summary']['total_pairs']) * 100
    
    # Log summary
    logger.info("=" * 60)
    logger.info("ADAPTER FILTERING SUMMARY")
    logger.info("=" * 60)
    logger.info(f"Wells processed: {len(all_stats['wells'])}/{len(wells)}")
    logger.info(f"Total read pairs: {all_stats['summary']['total_pairs']:,}")
    logger.info(f"Kept pairs: {all_stats['summary']['total_kept']:,}")
    logger.info(f"Removed dimers: {all_stats['summary']['total_removed']:,}")
    logger.info(f"Overall dimer rate: {all_stats['summary']['overall_dimer_rate']:.1f}%")
    
    # Save statistics
    if args.stats:
        # Create summary TSV
        summary_tsv = outdir / "adapter_filter_summary.tsv"
        with open(summary_tsv, 'w') as f:
            f.write("well_id\ttotal_pairs\tkept_pairs\tremoved_dimers\tdimer_rate\n")
            for well_id, stats in all_stats['wells'].items():
                f.write(f"{well_id}\t{stats['total_pairs']}\t{stats['kept_pairs']}\t"
                       f"{stats['removed_dimers']}\t{stats['dimer_rate']:.2f}\n")
        
        # Save detailed JSON statistics
        with open(args.stats, 'w') as f:
            json.dump(all_stats, f, indent=2)
        
        logger.info(f"Statistics saved to {args.stats}")
        logger.info(f"Summary table saved to {summary_tsv}")
    
    logger.info("Adapter filtering completed successfully")

if __name__ == "__main__":
    main()
