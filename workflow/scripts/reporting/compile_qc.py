#!/usr/bin/env python3
"""
Compile QC metrics from all pipeline stages.
Aggregates statistics from demultiplexing, deduplication, filtering, and alignment.
"""

import json
import pandas as pd
from pathlib import Path
import argparse
import logging
import sys
from typing import Dict, List
import re

def setup_logging(log_level: str = "INFO"):
    """Set up logging configuration."""
    logging.basicConfig(
        level=getattr(logging, log_level),
        format='%(asctime)s - %(name)s - %(levelname)s - %(message)s',
        handlers=[logging.StreamHandler(sys.stdout)]
    )
    return logging.getLogger(__name__)

def parse_args():
    parser = argparse.ArgumentParser(
        description="Compile QC metrics from all pipeline stages."
    )
    parser.add_argument('--results-dir', required=True,
                       help="Pipeline results directory")
    parser.add_argument('--output-json', required=True,
                       help="Output JSON file with all QC metrics")
    parser.add_argument('--output-tsv', required=True,
                       help="Output TSV file with summary table")
    parser.add_argument('--log-level', default='INFO',
                       choices=['DEBUG', 'INFO', 'WARNING', 'ERROR'],
                       help="Logging level")
    return parser.parse_args()

def load_json_stats(file_path: Path) -> Dict:
    """Load statistics from JSON file."""
    if file_path.exists():
        with open(file_path, 'r') as f:
            return json.load(f)
    return {}

def parse_alignment_stats(stats_file: Path) -> Dict:
    """Parse samtools stats output."""
    stats = {
        'total_reads': 0,
        'mapped_reads': 0,
        'properly_paired': 0,
        'average_quality': 0
    }
    
    if not stats_file.exists():
        return stats
    
    with open(stats_file, 'r') as f:
        for line in f:
            if line.startswith('SN\traw total sequences:'):
                stats['total_reads'] = int(line.split('\t')[2])
            elif line.startswith('SN\treads mapped:'):
                stats['mapped_reads'] = int(line.split('\t')[2])
            elif line.startswith('SN\treads properly paired:'):
                stats['properly_paired'] = int(line.split('\t')[2])
            elif line.startswith('SN\taverage quality:'):
                stats['average_quality'] = float(line.split('\t')[2])
    
    if stats['total_reads'] > 0:
        stats['mapping_rate'] = stats['mapped_reads'] / stats['total_reads']
    else:
        stats['mapping_rate'] = 0
    
    return stats

def compile_demux_stats(results_dir: Path, plates: List[str], logger) -> Dict:
    """Compile demultiplexing statistics."""
    logger.info("Compiling demultiplexing statistics...")
    
    demux_stats = {}
    for plate in plates:
        stats_file = results_dir / "01_demux" / plate / "demux_stats.json"
        if stats_file.exists():
            demux_stats[plate] = load_json_stats(stats_file)
        else:
            logger.warning(f"Demux stats not found for {plate}")
    
    return demux_stats

def compile_dedup_stats(results_dir: Path, plates: List[str], logger) -> Dict:
    """Compile deduplication statistics."""
    logger.info("Compiling deduplication statistics...")
    
    dedup_stats = {}
    for plate in plates:
        stats_file = results_dir / "02_dedup" / plate / "dedup_stats.json"
        if stats_file.exists():
            dedup_stats[plate] = load_json_stats(stats_file)
        else:
            logger.warning(f"Dedup stats not found for {plate}")
    
    return dedup_stats

def compile_filter_stats(results_dir: Path, plates: List[str], logger) -> Dict:
    """Compile adapter filtering statistics."""
    logger.info("Compiling filter statistics...")
    
    filter_stats = {}
    for plate in plates:
        stats_file = results_dir / "03_filtered" / plate / "filter_stats.json"
        if stats_file.exists():
            filter_stats[plate] = load_json_stats(stats_file)
        else:
            logger.warning(f"Filter stats not found for {plate}")
    
    return filter_stats

def compile_alignment_stats(results_dir: Path, plates: List[str], logger) -> Dict:
    """Compile alignment statistics."""
    logger.info("Compiling alignment statistics...")
    
    align_stats = {}
    for plate in plates:
        align_dir = results_dir / "04_aligned" / plate
        if not align_dir.exists():
            logger.warning(f"Alignment directory not found for {plate}")
            continue
        
        align_stats[plate] = {}
        for stats_file in align_dir.glob("*.stats.txt"):
            well = stats_file.stem.replace('.stats', '')
            align_stats[plate][well] = parse_alignment_stats(stats_file)
    
    return align_stats

def compile_aneufinder_stats(results_dir: Path, plates: List[str], logger) -> Dict:
    """Check AneuFinder completion status."""
    logger.info("Checking AneuFinder completion...")
    
    aneufinder_stats = {}
    for plate in plates:
        flag_file = results_dir / "05_aneufinder" / plate / "aneufinder_complete.flag"
        aneufinder_stats[plate] = {
            'completed': flag_file.exists(),
            'output_dir': str(results_dir / "05_aneufinder" / plate / "output")
        }
    
    return aneufinder_stats

def create_summary_table(all_stats: Dict, output_file: Path, logger):
    """Create summary TSV table."""
    logger.info("Creating summary table...")
    
    rows = []
    
    for plate in all_stats.get('plates', []):
        # Get stats for this plate
        demux = all_stats['demux'].get(plate, {})
        dedup = all_stats['dedup'].get(plate, {})
        filter_stats = all_stats['filter'].get(plate, {})
        align = all_stats['alignment'].get(plate, {})
        aneufinder = all_stats['aneufinder'].get(plate, {})
        
        # Calculate summary metrics
        row = {
            'plate': plate,
            'total_reads': demux.get('total_reads', 0),
            'assigned_reads': demux.get('assigned_reads', 0),
            'unassigned_reads': demux.get('unassigned_reads', 0)
        }
        
        # Add dedup stats
        if dedup and 'summary' in dedup:
            row['unique_reads'] = dedup['summary'].get('total_unique', 0)
            row['duplicate_reads'] = dedup['summary'].get('total_duplicates', 0)
            row['dedup_rate'] = dedup['summary'].get('overall_dedup_rate', 0)
        
        # Add filter stats
        if filter_stats and 'summary' in filter_stats:
            row['kept_pairs'] = filter_stats['summary'].get('total_kept', 0)
            row['dimer_pairs'] = filter_stats['summary'].get('total_removed', 0)
            row['dimer_rate'] = filter_stats['summary'].get('overall_dimer_rate', 0)
        
        # Add alignment stats (average across wells)
        if plate in align and align[plate]:
            total_mapped = sum(w.get('mapped_reads', 0) for w in align[plate].values())
            total_reads = sum(w.get('total_reads', 0) for w in align[plate].values())
            if total_reads > 0:
                row['avg_mapping_rate'] = (total_mapped / total_reads) * 100
            else:
                row['avg_mapping_rate'] = 0
            row['wells_aligned'] = len(align[plate])
        
        # Add AneuFinder status
        row['aneufinder_complete'] = aneufinder.get('completed', False)
        
        rows.append(row)
    
    # Create DataFrame and save
    df = pd.DataFrame(rows)
    df.to_csv(output_file, sep='\t', index=False)
    logger.info(f"Summary table saved to {output_file}")
    
    return df

def calculate_overall_metrics(all_stats: Dict) -> Dict:
    """Calculate overall pipeline metrics."""
    metrics = {
        'total_plates': len(all_stats.get('plates', [])),
        'total_reads_processed': 0,
        'total_unique_reads': 0,
        'overall_dedup_rate': 0,
        'overall_dimer_rate': 0,
        'overall_mapping_rate': 0,
        'plates_completed': 0
    }
    
    total_reads = 0
    total_duplicates = 0
    total_dimers = 0
    total_pairs = 0
    total_mapped = 0
    total_alignment_reads = 0
    
    for plate in all_stats.get('plates', []):
        # Demux stats
        if plate in all_stats.get('demux', {}):
            metrics['total_reads_processed'] += all_stats['demux'][plate].get('total_reads', 0)
        
        # Dedup stats
        if plate in all_stats.get('dedup', {}) and 'summary' in all_stats['dedup'][plate]:
            summary = all_stats['dedup'][plate]['summary']
            total_reads += summary.get('total_reads', 0)
            total_duplicates += summary.get('total_duplicates', 0)
            metrics['total_unique_reads'] += summary.get('total_unique', 0)
        
        # Filter stats
        if plate in all_stats.get('filter', {}) and 'summary' in all_stats['filter'][plate]:
            summary = all_stats['filter'][plate]['summary']
            total_pairs += summary.get('total_pairs', 0)
            total_dimers += summary.get('total_removed', 0)
        
        # Alignment stats
        if plate in all_stats.get('alignment', {}):
            for well_stats in all_stats['alignment'][plate].values():
                total_mapped += well_stats.get('mapped_reads', 0)
                total_alignment_reads += well_stats.get('total_reads', 0)
        
        # AneuFinder completion
        if plate in all_stats.get('aneufinder', {}) and all_stats['aneufinder'][plate].get('completed'):
            metrics['plates_completed'] += 1
    
    # Calculate rates
    if total_reads > 0:
        metrics['overall_dedup_rate'] = (total_duplicates / total_reads) * 100
    
    if total_pairs > 0:
        metrics['overall_dimer_rate'] = (total_dimers / total_pairs) * 100
    
    if total_alignment_reads > 0:
        metrics['overall_mapping_rate'] = (total_mapped / total_alignment_reads) * 100
    
    return metrics

def main():
    args = parse_args()
    logger = setup_logging(args.log_level)
    
    logger.info("Starting QC compilation...")
    
    results_dir = Path(args.results_dir)
    if not results_dir.exists():
        logger.error(f"Results directory not found: {results_dir}")
        sys.exit(1)
    
    # Find all plates
    plates = []
    demux_dir = results_dir / "01_demux"
    if demux_dir.exists():
        plates = [d.name for d in demux_dir.iterdir() if d.is_dir()]
    
    logger.info(f"Found {len(plates)} plates: {', '.join(plates)}")
    
    # Compile statistics from each stage
    all_stats = {
        'plates': plates,
        'demux': compile_demux_stats(results_dir, plates, logger),
        'dedup': compile_dedup_stats(results_dir, plates, logger),
        'filter': compile_filter_stats(results_dir, plates, logger),
        'alignment': compile_alignment_stats(results_dir, plates, logger),
        'aneufinder': compile_aneufinder_stats(results_dir, plates, logger)
    }
    
    # Calculate overall metrics
    all_stats['overall_metrics'] = calculate_overall_metrics(all_stats)
    
    # Save JSON output
    output_json = Path(args.output_json)
    output_json.parent.mkdir(parents=True, exist_ok=True)
    with open(output_json, 'w') as f:
        json.dump(all_stats, f, indent=2)
    logger.info(f"QC metrics saved to {output_json}")
    
    # Create summary table
    output_tsv = Path(args.output_tsv)
    summary_df = create_summary_table(all_stats, output_tsv, logger)
    
    # Log summary
    logger.info("=" * 60)
    logger.info("PIPELINE SUMMARY")
    logger.info("=" * 60)
    for key, value in all_stats['overall_metrics'].items():
        if isinstance(value, float):
            logger.info(f"{key}: {value:.2f}")
        else:
            logger.info(f"{key}: {value}")
    
    logger.info("QC compilation completed successfully!")

if __name__ == "__main__":
    main()
