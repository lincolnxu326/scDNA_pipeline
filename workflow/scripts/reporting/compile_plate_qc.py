#!/usr/bin/env python3
"""
Compile QC metrics for a single plate.
Generates summary statistics from all processing stages.
"""

import json
import argparse
import sys
from pathlib import Path
from typing import Dict
import pandas as pd

def parse_args():
    parser = argparse.ArgumentParser(
        description="Compile QC metrics for a single plate"
    )
    parser.add_argument('--plate', required=True,
                       help='Plate identifier')
    parser.add_argument('--output-txt', required=True,
                       help='Output text summary file')
    parser.add_argument('--output-json', required=True,
                       help='Output JSON metrics file')
    return parser.parse_args()

def load_json_safe(file_path: str) -> Dict:
    """Safely load JSON file."""
    path = Path(file_path)
    if path.exists():
        with open(path, 'r') as f:
            return json.load(f)
    return {}

def parse_samtools_stats(stats_file: Path) -> Dict:
    """Parse samtools stats output."""
    stats = {
        'total_reads': 0,
        'mapped_reads': 0,
        'properly_paired': 0,
        'mapping_rate': 0.0
    }
    
    if stats_file.exists():
        with open(stats_file, 'r') as f:
            for line in f:
                if line.startswith('SN\traw total sequences:'):
                    stats['total_reads'] = int(line.split('\t')[2])
                elif line.startswith('SN\treads mapped:'):
                    stats['mapped_reads'] = int(line.split('\t')[2])
                elif line.startswith('SN\treads properly paired:'):
                    stats['properly_paired'] = int(line.split('\t')[2])
        
        if stats['total_reads'] > 0:
            stats['mapping_rate'] = (stats['mapped_reads'] / stats['total_reads']) * 100
    
    return stats

def compile_plate_metrics(plate: str) -> Dict:
    """Compile all metrics for a plate."""
    metrics = {
        'plate': plate,
        'demux': {},
        'dedup': {},
        'filter': {},
        'alignment': {},
        'aneufinder': {}
    }
    
    # Load demux stats
    demux_file = Path(f"{plate}/demux/demux_stats.json")
    if demux_file.exists():
        metrics['demux'] = load_json_safe(demux_file)
    
    # Load dedup stats
    dedup_file = Path(f"{plate}/dedup/dedup_stats.json")
    if dedup_file.exists():
        metrics['dedup'] = load_json_safe(dedup_file)
    
    # Load filter stats
    filter_file = Path(f"{plate}/filtered/filter_stats.json")
    if filter_file.exists():
        metrics['filter'] = load_json_safe(filter_file)
    
    # Load alignment stats
    bam_dir = Path(f"{plate}/bam")
    if bam_dir.exists():
        alignment_stats = {}
        for stats_file in bam_dir.glob("*.stats.txt"):
            well = stats_file.stem.replace('.stats', '')
            alignment_stats[well] = parse_samtools_stats(stats_file)
        metrics['alignment'] = alignment_stats
    
    # Check AneuFinder completion
    aneufinder_flag = Path(f"{plate}/aneufinder/complete.flag")
    metrics['aneufinder']['completed'] = aneufinder_flag.exists()
    
    # Check for AneuFinder output
    aneufinder_dir = Path(f"{plate}/aneufinder/output")
    if aneufinder_dir.exists():
        model_dir = aneufinder_dir / "MODELS"
        if model_dir.exists():
            model_files = list(model_dir.glob("**/*.RData"))
            metrics['aneufinder']['num_models'] = len(model_files)
        
        # Check for plots
        plot_dir = aneufinder_dir / "PLOTS"
        if plot_dir.exists():
            metrics['aneufinder']['plots_generated'] = True
            plot_files = list(plot_dir.glob("*.pdf"))
            metrics['aneufinder']['num_plots'] = len(plot_files)
    
    return metrics

def generate_summary_text(metrics: Dict) -> str:
    """Generate human-readable summary text."""
    lines = []
    lines.append("=" * 70)
    lines.append(f"SINGLE-CELL SEQUENCING PIPELINE QC SUMMARY")
    lines.append(f"Plate: {metrics['plate']}")
    lines.append("=" * 70)
    lines.append("")
    
    # Demux summary
    if metrics['demux']:
        lines.append("DEMULTIPLEXING")
        lines.append("-" * 30)
        lines.append(f"Total reads: {metrics['demux'].get('total_reads', 0):,}")
        lines.append(f"Assigned reads: {metrics['demux'].get('assigned_reads', 0):,}")
        lines.append(f"Unassigned reads: {metrics['demux'].get('unassigned_reads', 0):,}")
        if metrics['demux'].get('total_reads', 0) > 0:
            assign_rate = (metrics['demux'].get('assigned_reads', 0) / 
                          metrics['demux'].get('total_reads', 1)) * 100
            lines.append(f"Assignment rate: {assign_rate:.1f}%")
        lines.append("")
    
    # Dedup summary
    if metrics['dedup'] and 'summary' in metrics['dedup']:
        lines.append("DEDUPLICATION")
        lines.append("-" * 30)
        summary = metrics['dedup']['summary']
        lines.append(f"Method: {metrics['dedup'].get('method', 'unknown')}")
        lines.append(f"Total reads: {summary.get('total_reads', 0):,}")
        lines.append(f"Unique reads: {summary.get('total_unique', 0):,}")
        lines.append(f"Duplicate reads: {summary.get('total_duplicates', 0):,}")
        lines.append(f"Duplication rate: {summary.get('overall_dedup_rate', 0):.1f}%")
        lines.append("")
    
    # Filter summary
    if metrics['filter'] and 'summary' in metrics['filter']:
        lines.append("ADAPTER FILTERING")
        lines.append("-" * 30)
        summary = metrics['filter']['summary']
        lines.append(f"Total pairs: {summary.get('total_pairs', 0):,}")
        lines.append(f"Kept pairs: {summary.get('total_kept', 0):,}")
        lines.append(f"Removed (dimers): {summary.get('total_removed', 0):,}")
        lines.append(f"Dimer rate: {summary.get('overall_dimer_rate', 0):.1f}%")
        lines.append("")
    
    # Alignment summary
    if metrics['alignment']:
        lines.append("ALIGNMENT")
        lines.append("-" * 30)
        total_reads = sum(w.get('total_reads', 0) for w in metrics['alignment'].values())
        mapped_reads = sum(w.get('mapped_reads', 0) for w in metrics['alignment'].values())
        lines.append(f"Wells processed: {len(metrics['alignment'])}")
        lines.append(f"Total reads: {total_reads:,}")
        lines.append(f"Mapped reads: {mapped_reads:,}")
        if total_reads > 0:
            lines.append(f"Overall mapping rate: {(mapped_reads/total_reads)*100:.1f}%")
        
        # Per-well stats
        mapping_rates = [w.get('mapping_rate', 0) for w in metrics['alignment'].values()]
        if mapping_rates:
            lines.append(f"Mean mapping rate: {sum(mapping_rates)/len(mapping_rates):.1f}%")
            lines.append(f"Min mapping rate: {min(mapping_rates):.1f}%")
            lines.append(f"Max mapping rate: {max(mapping_rates):.1f}%")
        lines.append("")
    
    # AneuFinder summary
    lines.append("ANEUFINDER ANALYSIS")
    lines.append("-" * 30)
    if metrics['aneufinder'].get('completed'):
        lines.append("Status: COMPLETED")
        if 'num_models' in metrics['aneufinder']:
            lines.append(f"Models generated: {metrics['aneufinder']['num_models']}")
        if metrics['aneufinder'].get('plots_generated'):
            lines.append(f"Plots generated: YES ({metrics['aneufinder'].get('num_plots', 0)} files)")
    else:
        lines.append("Status: NOT COMPLETED")
    lines.append("")
    
    # Summary statistics
    lines.append("=" * 70)
    lines.append("PIPELINE SUMMARY")
    lines.append("=" * 70)
    
    # Calculate yield
    if metrics['demux'] and metrics['filter']:
        initial_reads = metrics['demux'].get('total_reads', 0)
        final_pairs = metrics['filter'].get('summary', {}).get('total_kept', 0)
        if initial_reads > 0:
            yield_pct = (final_pairs * 2 / initial_reads) * 100  # *2 because pairs
            lines.append(f"Overall yield: {yield_pct:.1f}%")
    
    # Check for issues
    issues = []
    if metrics['demux'].get('assigned_reads', 0) / max(metrics['demux'].get('total_reads', 1), 1) < 0.7:
        issues.append("Low barcode assignment rate (<70%)")
    if metrics['dedup'].get('summary', {}).get('overall_dedup_rate', 0) > 50:
        issues.append("High duplication rate (>50%)")
    if metrics['filter'].get('summary', {}).get('overall_dimer_rate', 0) > 10:
        issues.append("High dimer contamination (>10%)")
    
    if issues:
        lines.append("")
        lines.append("POTENTIAL ISSUES:")
        for issue in issues:
            lines.append(f"  ⚠️  {issue}")
    else:
        lines.append("")
        lines.append("✅ All QC metrics within expected ranges")
    
    lines.append("")
    lines.append("=" * 70)
    
    return "\n".join(lines)

def main():
    args = parse_args()
    
    print(f"Compiling QC metrics for plate: {args.plate}")
    
    # Compile metrics
    metrics = compile_plate_metrics(args.plate)
    
    # Generate summary text
    summary_text = generate_summary_text(metrics)
    
    # Write text summary
    with open(args.output_txt, 'w') as f:
        f.write(summary_text)
    print(f"Text summary written to: {args.output_txt}")
    
    # Write JSON metrics
    with open(args.output_json, 'w') as f:
        json.dump(metrics, f, indent=2)
    print(f"JSON metrics written to: {args.output_json}")
    
    # Print summary to console
    print("\n" + summary_text)

if __name__ == "__main__":
    main()
