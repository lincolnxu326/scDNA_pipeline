#!/usr/bin/env python3
"""
Compile QC metrics across all plates and generate combined HTML report.
"""

import json
import pandas as pd
from pathlib import Path
import argparse
import sys
from typing import Dict, List
import matplotlib.pyplot as plt
import seaborn as sns
from datetime import datetime

def parse_args():
    parser = argparse.ArgumentParser(
        description="Compile QC metrics across all plates"
    )
    parser.add_argument('--output-dir', required=True,
                       help='Base output directory containing plate folders')
    parser.add_argument('--plates', required=True,
                       help='Comma-separated list of plate names')
    parser.add_argument('--output-html', required=True,
                       help='Output HTML report file')
    parser.add_argument('--output-tsv', required=True,
                       help='Output TSV summary file')
    return parser.parse_args()

def load_json_safe(file_path: Path) -> Dict:
    """Safely load JSON file."""
    if file_path.exists():
        try:
            with open(file_path, 'r') as f:
                return json.load(f)
        except:
            return {}
    return {}

def load_tsv_safe(file_path: Path) -> pd.DataFrame:
    """Safely load TSV file."""
    if file_path.exists():
        try:
            return pd.read_csv(file_path, sep='\t')
        except:
            return pd.DataFrame()
    return pd.DataFrame()

def compile_plate_metrics(output_dir: Path, plates: List[str]) -> Dict:
    """Compile metrics from all plates."""
    
    all_metrics = {
        'plates': {},
        'summary': {
            'total_plates': len(plates),
            'total_reads': 0,
            'total_wells': 0,
            'average_duplication_rate': 0,
            'average_mapping_rate': 0,
            'average_adapter_rate': 0
        }
    }
    
    dup_rates = []
    map_rates = []
    adapter_rates = []
    
    for plate in plates:
        plate_dir = output_dir / plate
        if not plate_dir.exists():
            print(f"Warning: Plate directory not found: {plate_dir}")
            continue
        
        plate_metrics = {
            'demux': {},
            'dedup': {},
            'filter': {},
            'alignment': {},
            'multiqc': False
        }
        
        # Load demux stats
        demux_file = plate_dir / "demux" / "demux_stats.json"
        demux_stats = load_json_safe(demux_file)
        if demux_stats:
            plate_metrics['demux'] = demux_stats
            all_metrics['summary']['total_reads'] += demux_stats.get('total_reads', 0)
            if 'wells' in demux_stats:
                all_metrics['summary']['total_wells'] += len(demux_stats['wells']) - 1  # Exclude unassigned
        
        # Load dedup stats
        dedup_file = plate_dir / "dedup" / "dedup_stats.json"
        dedup_stats = load_json_safe(dedup_file)
        if dedup_stats:
            plate_metrics['dedup'] = dedup_stats
            if 'summary' in dedup_stats and 'overall_dedup_rate' in dedup_stats['summary']:
                dup_rates.append(dedup_stats['summary']['overall_dedup_rate'])
        
        # Load filter stats
        filter_file = plate_dir / "filtered" / "filter_stats.json"
        filter_stats = load_json_safe(filter_file)
        if filter_stats:
            plate_metrics['filter'] = filter_stats
            if 'summary' in filter_stats and 'overall_dimer_rate' in filter_stats['summary']:
                adapter_rates.append(filter_stats['summary']['overall_dimer_rate'])
        
        # Check for alignment stats
        bam_dir = plate_dir / "bam"
        if bam_dir.exists():
            flagstat_files = list(bam_dir.glob("*.flagstat.txt"))
            if flagstat_files:
                mapping_rates = []
                for ffile in flagstat_files:
                    with open(ffile, 'r') as f:
                        for line in f:
                            if "mapped (" in line:
                                rate = float(line.split('(')[1].split('%')[0])
                                mapping_rates.append(rate)
                                break
                if mapping_rates:
                    avg_map_rate = sum(mapping_rates) / len(mapping_rates)
                    plate_metrics['alignment']['average_mapping_rate'] = avg_map_rate
                    map_rates.append(avg_map_rate)
        
        # Check for MultiQC report
        multiqc_file = plate_dir / "multiqc" / "multiqc_report.html"
        plate_metrics['multiqc'] = multiqc_file.exists()
        
        all_metrics['plates'][plate] = plate_metrics
    
    # Calculate summary statistics
    if dup_rates:
        all_metrics['summary']['average_duplication_rate'] = sum(dup_rates) / len(dup_rates)
    if map_rates:
        all_metrics['summary']['average_mapping_rate'] = sum(map_rates) / len(map_rates)
    if adapter_rates:
        all_metrics['summary']['average_adapter_rate'] = sum(adapter_rates) / len(adapter_rates)
    
    return all_metrics

def generate_html_report(metrics: Dict, output_file: Path):
    """Generate HTML report with QC metrics."""
    
    html_content = f"""
<!DOCTYPE html>
<html lang="en">
<head>
    <meta charset="UTF-8">
    <meta name="viewport" content="width=device-width, initial-scale=1.0">
    <title>Single-Cell DNA Pipeline - Combined QC Report</title>
    <style>
        body {{
            font-family: 'Segoe UI', Tahoma, Geneva, Verdana, sans-serif;
            margin: 0;
            padding: 20px;
            background: linear-gradient(135deg, #667eea 0%, #764ba2 100%);
            min-height: 100vh;
        }}
        .container {{
            max-width: 1400px;
            margin: 0 auto;
            background: white;
            border-radius: 10px;
            box-shadow: 0 10px 40px rgba(0,0,0,0.1);
            padding: 30px;
        }}
        h1 {{
            color: #2c3e50;
            border-bottom: 3px solid #667eea;
            padding-bottom: 10px;
            margin-bottom: 30px;
        }}
        h2 {{
            color: #34495e;
            margin-top: 30px;
            border-bottom: 1px solid #ecf0f1;
            padding-bottom: 5px;
        }}
        .summary-grid {{
            display: grid;
            grid-template-columns: repeat(auto-fit, minmax(200px, 1fr));
            gap: 20px;
            margin: 20px 0;
        }}
        .stat-card {{
            background: linear-gradient(135deg, #667eea 0%, #764ba2 100%);
            color: white;
            padding: 20px;
            border-radius: 8px;
            text-align: center;
            box-shadow: 0 4px 6px rgba(0,0,0,0.1);
        }}
        .stat-value {{
            font-size: 2em;
            font-weight: bold;
            margin: 10px 0;
        }}
        .stat-label {{
            font-size: 0.9em;
            opacity: 0.9;
        }}
        table {{
            width: 100%;
            border-collapse: collapse;
            margin: 20px 0;
        }}
        th {{
            background: #667eea;
            color: white;
            padding: 12px;
            text-align: left;
            font-weight: 600;
        }}
        td {{
            padding: 10px 12px;
            border-bottom: 1px solid #ecf0f1;
        }}
        tr:hover {{
            background: #f8f9fa;
        }}
        .good {{ color: #27ae60; font-weight: bold; }}
        .warning {{ color: #f39c12; font-weight: bold; }}
        .bad {{ color: #e74c3c; font-weight: bold; }}
        .timestamp {{
            text-align: right;
            color: #7f8c8d;
            font-size: 0.9em;
            margin-top: 30px;
        }}
        .badge {{
            display: inline-block;
            padding: 3px 8px;
            border-radius: 4px;
            font-size: 0.85em;
            font-weight: bold;
        }}
        .badge-success {{ background: #d4edda; color: #155724; }}
        .badge-warning {{ background: #fff3cd; color: #856404; }}
        .badge-danger {{ background: #f8d7da; color: #721c24; }}
        .progress-bar {{
            width: 100%;
            height: 20px;
            background: #ecf0f1;
            border-radius: 10px;
            overflow: hidden;
            margin: 5px 0;
        }}
        .progress-fill {{
            height: 100%;
            background: linear-gradient(90deg, #27ae60, #667eea);
            transition: width 0.3s ease;
        }}
    </style>
</head>
<body>
    <div class="container">
        <h1>🧬 Single-Cell DNA Sequencing Pipeline - Combined QC Report</h1>
        
        <div class="summary-grid">
            <div class="stat-card">
                <div class="stat-label">Total Plates</div>
                <div class="stat-value">{metrics['summary']['total_plates']}</div>
            </div>
            <div class="stat-card">
                <div class="stat-label">Total Reads</div>
                <div class="stat-value">{metrics['summary']['total_reads']:,.0f}</div>
            </div>
            <div class="stat-card">
                <div class="stat-label">Total Wells</div>
                <div class="stat-value">{metrics['summary']['total_wells']}</div>
            </div>
            <div class="stat-card">
                <div class="stat-label">Avg Duplication Rate</div>
                <div class="stat-value">{metrics['summary']['average_duplication_rate']:.1f}%</div>
            </div>
            <div class="stat-card">
                <div class="stat-label">Avg Mapping Rate</div>
                <div class="stat-value">{metrics['summary']['average_mapping_rate']:.1f}%</div>
            </div>
            <div class="stat-card">
                <div class="stat-label">Avg Adapter Rate</div>
                <div class="stat-value">{metrics['summary']['average_adapter_rate']:.1f}%</div>
            </div>
        </div>
        
        <h2>📊 Plate-by-Plate Summary</h2>
        <table>
            <thead>
                <tr>
                    <th>Plate</th>
                    <th>Total Reads</th>
                    <th>Assigned Reads</th>
                    <th>Duplication Rate</th>
                    <th>Adapter Rate</th>
                    <th>Mapping Rate</th>
                    <th>MultiQC</th>
                    <th>Status</th>
                </tr>
            </thead>
            <tbody>
    """
    
    # Add plate rows
    for plate, data in metrics['plates'].items():
        # Extract metrics
        total_reads = data.get('demux', {}).get('total_reads', 0)
        assigned_reads = data.get('demux', {}).get('assigned_reads', 0)
        dup_rate = data.get('dedup', {}).get('summary', {}).get('overall_dedup_rate', 0)
        adapter_rate = data.get('filter', {}).get('summary', {}).get('overall_dimer_rate', 0)
        map_rate = data.get('alignment', {}).get('average_mapping_rate', 0)
        has_multiqc = data.get('multiqc', False)
        
        # Determine status
        issues = []
        if assigned_reads / max(total_reads, 1) < 0.7:
            issues.append("Low assignment")
        if dup_rate > 50:
            issues.append("High duplication")
        if adapter_rate > 10:
            issues.append("High adapters")
        if map_rate < 70 and map_rate > 0:
            issues.append("Low mapping")
        
        if not issues:
            status = '<span class="badge badge-success">✓ Good</span>'
        elif len(issues) == 1:
            status = f'<span class="badge badge-warning">⚠ {issues[0]}</span>'
        else:
            status = f'<span class="badge badge-danger">✗ Multiple issues</span>'
        
        # Color coding for rates
        dup_class = "good" if dup_rate < 30 else "warning" if dup_rate < 50 else "bad"
        adapter_class = "good" if adapter_rate < 5 else "warning" if adapter_rate < 10 else "bad"
        map_class = "good" if map_rate > 80 else "warning" if map_rate > 70 else "bad"
        
        multiqc_status = "✓" if has_multiqc else "✗"
        
        html_content += f"""
                <tr>
                    <td><strong>{plate}</strong></td>
                    <td>{total_reads:,}</td>
                    <td>{assigned_reads:,}</td>
                    <td class="{dup_class}">{dup_rate:.1f}%</td>
                    <td class="{adapter_class}">{adapter_rate:.1f}%</td>
                    <td class="{map_class}">{map_rate:.1f}%</td>
                    <td>{multiqc_status}</td>
                    <td>{status}</td>
                </tr>
        """
    
    html_content += f"""
            </tbody>
        </table>
        
        <h2>🎯 Quality Thresholds</h2>
        <table>
            <tr>
                <th>Metric</th>
                <th>Good</th>
                <th>Warning</th>
                <th>Bad</th>
            </tr>
            <tr>
                <td>Barcode Assignment</td>
                <td class="good">&gt; 70%</td>
                <td class="warning">50-70%</td>
                <td class="bad">&lt; 50%</td>
            </tr>
            <tr>
                <td>Duplication Rate</td>
                <td class="good">&lt; 30%</td>
                <td class="warning">30-50%</td>
                <td class="bad">&gt; 50%</td>
            </tr>
            <tr>
                <td>Adapter Contamination</td>
                <td class="good">&lt; 5%</td>
                <td class="warning">5-10%</td>
                <td class="bad">&gt; 10%</td>
            </tr>
            <tr>
                <td>Mapping Rate</td>
                <td class="good">&gt; 80%</td>
                <td class="warning">70-80%</td>
                <td class="bad">&lt; 70%</td>
            </tr>
        </table>
        
        <div class="timestamp">
            Report generated: {datetime.now().strftime('%Y-%m-%d %H:%M:%S')}
        </div>
    </div>
</body>
</html>
    """
    
    with open(output_file, 'w') as f:
        f.write(html_content)

def generate_summary_tsv(metrics: Dict, output_file: Path):
    """Generate TSV summary table."""
    
    rows = []
    for plate, data in metrics['plates'].items():
        row = {
            'plate': plate,
            'total_reads': data.get('demux', {}).get('total_reads', 0),
            'assigned_reads': data.get('demux', {}).get('assigned_reads', 0),
            'unique_reads': data.get('dedup', {}).get('summary', {}).get('total_unique', 0),
            'duplication_rate': data.get('dedup', {}).get('summary', {}).get('overall_dedup_rate', 0),
            'kept_pairs': data.get('filter', {}).get('summary', {}).get('total_kept', 0),
            'adapter_rate': data.get('filter', {}).get('summary', {}).get('overall_dimer_rate', 0),
            'mapping_rate': data.get('alignment', {}).get('average_mapping_rate', 0),
            'multiqc_generated': data.get('multiqc', False)
        }
        rows.append(row)
    
    df = pd.DataFrame(rows)
    df.to_csv(output_file, sep='\t', index=False)

def main():
    args = parse_args()
    
    output_dir = Path(args.output_dir)
    plates = args.plates.split(',')
    
    print(f"Compiling QC metrics for {len(plates)} plates...")
    
    # Compile metrics
    metrics = compile_plate_metrics(output_dir, plates)
    
    # Generate HTML report
    generate_html_report(metrics, Path(args.output_html))
    print(f"HTML report generated: {args.output_html}")
    
    # Generate TSV summary
    generate_summary_tsv(metrics, Path(args.output_tsv))
    print(f"TSV summary generated: {args.output_tsv}")
    
    # Print summary to console
    print("\n" + "="*60)
    print("COMBINED QC SUMMARY")
    print("="*60)
    print(f"Total plates: {metrics['summary']['total_plates']}")
    print(f"Total reads: {metrics['summary']['total_reads']:,}")
    print(f"Total wells: {metrics['summary']['total_wells']}")
    print(f"Average duplication rate: {metrics['summary']['average_duplication_rate']:.1f}%")
    print(f"Average mapping rate: {metrics['summary']['average_mapping_rate']:.1f}%")
    print(f"Average adapter rate: {metrics['summary']['average_adapter_rate']:.1f}%")
    print("="*60)

if __name__ == "__main__":
    main()
