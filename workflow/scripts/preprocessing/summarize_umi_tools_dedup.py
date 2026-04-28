#!/usr/bin/env python3
"""
Summarize umi-tools BAM deduplication output for reporting.
"""

import argparse
import json
import logging
import subprocess
import sys
from pathlib import Path


def setup_logging(log_file: str = None):
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
        description="Summarize raw and umi-tools deduplicated BAM read counts."
    )
    parser.add_argument("--raw-dir", required=True,
                        help="Directory containing raw aligned BAMs")
    parser.add_argument("--dedup-dir", required=True,
                        help="Directory containing deduplicated BAMs")
    parser.add_argument("--outdir", required=True,
                        help="Output directory for dedup summary files")
    parser.add_argument("--method", default="directional",
                        help="umi-tools dedup method")
    parser.add_argument("--stats", required=True,
                        help="Output JSON statistics file")
    parser.add_argument("--log",
                        help="Log file path")
    parser.add_argument("--wells", nargs="+", required=True,
                        help="Well IDs to summarize")

    return parser.parse_args()


def count_primary_mapped_reads(bam_file: Path) -> int:
    command = [
        "samtools", "view",
        "-c",
        "-F", "2308",
        str(bam_file)
    ]
    result = subprocess.run(command, check=True, capture_output=True, text=True)
    return int(result.stdout.strip())


def main():
    args = parse_args()
    logger = setup_logging(args.log)

    raw_dir = Path(args.raw_dir)
    dedup_dir = Path(args.dedup_dir)
    outdir = Path(args.outdir)
    outdir.mkdir(parents=True, exist_ok=True)

    all_stats = {
        "method": args.method,
        "total_wells": len(args.wells),
        "wells": {},
        "summary": {
            "total_reads": 0,
            "total_unique": 0,
            "total_duplicates": 0,
            "total_filtered_short": 0,
            "overall_dedup_rate": 0
        }
    }

    logger.info("Summarizing umi-tools deduplication counts")

    for well_id in args.wells:
        raw_bam = raw_dir / f"{well_id}.bam"
        dedup_bam = dedup_dir / f"{well_id}.bam"

        if not raw_bam.exists() or not dedup_bam.exists():
            raise FileNotFoundError(
                f"Missing raw or deduplicated BAM for {well_id}: {raw_bam}, {dedup_bam}"
            )

        total_reads = count_primary_mapped_reads(raw_bam)
        unique_reads = count_primary_mapped_reads(dedup_bam)
        duplicate_reads = max(total_reads - unique_reads, 0)
        dedup_rate = (duplicate_reads / total_reads * 100) if total_reads else 0

        all_stats["wells"][well_id] = {
            "well_id": well_id,
            "total_reads": total_reads,
            "unique_reads": unique_reads,
            "duplicate_reads": duplicate_reads,
            "dedup_rate": dedup_rate
        }

        all_stats["summary"]["total_reads"] += total_reads
        all_stats["summary"]["total_unique"] += unique_reads
        all_stats["summary"]["total_duplicates"] += duplicate_reads

        logger.info(
            f"{well_id}: {unique_reads:,}/{total_reads:,} primary mapped reads retained "
            f"({dedup_rate:.1f}% duplicates)"
        )

    if all_stats["summary"]["total_reads"] > 0:
        all_stats["summary"]["overall_dedup_rate"] = (
            all_stats["summary"]["total_duplicates"] /
            all_stats["summary"]["total_reads"] * 100
        )

    summary_tsv = outdir / "dedup_summary.tsv"
    with open(summary_tsv, "w") as handle:
        handle.write("well_id\ttotal_reads\tunique_reads\tduplicate_reads\tdedup_rate\n")
        for well_id, stats in all_stats["wells"].items():
            handle.write(
                f"{well_id}\t{stats['total_reads']}\t{stats['unique_reads']}\t"
                f"{stats['duplicate_reads']}\t{stats['dedup_rate']:.2f}\n"
            )

    with open(args.stats, "w") as handle:
        json.dump(all_stats, handle, indent=2)

    logger.info(f"Dedup summary saved to {summary_tsv}")
    logger.info(f"Dedup stats saved to {args.stats}")


if __name__ == "__main__":
    main()
