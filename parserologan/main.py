import os
import sys
from pathlib import Path

import click

@click.command()
@click.argument('input_file', type=click.Path(exists=True))
@click.argument('hmm_file', type=click.Path(exists=True))
@click.option('--threads', type=int, default=None, help='Number of parallel downloads')
@click.option('--min_length', type=int, default=1000, help='Minimum contig length to consider')
@click.option('--evalue', type=float, default=1e-5, help='E-value threshold for HMM hits')
@click.option('--output_dir', type=click.Path(), default='parsero_logan', help='Output directory for results')
@click.option('--fix_circles', is_flag=True, default=False, help='Fix circular concatemers for sequences with DTRs')
@click.option('--log_file', type=click.Path(), default=None, help='Log file to write updates to instead of using the textual interface')
def main(input_file, hmm_file, threads, min_length, evalue, output_dir, fix_circles, log_file):
    if threads is None:
        threads = max(1, (os.cpu_count() or 1))

    with open(input_file) as f:
        sra_ids = [line.strip() for line in f if line.strip()]

    output_dir = Path(output_dir)
    output_dir.mkdir(parents=True, exist_ok=True)

    if log_file is None:
        from .run_textual import run_with_textual_interface
        run_with_textual_interface(sra_ids, threads, Path(hmm_file), min_length, evalue, output_dir, fix_circles)
    else:
        from .run_log_file import run_with_log_file
        run_with_log_file(sra_ids, threads, Path(hmm_file), min_length, evalue, output_dir, fix_circles, log_file)


if __name__ == '__main__':
    main()
