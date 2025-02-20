# main.py
import click
import multiprocessing
from pathlib import Path
import os
from .run_log_file import run_with_log_file
from .run_textual import run_with_textual_interface

@click.group()
def cli():
    """Logan Parsero: A tool for parsing and analyzing SRA data."""
    pass

@cli.command()
@click.argument('input_file', type=click.Path(exists=True))
@click.argument('hmm_file', type=click.Path(exists=True))
@click.option('--threads', type=int, default=multiprocessing.cpu_count(), help='Number of worker threads to use')
@click.option('--min-length', type=int, default=1000, help='Minimum length of sequences to process')
@click.option('--evalue', type=float, default=1e-10, help='E-value threshold for HMMER search')
@click.option('--output-dir', type=click.Path(), default=Path('results'), help='Directory to output results')
@click.option('--log-file', type=click.Path(), default=None, help='Log file to write updates to instead of using the textual interface')
@click.option('--save-all', is_flag=True, default=False, help='Save all results (fna, faa, gff) even if no HMM hits were found')

def main(input_file, hmm_file, threads, min_length, evalue, output_dir, log_file, save_all):
    sra_ids = [line.strip() for line in open(input_file).readlines() if line.strip()]
    os.makedirs(output_dir,exist_ok=True)
    
    if log_file:
        run_with_log_file(
            sra_ids,
            threads,
            hmm_file,
            min_length,
            evalue,
            output_dir,
            log_file,
            save_all
        )
    else:
        run_with_textual_interface(
            sra_ids,
            threads,
            hmm_file,
            min_length,
            evalue,
            output_dir,
            save_all
        )

if __name__ == "__main__":
    cli()
