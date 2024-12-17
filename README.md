# ParseroLogan

ParseroLogan is a tool that downloads, annotates, and performs HMM-based searches on SRA contigs coming from the [Logan database](https://github.com/IndexThePlanet/Logan) . It uses `pyrodigal-gv` for ORF finding, `pyhmmer` for HMM searching, and `textual` for a UI-based progress table. It can optionally log progress to a file for use on HPC systems where interactive UIs may not be practical.

## Features

- Download contigs from S3.
- Annotate contigs with ORFs and generate GFF and protein files.
- Run HMM searches to find specific protein hits.
- Display progress in a textual UI or continuously update a log file.

## Installation

You can install ParseroLogan using [pixi](https://github.com/jdevlieghere/pixi) or standard Python packaging tools:

```bash
# Using pixi
pixi install .
    ```

Or, if you prefer a standard Python install (e.g. in a virtual environment):

```bash
pip install .
```

## Usage

```bash
parserologan INPUT_FILE HMM_FILE [OPTIONS]
```

## Options

```bash
--threads: Number of parallel downloads (default: number of CPU cores)
--min_length: Minimum contig length (default: 1000)
--evalue: E-value threshold for HMM hits (default: 1e-5)
--output_dir: Directory to store results (default: parsero_logan)
--fix_circles: Fix circular concatemers for sequences with DTRs (default: False)
--log_file: Specify a path to a file to log progress instead of the textual UI (default: None)
```

## Example

```bash
parserologan sra_ids.txt profile.hmm --threads 4 --log_file progress.log
```