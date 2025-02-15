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
git clone https://github.com/pentamorfico/ParseroLogan.git
cd ParseroLogan
pixi install 
pixi run pip install 
```

Or, if you prefer a standard Python install (e.g. in a virtual environment). You will need to have these dependencies installed:

  * click
  * biopython
  * pyhmmer
  * pyrodigal-gv
  * textual
  * rich-click
  * awscli

```bash
pip install .
```

## Usage

```bash
parserologan INPUT_FILE HMM_FILE [OPTIONS]
```

## Options

```bash
Options:

  --threads INTEGER     Number of parallel downloads
  --min_length INTEGER  Minimum contig length to consider
  --evalue FLOAT        E-value threshold for HMM hits
  --output_dir PATH     Output directory for results
  --save_all            Save all results (fna, faa, gff) even if no HMM hits
                        were found
```

## Example

With pixi:

```bash
pixi run parserologan sra_ids.txt profile.hmm --log_file test.txt
```

Without pixi:

```bash
parserologan sra_ids.txt profile.hmm --log_file progress.log
```