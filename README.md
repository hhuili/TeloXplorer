![image](https://github.com/hhuili/TeloXplorer/blob/main/logo/logo.svg)
# TeloXplorer

[![License: MIT](https://img.shields.io/badge/License-MIT-yellow.svg)](https://opensource.org/licenses/MIT)
[![Conda Version](https://img.shields.io/conda/vn/bioconda/teloxplorer.svg?style=flat-square)](https://anaconda.org/bioconda/teloxplorer)

Teloxplorer is a pipeline for chromosome-specific telomere analysis from long-read sequencing data (Nanopore and PacBio).

## Table of Contents

- [Installation](#installation)
- [Quick Start](#quick-start)
  - Using Presets (`--preset`)
- [Telomere Methylation](#telomere-methylation)
- [Parameters](#parameters)
- [Outputs](#outputs)
- [License](#license)
- [Citation](#citation)

## Installation

```
# install with conda (Recommended)
conda install -c bioconda::teloxplorer
```

## Quick Start

Below is a standard command for running teloxplorer:
```
teloxplorer \
    --preset $species \
    --long-read-fastq $long_read_fastq \
    --ref-genome $ref_genome \
    --tel-repeats $tel_repeats \
    --mm2-preset "$mm2_preset" \
    --out-prefix $sample \
    --outdir $batch_id \
    --threads $threads \
    --plot -H 3 -W 10
```

**Using Presets (`--preset`)**

To simplify configuration, we provide optimized parameter presets for several common species.

**Recommended Presets:**

- `--preset human` (Human)
- `--preset yeast` (Yeast)
- `--preset mouse` (Mouse)
- `--preset arabidopsis` (Arabidopsis)

**For Other Species:**

If your species is not on this list, we recommend comparing its telomere repeat sequence to those of the species above. Choose the preset from the species with the most similar repeat sequence as a starting point.

For example, many mammals (e.g., cow, dog) share the same `TTAGGG` repeat as humans, making `--preset human` a good initial choice.

**Customization and Parameter Overrides**

The `--preset` flag automatically sets default values for many parameters. You can **override** any of these defaults by explicitly specifying the parameter in your command.

This is particularly useful for fine-tuning the analysis for "other species." The most common parameters to adjust are:

- `--genome_subtel_range`: Size of the subtelomeric region (bp) from the chromosome ends. Telomeres within this range are classified as `terminal`, while those outside are `internal`.
- `--min_tel_freq`: Sets the minimum frequency of telomeric repeats required to classify a read as telomeric.
- `--bloom_options`: Adjusts parameters for the BLOOM merging.

**Override Example:**

Imagine you want to use the human preset but require a stricter frequency threshold. You would run:

```
teloxplorer \
    --preset human \
    --min_tel_freq 0.4 \
    ... # other arguments
```
In this command, teloxplorer will load all settings from the human preset but will use 0.4 for min_tel_freq instead of the preset's default value.

## Telomere Methylation

```
telox-methyl \
    --preset $species \
    --tel-reads $sample.chromtel_length.tsv \
    --modBAM $bam \
    --ref-genome $ref_genome \
    --tel-repeats $tel_repeats \
    --out-prefix $sample \
    --outdir $batch_id \
    --threads $threads \
    --plot -W 8 -H 10
```

## Parameters

A detailed list of all command-line arguments: 
- `--long-read-fastq`: Input long-read FASTQ file (gzipped or plain)
- `--ref-genome`: Reference genome FASTA file
- `--outdir`: Output directory
- `--out-prefix`: Prefix for all output files
- `--preset`: Use a built-in parameter preset (human, yeast, etc.)
- `--tel-repeats`: Telomere repeat definition file (one per line: arm<TAB>repeat, regex supported)
- `--genome_subtel_range`: Size of the subtelomeric region (bp) from the chromosome ends. Telomeres within this range are classified as `terminal`, while those outside are `internal`. (Varies by preset)
- `--min_tel_freq`: Frequency threshold for telomere segment definition (Varies by preset)
- `--bloom_options`: BLOOM merging parameters (Varies by preset)	
- `--mm2-preset`: minimap2 preset (e.g., map-pb, map-ont)
- `--plot`:	Generate summary plots
- `-H, --height`: Plot height in inches
- `-W, --width`: Plot width in inches
- `--threads`: Number of threads to use

## Outputs

teloxplorer will create an output directory specified by --outdir. Key output files include:
- `$prefix.chromtel_length.tsv`: A tab-separated file containing chomosome-specific telomere length estimates.
- `$prefix.chromtel_summary.tsv`: A summary file of chomosome-specific telomere length with statistics (mean, median, etc.).
- `$prefix.chromtel.tel_length.pdf`: (--plot) A boxplot showing the telomere length accross chromosomes.
- `$prefix.neotel_length.tsv`: A tab-separated file containing neotelomere length estimates.
- `$prefix.minitel_length.tsv`: A tab-separated file containing minitelomere length estimates.
- `$prefix.chromtel.TVR.tsv`: (--find-TVR "yes")
- `$prefix.chromtel.TVR_summary.tsv`:

## Dependencies

- Minimap2 (https://github.com/lh3/minimap2)
- Python 3 (with the following auto-installed packages)
  - numpy
  - pysam
  - regex

## License
This project is licensed under the MIT License - see the LICENSE file for details.

## Citation







