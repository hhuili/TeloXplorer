![image](https://github.com/hhuili/TeloXplorer/blob/main/logo/logo.svg)
## TeloXplorer

[![Conda Version](https://img.shields.io/conda/vn/huihui_li/teloxplorer.svg?style=flat-square)](https://anaconda.org/huihui_li/teloxplorer)

Teloxplorer is a tool for chromosome-specific telomere analysis from long-read sequencing data (Nanopore and PacBio).

## Table of Contents

- [Installation](#installation)
- [Quick Start](#quick-start)
  - [Using Presets](#using-presets)
- [Commands](#commands)
  - [telox-asm](#telox-asm)
  - [telox-bulk](#telox-bulk)
  - [telox-length](#telox-length)
  - [telox-variants](#telox-variants)
  - [telox-methyl](#telox-methyl)
- [Parameters](#parameters)
- [Outputs](#outputs)
- [License](#license)
- [Citation](#citation)

## Installation

**Using Conda (Recommended)**

```
conda install -c bioconda::teloxplorer
```

## Quick Start

Below is a standard command for running teloxplorer:

```
teloxplorer --preset human -fq hg001.fq.gz -r chm13v2.0.fa -R tel_repeats.human.txt --mm2-preset "-as map-ont" -o hg001 --outdir hg001 -t 12 --plot
```

Parameters

### Using Presets

To simplify configuration, we provide optimized parameter presets for several common species.

**Recommended Presets:**

- `--preset human` (Human)
- `--preset yeast` (Yeast)
- `--preset mouse` (Mouse)
- `--preset arabidopsis` (Arabidopsis)

If your species is not on this list, we recommend comparing its telomere repeat sequence to those of the species above. Choose the preset from the species with the most similar repeat sequence as a starting point.

For example, many mammals (e.g., cow, dog) share the same `TTAGGG` repeat as humans, making `--preset human` a good initial choice.

**Preset Parameter Overrides:**

The `--preset` flag automatically sets default values for many parameters. You can **override** any of these defaults by explicitly specifying the parameter in your command.

This is particularly useful for fine-tuning the analysis for "other species." The most common parameters to adjust are:

- `--genome_subtel_range`: Size of the subtelomeric range (bp) from the chromosome ends. Telomeres within this range are classified as `terminal`, while those outside are `internal`.
- `--min_tel_freq`: Sets the minimum frequency of telomeric repeats required to classify a read as telomeric.
- `--bloom_options`: BLOOM merging options, see `BLOOM --help`.

```
# Override Example
teloxplorer --preset human --min_tel_freq 0.4 ...
```

## Commands

### telox-asm

Assembly telomere length and boundary detection.

```
telox-asm --preset $psecies -asm $asm -R $tel_repeats --outdir $outdir --out-prefix $prefix
```

### telox-bulk

Bulk telomere length eastimation (alignment free).

```
telox-bulk --preset $psecies -asm $asm -R $tel_repeats --outdir $outdir --out-prefix $prefix
```

### telox-length

Chromosome-specific telomere length estimation.

```
telox-length
```

### telox-variants

Chromosome-specific telomere variant repeat analysis

```
telox-variants
```

### telox-methyl

Chromosome-specific telomere methylation analysis

```
telox-methyl \
    --preset $species \
    --tel-reads $sample.chromtel_length.tsv \
    --modBAM $bam \
    -r $ref_genome \
    -R $tel_repeats \
    -o $sample \
    --outdir $batch_id \
    -t $threads \
    --plot
```

## Parameters

A detailed list of all command-line arguments: 
- `--long-read-fastq`: Input long-read FASTQ file (gzipped or plain)
- `-r, --ref-genome`: Reference genome FASTA file
- `--outdir`: Output directory
- `-o, --out-prefix`: Prefix for all output files
- `--preset`: Use a built-in parameter preset (human, yeast, etc.)
- `-R, --tel-repeats`: Telomere repeat definition file (one per line: arm<TAB>repeat, regex supported)
- `-B, --genome_subtel_range`: Size of the subtelomeric region (bp) from the chromosome ends. Telomeres within this range are classified as `terminal`, while those outside are `internal`. (Varies by preset)
- `--mm2-preset`: minimap2 preset (e.g., map-pb, map-ont)
- `-q, --min-mapq`: Minimum mapping quality [default: 20]
- `-f, --min_tel_freq`: Frequency threshold for telomere segment definition (Varies by preset)
- `-rl, --min-read-len`: Minimum read length (bp) [default: 1000]
- `--primary-merge`: Merge adjacent telomere segments pre-BLOOM [default: no]
- `--max-mismatch`: Mismatch tolerance for telomere re-labeling [default: 2]
- `--max-initial-offset`: Max non-telomere length at read start [default: 200]
- `-tl, --min-tel-len`: Minimum telomere length [default: 100]
- `-sl, --min-subtel-len`: Minimum sub-telomere length [default: 200]
- `--bloom_options`: BLOOM merging parameters (Varies by preset)
- `--plot`: Generate telomere length plot
- `-H, --height`: Plot height in inches
- `-W, --width`: Plot width in inches
- `-k, --kmers`: K-mer sizes for repeat unit identification
- `-t, --threads`: Number of threads to use

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
  - pandas
  - pysam
  - regex
  - pyfastx
  - natsort
  - plotnine

## License

This project is licensed under the MIT License - see the LICENSE file for details.

## Citation












