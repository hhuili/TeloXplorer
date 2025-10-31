![image](https://github.com/hhuili/TeloXplorer/blob/main/logo/logo.svg)
## TeloXplorer

[![Conda Version](https://img.shields.io/conda/vn/huihui_li/teloxplorer.svg?style=flat-square)](https://anaconda.org/huihui_li/teloxplorer)

Teloxplorer is a tool for chromosome-specific telomere analysis from long-read sequencing data (Nanopore and PacBio).

## Table of Contents

- [Installation](#installation)
- [Teloxplorer](#teloxplorer)
  - [Using presets](#using-presets)
  - [Parameters](#parameters)
  - [Outputs](#outputs)
- [Subcommands](#subcommands)
- [License](#license)
- [Citation](#citation)

## Installation

**Using Conda (Recommended)**

```
conda install -c bioconda::teloxplorer
```

## Teloxplorer

Below is a standard command for running teloxplorer:

```
teloxplorer --preset human \
  -fq hg001.fq.gz \
  -r chm13v2.0.fa \
  -R tel_repeats.human.txt \
  --mm2-preset "-ax map-ont" \
  -o hg001
  --outdir hg001 \
  -t 12 \
  --plot -H 3 -W 10
```

Parameters

### Using presets

To simplify configuration, we provide optimized parameter presets for several common species.

**Recommended presets:**

- `--preset human` (Human)
- `--preset yeast` (Yeast)
- `--preset mouse` (Mouse)
- `--preset arabidopsis` (Arabidopsis)

If your species is not on this list, we recommend comparing its telomere repeat sequence to those of the species above. Choose the preset from the species with the most similar repeat sequence as a starting point.

For example, many mammals (e.g., cow, dog) share the same `TTAGGG` repeat as humans, making `--preset human` a good initial choice.

**Preset parameter overrides:**

The `--preset` flag automatically sets default values for many parameters. You can **override** any of these defaults by explicitly specifying the parameter in your command.

This is particularly useful for fine-tuning the analysis for "other species." The most common parameters to adjust are:

- `--genome_subtel_range`: Size of the subtelomeric range (bp) from the chromosome ends. Telomeres within this range are classified as `terminal`, while those outside are `internal`.
- `--min_tel_freq`: Sets the minimum frequency of telomeric repeats required to classify a read as telomeric.
- `--bloom_options`: BLOOM merging options, see `BLOOM --help`.

```
# override example
teloxplorer --preset human --min_tel_freq 0.4 ...
```

### Parameters

A detailed list of all command-line arguments: 
```

```

### Outputs

teloxplorer will create an output directory specified by --outdir. Key output files include:

- `$prefix.chromtel_length.tsv`: A tab-separated file containing chomosome-specific telomere length estimates.
- `$prefix.chromtel_summary.tsv`: A summary file of chomosome-specific telomere length with statistics (mean, median, etc.).
- `$prefix.chromtel.tel_length.pdf`: (--plot) A boxplot showing the telomere length accross chromosomes.
- `$prefix.neotel_length.tsv`: A tab-separated file containing neotelomere length estimates.
- `$prefix.minitel_length.tsv`: A tab-separated file containing minitelomere length estimates.
- `$prefix.chromtel.TVR.tsv`: (--find-TVR "yes")
- `$prefix.chromtel.TVR_summary.tsv`:

## Subcommands

|Command                                                                                |Function                                                                 |Input                 |Output                                 |
|:--------------------------------------------------------------------------------------|:------------------------------------------------------------------------|:---------------------|:--------------------------------------|
|[telox-asm](https://github.com/hhuili/TeloXplorer/usage/#telox-asm)                    |Assembly telomere length estimation and boundary detection               |Assembly (fasta)      |Teloemre length/boundary (chromosome)  |
|[telox-bulk](https://github.com/hhuili/TeloXplorer/usage/#telox-bulk)                  |Bulk telomere length eastimation (alignment free)                        |Long reads (fastq)    |Teloemre length (bulk)                 |
|[telox-length](https://github.com/hhuili/TeloXplorer/usage/#telox-length)              |Chromosome-specific telomere length estimation                           |Long reads (fastq)    |Teloemre length (chromosome)           |
|[telox-variants](https://github.com/hhuili/TeloXplorer/usage/#telox-variants)          |Chromosome-specific telomere variant repeat analysis                     |Teloemre seq (fasta)  |Teloemre variants (chromosome)         |
|[telox-methyl](https://github.com/hhuili/TeloXplorer/usage/#telox-methyl)              |Chromosome-specific telomere methylation analysis                        |modBAM (bam)          |Teloemre bedMethyl (chromosome)        |
|[length-plot](https://github.com/hhuili/TeloXplorer/usage/#length-plot)                |Telomere length plot across chromosomes                                  |Teloemre length (tsv) |plot (pdf)                             |
|[methyl-plot](https://github.com/hhuili/TeloXplorer/usage/#methyl-plot)                |Methylation levels plot flanking telomere-subtelomere boundaries         |bedMethyl (bed)       |plot (pdf)                             |

### telox-asm

Assembly telomere length and boundary detection.

```
telox-asm --preset human -asm chm13v2.0.fa -R tel_repeats.human.txt --outdir chm13 -o chm13
```

### telox-bulk

Bulk telomere length eastimation (alignment free).

```
telox-bulk --preset human -fq hg001.fq.gz -R tel_repeats.human.txt --outdir hg001 -o hg001
```

### telox-length

Chromosome-specific telomere length estimation.

```
telox-length --preset human --bam hg001.sort.bam -R tel_repeats.human.txt --outdir hg001 -o hg001
```

### telox-variants

Chromosome-specific telomere variant repeat analysis

```
telox-variants -fa hg001.chromtel_seq.fa -R tel_repeats.human.txt -kmers 5,6,7 --outdir hg001 -o hg001
```

### telox-methyl

Chromosome-specific telomere methylation analysis

```
telox-methyl \
    --preset human \
    --tel-reads hg001.chromtel_length.tsv \
    --modBAM hg001.5mC_5hmC.pass.bam \
    -r chm13v2.0.fa \
    -R tel_repeats.human.txt \
    -o hg001
    --outdir hg001 \
    -t 12
    --plot -H 8 -W 8
```

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














