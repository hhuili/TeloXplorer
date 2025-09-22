![image](https://github.com/hhuili/TeloXplorer/blob/main/logo/logo.svg)
# Teloxplorer
Teloxplorer is a pipeline for chromosome-specific telomere analysis from long-read sequencing data.

# Installation
```
# install with conda
conda install -c bioconda teloxplorer

# 
```

# Quick Start
```
# for human telomere
teloxplorer \
    --preset Human \
    --long-read-fastq human_ont_test.fastq.gz \
    --ref-genome ref.fasta \
    --tel-repeats tel_repeats.Human.txt \
    --minimap2-preset "-ax map-ont" \
    --out-prefix ONT_demo \
    --outdir Human_telo \
    --threads 4

# for mouse telomere
teloxplorer \
    --preset Mouse \
    --long-read-fastq mouse_ont_test.fastq.gz \
    --ref-genome ref.fasta \
    --tel-repeats tel_repeats.Mouse.txt \
    --minimap2-preset "-ax map-ont" \
    --out-prefix mouse_ont_test \
    --outdir Mouse_telo \
    --threads 4

# for yeast telomere
teloxplorer \
    --preset Yeast \
    --long-read-fastq yeast_ont_test.fastq.gz \
    --ref-genome ref.fasta \
    --tel-repeats tel_repeats.Yeast.txt \
    --minimap2-preset "-ax map-ont" \
    --out-prefix yeast_ont_test \
    --outdir Yeast_telo \
    --threads 4
```

# Docs

# Dependencies
- Minimap2 (https://github.com/lh3/minimap2)
- Python 3 (with the following auto-installed packages)
  - numpy
  - pysam
  - regex
# Citation

# Acknowledgments

