# Copyright (C) 2025 Huihui Li <hhui.li@outlook.com>. Licensed under GNU GPL v3.0.

from pathlib import Path
from dataclasses import dataclass
from typing import Optional, Union
import shlex

from . import telox_utils
from . import logger as log_utils

@dataclass
class AlignConfig:
    motif: Optional[str] = "TTAGGG"
    motif_C: Optional[str] = None
    motif_G: Optional[str] = None
    min_motif_copy: int = 5
    min_read_qual: float = 0.0
    mm2_opts: str = "-ax map-ont"
    min_mapq: int = 20
    debug: bool = False

    @classmethod
    def from_params(cls, preset_name=None, **kwargs):
        return telox_utils.apply_preset_and_overrides(cls, preset_name, **kwargs)

class TelomereAligner:
    def __init__(self, fastq_input, ref_genome, prefix, outdir,
                 threads, mm2_options, min_mapping_quality, logger):

        self.fastq_input = Path(fastq_input)
        self.ref_genome = Path(ref_genome)
        self.prefix = prefix
        self.outdir = Path(outdir)
        self.threads = str(threads)
        self.mm2_options = mm2_options
        self.min_mapping_quality = str(min_mapping_quality)
        self.max_secondary_alignments = '2'
        self.logger = logger

        self.minimap2_exec = 'minimap2'
        self.samtools_exec = 'samtools'

        self.bam_raw_path = self.outdir / f"{self.prefix}.bam"
        self.bam_sorted_path = self.outdir / f"{self.prefix}.sort.bam"
        self.stats_path = self.outdir / f"{self.prefix}.samstat"

        if not self.outdir.exists():
            self.outdir.mkdir(parents=True, exist_ok=True)

    def _run_mapping(self):
        command = (
            f"{self.minimap2_exec} -t {self.threads} {self.mm2_options} "
            f"-N {self.max_secondary_alignments} -Y {self.ref_genome} {self.fastq_input} | "
            f"{self.samtools_exec} view -bS -q {self.min_mapping_quality} -o {self.bam_raw_path}"
        )
        self.logger.run_cmd(command, shell=True)

    def _sort_bam(self):
        command = [
            self.samtools_exec, 'sort', 
            '-o', str(self.bam_sorted_path), 
            str(self.bam_raw_path)
        ]
        self.logger.run_cmd(command)

    def _index_bam(self):
        command = [self.samtools_exec, 'index', str(self.bam_sorted_path)]
        self.logger.run_cmd(command)

    def _generate_stats(self):
        command = (
            f"{self.samtools_exec} flagstat {self.bam_sorted_path} > {self.stats_path}"
        )
        self.logger.run_cmd(command, shell=True)

    def _cleanup(self):
        if self.bam_raw_path.exists():
            self.bam_raw_path.unlink()

    def run(self):
        self._run_mapping()
        self._sort_bam()
        self._index_bam()
        self._generate_stats()
        self._cleanup()

def run_telogrep(
    fastq_input: Optional[Union[str, Path]], 
    motif: str, 
    config: AlignConfig,
    fastq_output: Union[str, Path], 
    modbam_input: Optional[Union[str, Path]] = None,
    bam_input: Optional[Union[str, Path]] = None,
    threads: int = 4,
    logger = None
):
    if logger is None:
        logger = log_utils.get_logger()

    fastq_output = str(fastq_output)
    
    if not fastq_input and not modbam_input and not bam_input:
        raise ValueError("No input provided (fastq, bam, or modbam).")

    if not motif:
        raise ValueError("Telomere motif is not defined! Please check -m or --preset options.")

    motif_dict = telox_utils.parse_motifs(
        config.motif, 
        getattr(config, 'motif_C', None), 
        getattr(config, 'motif_G', None)
    )
    search_str_L = motif_dict["L"]
    search_str_R = motif_dict["R"]

    telogrep_args = []
    telogrep_args.extend(["-q", str(config.min_read_qual)])
    telogrep_args.extend(["-o", fastq_output])
    telogrep_args.extend(["-t", str(threads)])

    if telox_utils.is_dna_sequence(search_str_L):
        telogrep_args.extend([
            "-L", search_str_L,
            "-R", search_str_R,
            "-n", str(config.min_motif_copy)
        ])
    else:
        telogrep_args.extend([
            "-L", search_str_L,
            "-R", search_str_R,
            "-n", str(config.min_motif_copy),
            "--regex"
        ])

    bam_source = modbam_input if modbam_input else bam_input

    if bam_source:
        bam_source = str(bam_source)
        telogrep_args_str = " ".join(shlex.quote(str(arg)) for arg in telogrep_args)

        pipeline_cmd = (
            f"samtools fastq -@ {threads} {bam_source} | "
            f"telogrep -i - {telogrep_args_str}"
        )

        logger.run_cmd(pipeline_cmd, shell=True)

    else:
        fastq_input = str(fastq_input)
        cmd_list = ["telogrep", "-i", fastq_input] + telogrep_args
        logger.run_cmd(cmd_list, shell=True)


def align_reads(
    fastq_input: Union[str, Path],
    ref_genome: Union[str, Path],
    outdir: Union[str, Path],
    prefix: str,
    config: AlignConfig,
    threads: int = 4,
    logger = None
):
    logger = logger or log_utils.get_logger()
    outdir = Path(outdir)
    outdir.mkdir(parents=True, exist_ok=True)

    mapper = TelomereAligner(
        fastq_input=str(fastq_input),
        ref_genome=str(ref_genome),
        mm2_options=config.mm2_opts,
        min_mapping_quality=config.min_mapq,
        outdir=str(outdir),
        prefix=prefix,
        threads=threads,
        logger=logger
    )
    
    mapper.run()
