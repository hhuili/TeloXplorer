# Copyright (C) 2025 Huihui Li <hhui.li@outlook.com>. Licensed under GNU GPL v3.0.

import os
os.environ["OMP_NUM_THREADS"] = "1"
os.environ["OPENBLAS_NUM_THREADS"] = "1"
os.environ["MKL_NUM_THREADS"] = "1"
os.environ["VECLIB_MAXIMUM_THREADS"] = "1"
os.environ["NUMEXPR_NUM_THREADS"] = "1"

import sys
import re
import logging
from pathlib import Path
from typing import Optional, Annotated, Literal
from dataclasses import fields, is_dataclass
import typer
import click

from . import __version__
from . import logger as log_utils

logging.getLogger('fontTools').setLevel(logging.WARNING)
logging.getLogger('matplotlib').setLevel(logging.WARNING)

app = typer.Typer(
    name="telox",
    help="TeloXplorer: chromosome-end-resolved telomere analysis for long-read sequencing data",
    add_completion=False,
    no_args_is_help=True,
    pretty_exceptions_show_locals=False,
    add_help_option=False
)

def help_callback(ctx: typer.Context, value: bool):
    if value:
        typer.echo(ctx.get_help())
        raise typer.Exit()

def version_callback(value: bool):
    if value:
        typer.echo(f"teloxplorer v{__version__}")
        raise typer.Exit()

class BpSizeType(click.ParamType):
    name = "size"
    def convert(self, value, param, ctx):
        if value is None: return None
        if isinstance(value, int): return value

        size_str = str(value).strip().upper()
        match = re.match(r'^([0-9\.]+)\s*(K|KB|M|MB|G|GB|B|BP)?$', size_str)
        if not match:
            self.fail(f"Invalid size format: '{value}'. Expected format like '500k' or '0.5M'.", param, ctx)

        number_part = float(match.group(1))
        unit_part = match.group(2) or ""

        if unit_part.startswith('K'): multiplier = 10**3
        elif unit_part.startswith('M'): multiplier = 10**6
        elif unit_part.startswith('G'): multiplier = 10**9
        else: multiplier = 1

        return int(number_part * multiplier)

BP_SIZE = BpSizeType()


def parse_haplotype_filter(value: Optional[str]):
    if value is not None:
        from .viz_tvr_hap import HaplotypeFilter

        try:
            return HaplotypeFilter.parse(value)
        except (TypeError, ValueError) as exc:
            raise typer.BadParameter(str(exc)) from exc
    return None


# ==========================================
# OPTION DEFINITIONS
# ==========================================

HelpOpt = Annotated[Optional[bool], typer.Option("--help", "-h", help="Show this help and exit", callback=help_callback, is_eager=True, rich_help_panel="General")]
VersionOpt = Annotated[Optional[bool], typer.Option("--version", "-v", help="Show version and exit", callback=version_callback, is_eager=True, rich_help_panel="General")]

# --- Shared Options ---
PresetOpt = Annotated[Literal["human", "human-r9", "mouse", "yeast", "arabidopsis", "arabidopsis-r9"], typer.Option("--preset", "-x", help="Species preset", rich_help_panel="General")]
AsmPresetOpt = Annotated[Literal["human", "mouse", "yeast", "arabidopsis"], typer.Option("--preset", "-x", help="Species preset", rich_help_panel="General")]
OutDirOpt = Annotated[Path, typer.Option("--outdir", "-o", metavar="DIR", help="Output directory", rich_help_panel="General")]
PrefixOpt = Annotated[str, typer.Option("--prefix", "-p", metavar="STR", help="Output filename prefix", rich_help_panel="General")]
ThreadsOpt = Annotated[int, typer.Option("--threads", "-t", metavar="INT", help="Number of threads", rich_help_panel="General")]
DebugOpt = Annotated[bool, typer.Option("--debug", help="Enable debug logging and retain temporary intermediates", rich_help_panel="General")]
MiniLogOpt = Annotated[bool, typer.Option("--mini-log", help="Use minimalist logging style", rich_help_panel="General", show_default=True)]
ChunkSizeOpt = Annotated[int, typer.Option("--chunk-size", "-c", metavar="INT", help="Number of reads processed per chunk", rich_help_panel="General")]
StartStepOpt = Annotated[int, typer.Option("--start-step", min=1, max=5, metavar="INT", help="Resume at step 1-5: extract, align, length, TVR, or methyl", rich_help_panel="General")]

FastqOpt = Annotated[Optional[Path], typer.Option("--fastq", "-i", metavar="FILE", help="Input long-read FASTQ file (supports .gz)", rich_help_panel="General")]
ReadBAMOpt = Annotated[Optional[Path], typer.Option("--bam", "-b", metavar="FILE", help="Input unaligned long-read BAM file", rich_help_panel="General")]
ModBAMOpt = Annotated[Optional[Path], typer.Option("--modbam", metavar="FILE", help="Input modified BAM file (MM/ML tags)", rich_help_panel="General")]
RefGenomeOpt = Annotated[Path, typer.Option("--ref", "-r", metavar="FILE", help="Reference genome FASTA file", rich_help_panel="General")]
ChromMapOpt = Annotated[Optional[Path], typer.Option("--chrom-map", metavar="FILE", help="TSV mapping reference contigs via chrom, base_chrom, and chrom_phase columns", rich_help_panel="General")]
MotifOpt = Annotated[Optional[str], typer.Option("--motif", "-m", metavar="STR", help="Override the preset telomere motif or regex", rich_help_panel="General")]

# --- Align Options ---
MM2OptsOpt = Annotated[str, typer.Option("--mm2-opts", metavar="STR", help="Alignment options passed to minimap2", rich_help_panel="Align")]
MinReadQOpt = Annotated[float, typer.Option("--min-read-qual", metavar="FLOAT", help="Minimum average Phred read quality", rich_help_panel="Align")]
MinMAPQOpt = Annotated[int, typer.Option("--min-mapq", "-q", metavar="INT", help="Minimum mapping quality", rich_help_panel="Align")]

# --- Length Options ---
AlignedBAMOpt = Annotated[Path, typer.Option("--bam", "-b", metavar="FILE", help="Input aligned BAM file", rich_help_panel="General")]
TerminalRangeOpt = Annotated[Optional[int], typer.Option("--terminal-range", "-T", click_type=BP_SIZE, metavar="SIZE", help="[Preset] Terminal-region span from each chromosome end; accepts bp suffixes (e.g. 500k)", rich_help_panel="Length")]
MinReadLenOpt = Annotated[int, typer.Option("--min-read-len", metavar="INT", help="Minimum read length (bp)", rich_help_panel="Length")]
BaselineDensityOpt = Annotated[Optional[float], typer.Option("--baseline-density", "-d", metavar="FLOAT", help="[Preset] Minimum telomeric-motif density used to classify a segment", rich_help_panel="Length")]
FuzzyMismatchOpt = Annotated[Optional[int], typer.Option("--fuzzy-mismatch", metavar="INT", help="[Preset] Maximum mismatches per motif for fuzzy density", rich_help_panel="Length")]
MaxOffsetOpt = Annotated[Optional[int], typer.Option("--max-offset", metavar="INT", help="[Preset] Maximum distal non-telomeric offset in bp; -1 disables the limit", rich_help_panel="Length")]
MinTelLenOpt = Annotated[Optional[int], typer.Option("--min-tel-len", metavar="INT", help="[Preset] Minimum estimated telomere length in bp", rich_help_panel="Length")]
MinAnchorLenOpt = Annotated[Optional[int], typer.Option("--min-anchor-len", metavar="INT", help="[Preset] Minimum subtelomeric alignment length in bp", rich_help_panel="Length")]
PreMergeOpt = Annotated[Optional[Literal["yes", "no"]], typer.Option("--pre-merge", help="[Preset] Merge adjacent motif segments before boundary detection", rich_help_panel="Length")]
BloomPwidthOpt = Annotated[Optional[int], typer.Option("--bloom-pwidth", metavar="INT", help="[Preset] Horizontal smoothing tolerance in bp", rich_help_panel="Length")]
BloomYdisOpt = Annotated[Optional[float], typer.Option("--bloom-ydis", metavar="FLOAT", help="[Preset] Vertical density-difference tolerance", rich_help_panel="Length")]
MinTelQOpt = Annotated[float, typer.Option("--min-tel-qual", metavar="FLOAT", help="Minimum average Phred quality of the inferred telomeric tract. Lower-quality sequences are excluded from downstream TVR analysis", rich_help_panel="Length")]
PlotLengthOpt = Annotated[bool, typer.Option("--plot-length", help="Enable telomere length plot", rich_help_panel="Length", show_default=True)]

# --- ASM Options ---
ASMOpt = Annotated[Path, typer.Option("--assembly", "-a", metavar="FILE", help="Input assembly FASTA", rich_help_panel="General")]
TerminalLengthOpt = Annotated[Optional[int], typer.Option("--terminal-length", click_type=BP_SIZE, metavar="SIZE", help="[Preset] Sequence span extracted from each chromosome end; accepts bp suffixes (e.g. 500k)", rich_help_panel="General")]
MaskedFastaOpt = Annotated[Optional[Path], typer.Option("--masked-fasta", metavar="FILE", help="Output assembly FASTA with masked telomeric regions", rich_help_panel="General")]

# --- TVR Options ---
TelSeqOpt = Annotated[Path, typer.Option("--telseq", "-i", metavar="FILE", help="Input telomere sequence fasta (telox-length output)", rich_help_panel="General")]
KmersOpt = Annotated[Optional[str], typer.Option("--kmers", "-k", metavar="STR", help="[Preset] Comma-separated k-mer lengths for TVR decomposition", rich_help_panel="TVR")]
MinClusterSizeOpt = Annotated[Optional[int], typer.Option("--min-cluster-size", metavar="INT", help="[Preset] Minimum reads required for a valid cluster", rich_help_panel="TVR")]
ClusterEpsilonOpt = Annotated[Optional[float], typer.Option("--cluster-epsilon", metavar="FLOAT", help="[Preset] Distance threshold for HDBSCAN sub-cluster merging", rich_help_panel="TVR")]
MaxDepthOpt = Annotated[int, typer.Option("--max-depth", metavar="INT", help="Maximum read depth per cluster used for POA consensus building", rich_help_panel="TVR")]
IsAsmOpt = Annotated[bool, typer.Option("--assembly-mode", help="Accept assembly-derived haplotypes", rich_help_panel="TVR")]

# --- Methyl Options ---
SubtelExtensionOpt = Annotated[int, typer.Option("--subtel-extension", "-e", click_type=BP_SIZE, metavar="SIZE", help="Distance extending inward into subtelomere, accepts bp suffixes (e.g., 10k)", rich_help_panel="Methyl", show_default="10k")]
CThresholdOpt = Annotated[float, typer.Option("--c-threshold", metavar="FLOAT", help="Minimum probability for unmodified Cytosine", rich_help_panel="Methyl")]
ModMThresholdOpt = Annotated[float, typer.Option("--mod-m-threshold", metavar="FLOAT", help="Minimum probability for 5mC (m)", rich_help_panel="Methyl")]
ModHThresholdOpt = Annotated[float, typer.Option("--mod-h-threshold", metavar="FLOAT", help="Minimum probability for 5hmC (h)", rich_help_panel="Methyl")]

# --- Plotting General Options ---
PlotMotifOpt = Annotated[Optional[str], typer.Option("--motif", "-m", metavar="STR", help="Canonical motif treated as telomeric background", rich_help_panel="General")]
PlotWidthOpt = Annotated[Optional[float], typer.Option("--width", "-W", metavar="FLOAT", help="Figure width in inches", rich_help_panel="Plotting")]
PlotHeightOpt = Annotated[Optional[float], typer.Option("--height", "-H", metavar="FLOAT", help="Figure height in inches", rich_help_panel="Plotting")]
FigureFormatOpt = Annotated[Optional[Literal["pdf", "png", "svg"]], typer.Option("--format", "--output-format", metavar="FORMAT", help="Figure format: pdf (default), png, or svg.", rich_help_panel="Plotting")]
ChromsOpt = Annotated[Optional[str], typer.Option("--chroms", metavar="STR", help="Comma-separated list of chromosomes to plot", rich_help_panel="Plotting")]
XlimStrOpt = Annotated[Optional[str], typer.Option("--xlim", metavar="MIN_BP,MAX_BP|none", help="X-axis limits as MIN_BP,MAX_BP relative to the T-S boundary; 'none' uses each panel's full data span", rich_help_panel="Plotting")]
XlimFloatOpt = Annotated[Optional[str], typer.Option("--xlim", metavar="BP|none", help="Maximum distance from the T-S boundary to plot, in bp; 'none' uses each panel's full haplotype span", rich_help_panel="Plotting")]
TVRSampleSheetOpt = Annotated[Optional[Path], typer.Option("--sample-sheet", "-s", metavar="FILE", help="Two-column SAMPLE TVR_FILE sheet", rich_help_panel="Plotting")]
HapSampleSheetOpt = Annotated[Optional[Path], typer.Option("--sample-sheet", "-s", metavar="FILE", help="Two-column SAMPLE CONSENSUS_FILE sheet", rich_help_panel="Plotting")]
SampleLabelOpt = Annotated[Optional[str], typer.Option("--sample-label", metavar="STR", help="Title prefix for generated plots", rich_help_panel="Plotting")]
HapRasterizeOpt = Annotated[Optional[bool], typer.Option("--rasterize/--no-rasterize", help="Force or disable rasterization; by default only heatmaps with at least 1,000 haplotypes are rasterized", rich_help_panel="Plotting")]
LengthRasterizeOpt = Annotated[Optional[bool], typer.Option("--rasterize/--no-rasterize",help="Force or disable point rasterization; otherwise rasterize automatically above 20,000 reads",rich_help_panel="Length")]

# --- Plot-Length Options ---
LengthTSVOpt = Annotated[Path, typer.Option("--input", "-i", metavar="FILE", help="Input telomere length file", rich_help_panel="Plotting")]
ArmColorsOpt = Annotated[Optional[str], typer.Option("--arm-colors", metavar="STR", help="Custom arm colors (e.g. '#45732b,#f97644')", rich_help_panel="Plotting")]

# --- Plot-Methyl Options ---
ModsOpt = Annotated[Optional[Path], typer.Option("--mods", metavar="FILE", help="Input telomere modification file", rich_help_panel="Plotting")]
BinSizeOpt = Annotated[Optional[int], typer.Option("--bin-size", metavar="INT", help="Zero-anchored methylation bin width in bp", rich_help_panel="Plotting")]
ShowHeatmapOpt = Annotated[Optional[bool], typer.Option("--show-heatmap/--no-heatmap", help="Show the per-arm methylation heatmap", rich_help_panel="Plotting")]
HeatmapCmapOpt = Annotated[Optional[str], typer.Option("--heatmap-cmap", metavar="STR", help="Matplotlib colormap (e.g. 'Reds' or 'viridis')", rich_help_panel="Plotting")]
HeatmapSmoothWindowOpt = Annotated[Optional[int], typer.Option("--smooth-window", metavar="INT", help="Gaussian smoothing window in bp", rich_help_panel="Plotting")]
ShowSDOpt = Annotated[Optional[bool], typer.Option("--show-sd/--no-sd", help="Show chromosome-end standard deviation shading", rich_help_panel="Plotting")]
MinBinReadsOpt = Annotated[Optional[int], typer.Option("--min-bin-reads", metavar="INT", help="Minimum contributing reads per chromosome-arm bin", rich_help_panel="Plotting")]
MinBinArmsOpt = Annotated[Optional[int], typer.Option("--min-bin-arms", metavar="INT", help="Minimum chromosome arms required for an aggregate bin", rich_help_panel="Plotting")]

# --- Plot-Reads Options ---
TVROpt = Annotated[Optional[Path], typer.Option("--tvr", metavar="FILE", help="Input TVR file", rich_help_panel="Plotting")]
ConsensusOpt = Annotated[Optional[Path], typer.Option("--consensus", metavar="FILE", help="Input TVR haplotype consensus file", rich_help_panel="Plotting")]
TopTVROpt = Annotated[Optional[int], typer.Option("--top-tvr", metavar="INT", help="Maximum TVR motifs colored per sample and arm", rich_help_panel="Plotting")]
TVRColorsOpt = Annotated[Optional[Path], typer.Option("--tvr-colors", metavar="FILE", help="MOTIF<TAB>COLOR file or a generated *.tvr_colors.tsv map", rich_help_panel="Plotting")]
ReadThicknessOpt = Annotated[Optional[float], typer.Option("--read-thickness", metavar="FLOAT", help="Read-track thickness as a fraction of row spacing (0, 1]", rich_help_panel="Plotting")]
ConsensusThicknessOpt = Annotated[Optional[float], typer.Option("--consensus-thickness", metavar="FLOAT", help="Positive consensus-track thickness relative to row spacing; gaps expand above 1", rich_help_panel="Plotting")]
ShowOutliersOpt = Annotated[Optional[bool], typer.Option("--show-outliers", help="Show reads assigned to the 'Outlier' cluster", rich_help_panel="Plotting")]
ShowUnmethylatedOpt = Annotated[Optional[bool], typer.Option("--show-unmethylated", help="Show unmethylated cytosine calls", rich_help_panel="Plotting")]
ReadsPoolChromsOpt = Annotated[Optional[bool], typer.Option("--pool-chroms", help="Pool selected chromosomes by phase and arm", rich_help_panel="Plotting")]

# --- Plot-TVR-Hap Options ---
NoHeatmapOpt = Annotated[Optional[bool], typer.Option("--no-heatmap", help="Hide the TVR similarity heatmap", rich_help_panel="Plotting")]
ClusterRowsOpt = Annotated[Optional[bool], typer.Option("--cluster-rows/--no-cluster-rows", help="Cluster TVR haplotypes hierarchically", rich_help_panel="Plotting")]
ClusterProximalBpOpt = Annotated[Optional[int], typer.Option("--cluster-proximal-bp", metavar="INT", help="Proximal span used for pairwise haplotype similarity, in bp", rich_help_panel="Plotting")]
RowLabelColumnsOpt = Annotated[Optional[str], typer.Option("--row-label-columns", metavar="COLUMNS", help="Comma-separated columns for row labels. Default: 'auto' ('chrom,arm' / 'sample_id,chrom,arm'). Use 'none' to hide row labels.", rich_help_panel="Plotting")]
AnnotationColumnsOpt = Annotated[Optional[str], typer.Option("--annotation-columns", metavar="COLUMNS", help="Comma-separated columns from the haplotype input (e.g. 'chrom,chrom_phase'); default: 'chrom'; use 'none' to hide input-derived annotation tracks", rich_help_panel="Plotting")]
AnnoFileOpt = Annotated[Optional[Path], typer.Option("--annotation-file", metavar="FILE", help="Annotation table with an ID column followed by categorical fields", rich_help_panel="Plotting")]
AnnoColorsOpt = Annotated[Optional[Path], typer.Option("--annotation-colors", metavar="FILE", help="Headerless VALUE<TAB>COLOR file that overrides categorical colors", rich_help_panel="Plotting")]
LegendNcolOpt = Annotated[Optional[int], typer.Option("--legend-ncol", metavar="INT", help="Legend columns, 0 selects automatically", rich_help_panel="Plotting")]
WriteSimilarityScoreOpt = Annotated[bool, typer.Option("--write-similarity-score", help="Export pairwise haplotype similarity scores as TSV", rich_help_panel="Plotting")]
TrackThicknessOpt = Annotated[Optional[float], typer.Option("--track-thickness", metavar="FLOAT", help="Track thickness as a fraction of row spacing (0, 1]", rich_help_panel="Plotting")]
HapPoolChromsOpt = Annotated[Optional[bool], typer.Option("--pool-chroms", help="Pool selected chromosomes within each arm in multi-sample plots", rich_help_panel="Plotting")]
PoolArmsOpt = Annotated[Optional[bool], typer.Option("--pool-arms", help="Pool p/q or L/R arms into one plot", rich_help_panel="Plotting")]
HapFilterOpt = Annotated[Optional[str], typer.Option(
    "--hap-filter",
    callback=parse_haplotype_filter,
    metavar="HPRC2|MAX,SINGLE,MULTI",
    help=(
        "Chromosome-end haplotype filter: HPRC2 equals 2,0.4,0.2; "
        "custom values are maximum haplotypes, single-haplotype minimum "
        "frequency, and multi-haplotype minimum frequency"
    ),
    rich_help_panel="Haplotype filtering",
)]


# ==========================================
# HELPER FUNCTIONS
# ==========================================
def initialize_env(ctx: typer.Context, outdir: Path, prefix: str, debug: bool):
    ctx.ensure_object(dict)
    ctx.obj["outdir"] = outdir
    ctx.obj["prefix"] = prefix
    ctx.obj["debug"] = debug
    if not outdir.exists():
        outdir.mkdir(parents=True, exist_ok=True)

def require_exactly_one(**options):
    provided = [name for name, value in options.items() if value is not None]
    if len(provided) == 1:
        return
    names = ", ".join(options)
    typer.secho(
        f"Error: Provide exactly one of: {names}.",
        fg=typer.colors.RED,
        err=True,
    )
    raise typer.Exit(code=1)

def build_config(config_cls, preset_name: str, ctx: typer.Context):
    if not is_dataclass(config_cls):
        raise ValueError(f"{config_cls} must be a dataclass")

    valid_keys = {f.name for f in fields(config_cls)}
    params = dict(ctx.params)
    if "canonical_motif" in valid_keys and params.get("motif") is not None:
        params["canonical_motif"] = params["motif"]
    filtered_kwargs = {
        k: v for k, v in params.items()
        if k in valid_keys and v is not None
    }

    return config_cls.from_params(preset_name=preset_name, **filtered_kwargs)


def get_input_str(fastq: Optional[Path], bam: Optional[Path], modbam: Optional[Path]) -> str:
    if modbam: return f"modbam={modbam.name}"
    if fastq: return f"fastq={fastq.name}"
    if bam: return f"bam={bam.name}"
    return "Unknown"

# ==========================================
# MAIN CALLBACK
# ==========================================
@app.callback()
def main_callback(ctx: typer.Context, help_flag: HelpOpt = None, version_flag: VersionOpt = None,):
    """
    TeloXplorer: chromosome-end-resolved telomere analysis tool
    """
    pass

# ==========================================
# COMMAND: RUN (Pipeline)
# ==========================================
@app.command()
def run(
    ctx: typer.Context,
    preset: PresetOpt = "human",
    fastq: FastqOpt = None,
    modbam: ModBAMOpt = None,
    bam_reads: ReadBAMOpt = None,
    ref: RefGenomeOpt = ...,
    motif: MotifOpt = None,
    min_read_qual: MinReadQOpt = 0.0,
    mm2_opts: MM2OptsOpt = "-ax map-ont",
    min_mapq: MinMAPQOpt = 20,
    min_read_len: MinReadLenOpt = 1000,
    terminal_range: TerminalRangeOpt = None,
    baseline_density: BaselineDensityOpt = None,
    fuzzy_mismatch: FuzzyMismatchOpt = None,
    max_offset: MaxOffsetOpt = None,
    min_tel_len: MinTelLenOpt = None,
    min_anchor_len: MinAnchorLenOpt = None,
    pre_merge: PreMergeOpt = None,
    bloom_pwidth: BloomPwidthOpt = None,
    bloom_ydis: BloomYdisOpt = None,
    min_tel_qual: MinTelQOpt = 0.0,
    plot_length: PlotLengthOpt = False,
    width: Annotated[float, typer.Option("--width", "-W", help="Figure width in inches", rich_help_panel="Length")] = 8.0,
    height: Annotated[float, typer.Option("--height", "-H", help="Figure height in inches", rich_help_panel="Length")] = 3.0,
    kmers: KmersOpt = None,
    max_depth: MaxDepthOpt = 150,
    min_cluster_size: MinClusterSizeOpt = None,
    cluster_epsilon: ClusterEpsilonOpt = None,
    subtel_extension: SubtelExtensionOpt = 10000,
    c_threshold: CThresholdOpt = 0.66,
    mod_m_threshold: ModMThresholdOpt = 0.90,
    mod_h_threshold: ModHThresholdOpt = 0.90,
    outdir: OutDirOpt = Path("."),
    prefix: PrefixOpt = "telox",
    chunk_size: ChunkSizeOpt = 2000,
    threads: ThreadsOpt = 4,
    start_step: StartStepOpt = 1,
    mini_log: MiniLogOpt = False,
    debug: DebugOpt = False,
    help_flag: HelpOpt = None,
):
    """Run the complete end-to-end telomere pipeline"""
    initialize_env(ctx, outdir, prefix, debug)
    require_exactly_one(**{"--fastq": fastq, "--bam": bam_reads, "--modbam": modbam})

    from . import telox_align, telox_length, telox_tvr, telox_methyl, viz_length

    align_config = build_config(telox_align.AlignConfig, preset, ctx)
    length_config = build_config(telox_length.LengthConfig, preset, ctx)
    tvr_config = build_config(telox_tvr.TVRConfig, preset, ctx)
    plot_config = build_config(viz_length.LengthPlotter, None, ctx)
    methyl_config = build_config(telox_methyl.MethylConfig, None, ctx)

    global_total_steps = 5 if modbam else 4
    logger = log_utils.init_logger(outdir, prefix, debug, mini_log, total_steps=global_total_steps, version=__version__, start_step=1)

    logger.log_header(cmd_args=sys.argv, input_file=get_input_str(fastq, bam_reads, modbam), ref_file=ref.name, preset_name=preset, motif=align_config.motif)

    dirs = {k: outdir / f"0{i+1}_{k}" for i, k in enumerate(["reads", "align", "length", "tvr", "methyl"])}
    for d in dirs.values(): d.mkdir(parents=True, exist_ok=True)
    tel_like_reads = dirs["reads"] / f'{prefix}.tel_like.fastq'

    if start_step <= 1:
        logger.start_step("Extract telomere-like reads")
        telox_align.run_telogrep(fastq_input=fastq, modbam_input=modbam, bam_input=bam_reads, motif=align_config.motif, fastq_output=tel_like_reads, config=align_config, threads=threads, logger=logger)
        logger.end_step()
    else: logger.skip_step("Extract telomere-like reads")

    if start_step <= 2:
        logger.start_step("Align to the reference")
        telox_align.align_reads(fastq_input=tel_like_reads, ref_genome=ref, outdir=dirs["align"], prefix=prefix, config=align_config, threads=threads, logger=logger)
        logger.end_step()
    else: logger.skip_step("Align to the reference")

    if start_step <= 3:
        aligned_bam = dirs["align"] / f"{prefix}.sort.bam"
        if not aligned_bam.exists():
            typer.secho(f"Error: BAM file not found: {aligned_bam}", fg=typer.colors.RED, err=True)
            raise typer.Abort()
        logger.start_step("Estimate telomere lengths")
        telox_length.estimate_length(bam_file=aligned_bam, ref_genome=ref, outdir=dirs["length"], prefix=prefix, config=length_config, plot_config=plot_config, threads=threads, logger=logger)
        logger.end_step()
    else: logger.skip_step("Estimate telomere lengths")

    if start_step <= 4:
        telomere_seq = dirs["length"] / f"{prefix}.chromtel.seq.fasta"
        if telomere_seq.exists():
            logger.start_step("Profile TVRs")
            telox_tvr.find_tvr(telseq=telomere_seq, outdir=dirs["tvr"], prefix=f"{prefix}.chromtel", config=tvr_config, threads=threads, logger=logger)
            logger.end_step()
    else: logger.skip_step("Profile TVRs")

    if start_step <= 5 and modbam:
        telomere_fasta = dirs["length"] / f"{prefix}.chromtel.seq.fasta"
        if telomere_fasta.exists():
            logger.start_step("Profile modification calls")
            telox_methyl.find_methyl(modbam=modbam, telseq=telomere_fasta, outdir=dirs["methyl"], prefix=f"{prefix}.chromtel", config=methyl_config, threads=threads, logger=logger)
            logger.end_step()

    logger.success("teloxplorer")


# ==========================================
# COMMAND: ALIGN
# ==========================================
@app.command()
def align(
    ctx: typer.Context,
    preset: PresetOpt = "human",
    fastq: FastqOpt = None,
    bam_reads: ReadBAMOpt = None,
    modbam: ModBAMOpt = None,
    ref: RefGenomeOpt = ...,
    motif: MotifOpt = None,
    min_read_qual: MinReadQOpt = 0.0,
    mm2_opts: MM2OptsOpt = "-ax map-ont",
    min_mapq: MinMAPQOpt = 20,
    outdir: OutDirOpt = Path("."),
    prefix: PrefixOpt = "telox",
    threads: ThreadsOpt = 4,
    debug: DebugOpt = False,
    help_flag: HelpOpt = None,
):
    """Extract telomere-like reads and align to the reference genome"""
    initialize_env(ctx, outdir, prefix, debug)

    require_exactly_one(**{"--fastq": fastq, "--bam": bam_reads, "--modbam": modbam})

    from . import telox_align
    logger = log_utils.init_logger(outdir, prefix, debug, False, total_steps=2, version=__version__, start_step=1)
    config = build_config(telox_align.AlignConfig, preset, ctx)

    logger.log_header(cmd_args=sys.argv, input_file=get_input_str(fastq, bam_reads, modbam), ref_file=ref.name, preset_name=preset, motif=config.motif)

    read_dir, align_dir = outdir / "01_reads", outdir / "02_align"
    read_dir.mkdir(parents=True, exist_ok=True); align_dir.mkdir(parents=True, exist_ok=True)
    tel_like_reads = read_dir / f'{prefix}.tel_like.fastq'

    logger.start_step("Extract telomeric reads")
    telox_align.run_telogrep(fastq_input=fastq, modbam_input=modbam, bam_input=bam_reads, motif=config.motif, fastq_output=tel_like_reads, config=config, threads=threads, logger=logger)
    logger.end_step()

    logger.start_step("Align to reference")
    telox_align.align_reads(fastq_input=tel_like_reads, ref_genome=ref, outdir=align_dir, prefix=prefix, config=config, threads=threads, logger=logger)
    logger.end_step()
    logger.finish()


# ==========================================
# COMMAND: LENGTH
# ==========================================
@app.command()
def length(
    ctx: typer.Context,
    preset: PresetOpt = "human",
    bam: AlignedBAMOpt = ...,
    ref: RefGenomeOpt = ...,
    motif: MotifOpt = None,
    terminal_range: TerminalRangeOpt = None,
    min_read_len: MinReadLenOpt = 1000,
    min_tel_qual: MinTelQOpt = 0.0,
    baseline_density: BaselineDensityOpt = None,
    fuzzy_mismatch: FuzzyMismatchOpt = None,
    max_offset: MaxOffsetOpt = None,
    min_tel_len: MinTelLenOpt = None,
    min_anchor_len: MinAnchorLenOpt = None,
    pre_merge: PreMergeOpt = None,
    bloom_pwidth: BloomPwidthOpt = None,
    bloom_ydis: BloomYdisOpt = None,
    chunk_size: ChunkSizeOpt = 2000,
    plot_length: PlotLengthOpt = False,
    width: Annotated[float, typer.Option("--width", "-W", help="Figure width in inches", rich_help_panel="Length")] = 8.0,
    height: Annotated[float, typer.Option("--height", "-H", help="Figure height in inches", rich_help_panel="Length")] = 3.0,
    outdir: OutDirOpt = Path("."),
    prefix: PrefixOpt = "telox",
    threads: ThreadsOpt = 4,
    debug: DebugOpt = False,
    help_flag: HelpOpt = None,
):
    """Estimate telomere lengths from alignments"""
    initialize_env(ctx, outdir, prefix, debug)
    from . import telox_length, viz_length

    logger = log_utils.init_logger(outdir, prefix, debug, False, total_steps=1, version=__version__, start_step=1)
    length_config = build_config(telox_length.LengthConfig, preset, ctx)
    plot_config = build_config(viz_length.LengthPlotter, None, ctx)

    logger.log_header(cmd_args=sys.argv, input_file=bam, ref_file=ref.name, preset_name=preset, motif=length_config.motif)
    telox_length.estimate_length(bam_file=bam, ref_genome=ref, outdir=outdir, prefix=prefix, config=length_config, plot_config=plot_config, threads=threads, logger=logger)
    logger.finish()


# ==========================================
# COMMAND: TVR
# ==========================================
@app.command()
def tvr(
    ctx: typer.Context,
    preset: PresetOpt = "human",
    telseq: TelSeqOpt = ...,
    motif: MotifOpt = None,
    kmers: KmersOpt = None,
    min_cluster_size: MinClusterSizeOpt = None,
    max_depth: MaxDepthOpt = 150,
    cluster_epsilon: ClusterEpsilonOpt = None,
    is_asm: IsAsmOpt = False,
    outdir: OutDirOpt = Path("."),
    prefix: PrefixOpt = "telox",
    chunk_size: ChunkSizeOpt = 2000,
    threads: ThreadsOpt = 4,
    debug: DebugOpt = False,
    help_flag: HelpOpt = None,
):
    """Profile TVR composition and build haplotype consensus"""
    initialize_env(ctx, outdir, prefix, debug)
    from . import telox_tvr

    logger = log_utils.init_logger(outdir, prefix, debug, False, total_steps=1, version=__version__, start_step=1)
    tvr_config = build_config(telox_tvr.TVRConfig, preset, ctx)

    logger.log_header(cmd_args=sys.argv, input_file=telseq, preset_name=preset)
    telox_tvr.find_tvr(telseq=telseq, outdir=outdir, prefix=prefix, config=tvr_config, threads=threads, logger=logger)
    logger.finish()


# ==========================================
# COMMAND: METHYL
# ==========================================
@app.command()
def methyl(
    ctx: typer.Context,
    modbam: ModBAMOpt = ...,
    telseq: TelSeqOpt = ...,
    subtel_extension: SubtelExtensionOpt = 10000,
    c_threshold: CThresholdOpt = 0.66,
    mod_m_threshold: ModMThresholdOpt = 0.90,
    mod_h_threshold: ModHThresholdOpt = 0.90,
    chunk_size: ChunkSizeOpt = 2000,
    outdir: OutDirOpt = Path("."),
    prefix: PrefixOpt = "telox",
    threads: ThreadsOpt = 4,
    debug: DebugOpt = False,
    help_flag: HelpOpt = None,
):
    """Map 5mC/5hmC modification calls relative to T-S boundaries"""
    initialize_env(ctx, outdir, prefix, debug)
    from . import telox_methyl

    logger = log_utils.init_logger(outdir, prefix, debug, False, total_steps=1, version=__version__, start_step=1)
    methyl_config = build_config(telox_methyl.MethylConfig, None, ctx)

    logger.log_header(cmd_args=sys.argv, input_file=modbam)
    telox_methyl.find_methyl(modbam=modbam, telseq=telseq, outdir=outdir, prefix=prefix, config=methyl_config, threads=threads, logger=logger)
    logger.finish()


# ==========================================
# COMMAND: ASM
# ==========================================
@app.command()
def asm(
    ctx: typer.Context,
    preset: AsmPresetOpt = "human",
    assembly: ASMOpt = ...,
    motif: MotifOpt = None,
    terminal_length: TerminalLengthOpt = None,
    baseline_density: BaselineDensityOpt = None,
    fuzzy_mismatch: FuzzyMismatchOpt = None,
    max_offset: MaxOffsetOpt = None,
    pre_merge: PreMergeOpt = None,
    bloom_pwidth: BloomPwidthOpt = None,
    bloom_ydis: BloomYdisOpt = None,
    masked_assembly: MaskedFastaOpt = None,
    outdir: OutDirOpt = Path("."),
    prefix: PrefixOpt = "telox",
    debug: DebugOpt = False,
    help_flag: HelpOpt = None,
):
    """Estimate telomere lengths directly from genome assemblies"""
    initialize_env(ctx, outdir, prefix, debug)
    from . import telox_asm

    logger = log_utils.init_logger(outdir, prefix, debug, False, total_steps=1, version=__version__, start_step=1)
    asm_config = build_config(telox_asm.AsmConfig, preset, ctx)

    logger.log_header(cmd_args=sys.argv, input_file=assembly, preset_name=preset, motif=asm_config.motif)
    telox_asm.estimate_asm_length(assembly_input=assembly, outdir=outdir, prefix=prefix, config=asm_config, masked_assembly=str(masked_assembly) if masked_assembly else None)
    logger.finish()


# ==========================================
# COMMAND: READS
# ==========================================
@app.command()
def reads(
    ctx: typer.Context,
    preset: PresetOpt = "human",
    fastq: FastqOpt = None,
    bam_reads: ReadBAMOpt = None,
    motif: MotifOpt = None,
    min_read_len: MinReadLenOpt = 1000,
    min_read_qual: MinReadQOpt = 0.0,
    baseline_density: BaselineDensityOpt = None,
    fuzzy_mismatch: FuzzyMismatchOpt = None,
    max_offset: MaxOffsetOpt = None,
    min_tel_len: MinTelLenOpt = None,
    min_anchor_len: MinAnchorLenOpt = None,
    chunk_size: ChunkSizeOpt = 2000,
    pre_merge: PreMergeOpt = None,
    bloom_pwidth: BloomPwidthOpt = None,
    bloom_ydis: BloomYdisOpt = None,
    start_step: Annotated[int, typer.Option("--start-step", min=1, max=2, help="Start step (1:extract, 2:length)", rich_help_panel="General")] = 1,
    outdir: OutDirOpt = Path("."),
    prefix: PrefixOpt = "telox",
    debug: DebugOpt = False,
    threads: ThreadsOpt = 4,
    help_flag: HelpOpt = None,
):
    """Estimate telomere lengths directly from unaligned long reads"""
    initialize_env(ctx, outdir, prefix, debug)
    require_exactly_one(**{"--fastq": fastq, "--bam": bam_reads})
    from . import telox_align, telox_reads

    logger = log_utils.init_logger(outdir, prefix, debug, False, total_steps=2, version=__version__, start_step=1)
    read_dir, bulk_dir = outdir / "01_reads", outdir / "02_length"
    read_dir.mkdir(parents=True, exist_ok=True); bulk_dir.mkdir(parents=True, exist_ok=True)

    align_config = build_config(telox_align.AlignConfig, preset, ctx)
    bulk_config = build_config(telox_reads.BulkConfig, preset, ctx)

    logger.log_header(cmd_args=sys.argv, input_file=get_input_str(fastq, bam_reads, None), preset_name=preset, motif=bulk_config.motif)
    tel_like_reads = read_dir / f'{prefix}.tel_like.fastq'

    if start_step <= 1:
        logger.start_step("Extract telomeric reads")
        telox_align.run_telogrep(fastq_input=fastq, bam_input=bam_reads, motif=align_config.motif, fastq_output=tel_like_reads, config=align_config, threads=threads, logger=logger)
        logger.end_step()
    else: logger.skip_step("Extract telomeric reads")

    if start_step <= 2:
        logger.start_step("Estimate bulk telomere length")
        telox_reads.estimate_bulk_length(fastq_input=tel_like_reads, outdir=bulk_dir, prefix=prefix, config=bulk_config, threads=threads, logger=logger)
        logger.end_step()

    logger.finish()


# ==========================================
# COMMAND: PLOT-LENGTH
# ==========================================
@app.command(name="plot-length")
def plot_length(
    ctx: typer.Context,
    length_input: LengthTSVOpt = ...,
    chroms: ChromsOpt = None,
    arm_colors: ArmColorsOpt = None,
    width: PlotWidthOpt = 8.0,
    height: PlotHeightOpt = 3.0,
    rasterize: LengthRasterizeOpt = None,
    output_format: FigureFormatOpt = None,
    outdir: OutDirOpt = Path("."),
    prefix: PrefixOpt = "telox",
    help_flag: HelpOpt = None,
):
    """Plot chromosome-end telomere length distributions"""
    initialize_env(ctx, outdir, prefix, False)
    from . import viz_length

    logger = log_utils.init_logger(outdir, prefix, False, False, total_steps=1, version=__version__, start_step=1)
    length_config = build_config(viz_length.LengthPlotter, None, ctx)

    logger.log_header(cmd_args=sys.argv, input_file=length_input)
    viz_length.plot_length(length_input=length_input, outdir=outdir, prefix=prefix, config=length_config, logger=logger)
    logger.finish()


# ==========================================
# COMMAND: PLOT-READS
# ==========================================
@app.command(name="plot-reads")
def plot_reads(
    ctx: typer.Context,
    motif: PlotMotifOpt = None,
    tvr_file: TVROpt = None,
    consensus_file: ConsensusOpt = None,
    sample_sheet: TVRSampleSheetOpt = None,
    mods_file: ModsOpt = None,
    top_tvr: TopTVROpt = None,
    chroms: ChromsOpt = None,
    pool_chroms: ReadsPoolChromsOpt = None,
    tvr_colors: TVRColorsOpt = None,
    xlim: XlimStrOpt = None,
    read_thickness: ReadThicknessOpt = None,
    consensus_thickness: ConsensusThicknessOpt = None,
    legend_ncol: LegendNcolOpt = None,
    show_outliers: ShowOutliersOpt = None,
    show_unmethylated_c: ShowUnmethylatedOpt = None,
    sample_label: SampleLabelOpt = None,
    width: PlotWidthOpt = None,
    height: PlotHeightOpt = None,
    output_format: FigureFormatOpt = None,
    outdir: OutDirOpt = Path("."),
    prefix: PrefixOpt = "telox",
    help_flag: HelpOpt = None,
):
    """Plot single-molecule TVR architectures and methylation states"""
    initialize_env(ctx, outdir, prefix, False)
    require_exactly_one(**{"--tvr": tvr_file, "--sample-sheet": sample_sheet})
    if sample_sheet is not None and (mods_file is not None or consensus_file is not None):
        typer.secho(
            "Error: --sample-sheet mode supports TVR files only; "
            "--mods and --consensus are not supported.",
            fg=typer.colors.RED,
            err=True,
        )
        raise typer.Exit(code=1)

    from . import viz_reads

    logger = log_utils.init_logger(outdir, prefix, False, False, total_steps=1, version=__version__, start_step=1)
    user_input_path = sample_sheet if sample_sheet else tvr_file
    logger.log_header(cmd_args=sys.argv, input_file=user_input_path)

    config = build_config(viz_reads.TVRReadPlotter, None, ctx)
    viz_reads.plot_tvr_reads(
        tvr_file=tvr_file,
        sample_sheet=sample_sheet,
        mods_file=mods_file,
        consensus_file=consensus_file,
        outdir=outdir,
        prefix=prefix,
        config=config,
        logger=logger,
    )
    logger.finish()


# ==========================================
# COMMAND: PLOT-TVR-HAP
# ==========================================
@app.command(name="plot-tvr-hap")
def plot_tvr_hap(
    ctx: typer.Context,
    motif: PlotMotifOpt = None,
    consensus_file: ConsensusOpt = None,
    sample_sheet: HapSampleSheetOpt = None,
    top_tvr: TopTVROpt = None,
    chroms: ChromsOpt = None,
    pool_arms: PoolArmsOpt = None,
    pool_chroms: HapPoolChromsOpt = None,
    tvr_colors: TVRColorsOpt = None,
    hap_filter: HapFilterOpt = None,
    xlim: XlimFloatOpt = None,
    track_thickness: TrackThicknessOpt = None,
    cluster_rows: ClusterRowsOpt = None,
    cluster_proximal_bp: ClusterProximalBpOpt = None,
    no_heatmap: NoHeatmapOpt = None,
    heatmap_cmap: HeatmapCmapOpt = None,
    row_label_columns: RowLabelColumnsOpt = None,
    annotation_columns: AnnotationColumnsOpt = None,
    annotation_file: AnnoFileOpt = None,
    annotation_colors: AnnoColorsOpt = None,
    legend_ncol: LegendNcolOpt = None,
    write_similarity_score: WriteSimilarityScoreOpt = False,
    rasterize: HapRasterizeOpt = None,
    width: PlotWidthOpt = None,
    height: PlotHeightOpt = None,
    sample_label: SampleLabelOpt = None,
    output_format: FigureFormatOpt = None,
    outdir: OutDirOpt = Path("."),
    prefix: PrefixOpt = "telox",
    help_flag: HelpOpt = None,
):
    """Plot TVR haplotype consensus, clustering, and similarity heatmap"""
    initialize_env(ctx, outdir, prefix, False)
    require_exactly_one(**{"--consensus": consensus_file, "--sample-sheet": sample_sheet})

    from . import viz_tvr_hap

    logger = log_utils.init_logger(outdir, prefix, False, False, total_steps=1, version=__version__, start_step=1)
    user_input_path = sample_sheet if sample_sheet else consensus_file
    logger.log_header(cmd_args=sys.argv, input_file=user_input_path)

    config = build_config(viz_tvr_hap.TVRHapPlotter, None, ctx)
    viz_tvr_hap.plot_tvr_haps(consensus_file=consensus_file, sample_sheet=sample_sheet, outdir=outdir, prefix=prefix, config=config, write_similarity_score=write_similarity_score, logger=logger)
    logger.finish()


# ==========================================
# COMMAND: PLOT-METHYL
# ==========================================
@app.command(name="plot-methyl")
def plot_methyl(
    ctx: typer.Context,
    mods_file: ModsOpt = None,
    xlim: XlimStrOpt = None,
    bin_size: BinSizeOpt = None,
    smooth_window: HeatmapSmoothWindowOpt = None,
    show_sd: ShowSDOpt = None,
    show_heatmap: ShowHeatmapOpt = None,
    heatmap_cmap: HeatmapCmapOpt = None,
    min_bin_reads: MinBinReadsOpt = None,
    min_bin_arms: MinBinArmsOpt = None,
    width: PlotWidthOpt = None,
    height: PlotHeightOpt = None,
    output_format: FigureFormatOpt = None,
    outdir: OutDirOpt = Path("."),
    prefix: PrefixOpt = "telox",
    help_flag: HelpOpt = None,
):
    """Plot T-S boundary methylation heatmaps and aggregate 5mC frequency curves"""
    initialize_env(ctx, outdir, prefix, False)
    if mods_file is None:
        raise typer.BadParameter("--mods is required")
    from . import viz_methyl

    logger = log_utils.init_logger(outdir, prefix, False, False, total_steps=1, version=__version__, start_step=1)
    logger.log_header(cmd_args=sys.argv, input_file=mods_file)

    config = build_config(viz_methyl.MethylPlotter, None, ctx)
    viz_methyl.plot_arm_methyl(mods_file=mods_file, outdir=outdir, prefix=prefix, config=config, logger=logger)
    logger.finish()


if __name__ == "__main__":
    app()
