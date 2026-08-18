# Copyright (C) 2025 Huihui Li <hhui.li@outlook.com>. Licensed under GNU GPL v3.0.
from __future__ import annotations

from collections import Counter, defaultdict
from dataclasses import dataclass
from pathlib import Path
from typing import Iterable, Sequence

import numpy as np
import pandas as pd
from natsort import natsort_keygen, natsorted

from . import logger as log_utils
from .viz_core import (
    INK,
    MISSING_COLOR,
    STYLE,
    arm_sort_key,
    categorical_color_map,
    chrom_phase_key,
    chrom_phase_label,
    display_chrom,
    format_kb_axis,
    nature_rc,
    normalize_figure_format,
    normalize_chrom_phase,
    plot_intervals,
    require_schema,
    save_figure,
)
import matplotlib as mpl
import matplotlib.pyplot as plt


@dataclass(frozen=True, slots=True)
class MethylPlotter:
    xlim: str | Sequence[int] = "-5000,5000"
    bin_size: int = 50
    smooth_window: int = 500
    min_bin_reads: int = 5
    min_bin_arms: int = 3
    show_sd: bool = False
    show_heatmap: bool = True
    width: float = STYLE.inches(STYLE.double_column_mm)
    height: float | None = None
    heatmap_cmap: str = "Reds"
    output_format: str = "pdf"

    @classmethod
    def from_params(cls, **kwargs):
        return cls(**{k: v for k, v in kwargs.items() if k in cls.__dataclass_fields__ and v is not None})


    def resolved_xlim(self) -> tuple[int, int]:
        value = self.xlim
        if isinstance(value, str):
            parts = [part.strip() for part in value.split(",") if part.strip()]
        else:
            parts = list(value)
        if len(parts) != 2:
            raise ValueError("xlim must contain exactly two coordinates: XMIN,XMAX.")
        try:
            xmin, xmax = (int(part) for part in parts)
        except (TypeError, ValueError) as exc:
            raise ValueError("xlim coordinates must be integers.") from exc
        return xmin, xmax

    def validate(self) -> None:
        xmin, xmax = self.resolved_xlim()
        if xmin >= xmax:
            raise ValueError("xlim requires XMIN < XMAX.")
        if self.bin_size <= 0:
            raise ValueError("bin_size must be > 0.")
        if xmin % self.bin_size or xmax % self.bin_size:
            raise ValueError(
                "xlim boundaries must be integer multiples of bin_size so bins remain anchored at 0."
            )
        if self.smooth_window <= 0:
            raise ValueError("smooth_window must be > 0.")
        if self.smooth_window < self.bin_size:
            raise ValueError("smooth_window must be >= bin_size.")
        if self.min_bin_reads < 1:
            raise ValueError("min_bin_reads must be >= 1.")
        if self.min_bin_arms < 1:
            raise ValueError("min_bin_arms must be >= 1.")
        if self.width <= 0:
            raise ValueError("width must be > 0.")
        if self.height is not None and self.height <= 0:
            raise ValueError("height must be > 0.")
        normalize_figure_format(self.output_format)
        mpl.colormaps[self.heatmap_cmap]


@dataclass(slots=True)
class BinCounts:
    valid_calls: int = 0
    methylated_calls: int = 0
    contributing_reads: int = 0


@dataclass(frozen=True, slots=True)
class PhasePresentation:
    label: str | None
    curve_color: str
    row_label_color: str


REQUIRED_COLUMNS = ("read_id","chrom","chrom_phase","arm","mods_m","mods_h","valid_sites",)
NATURAL_SORT_KEY = natsort_keygen()
UNPHASED_CURVE_COLOR = "#f06449"
COLORBAR_WIDTH = 0.030


def parse_positions(value: object, *, column: str, read_id: str) -> tuple[int, ...]:
    if value is None or pd.isna(value):
        return ()

    text = str(value).strip()
    if text in {"", "."}:
        return ()

    positions: list[int] = []
    for token in text.split(","):
        token = token.strip()
        if not token or token == ".":
            continue
        if ":" in token:
            raise ValueError(
                f"{column} contains unsupported position:probability token {token!r} "
                f"for read {read_id!r}; regenerate the table with the current telox_methyl pipeline."
            )
        try:
            positions.append(int(token))
        except ValueError as exc:
            raise ValueError(f"Unparseable position {token!r} in {column} for read {read_id!r}.") from exc

    if len(positions) != len(set(positions)):
        raise ValueError(f"Duplicate positions found in {column} for read {read_id!r}.")
    return tuple(positions)


def validate_call_sets(
    *,
    read_id: str,
    valid_sites: set[int],
    mods_m: set[int],
    mods_h: set[int],
) -> None:
    extra_m = mods_m - valid_sites
    extra_h = mods_h - valid_sites
    overlap = mods_m & mods_h

    if extra_m:
        raise ValueError(
            f"Read {read_id!r}: mods_m contains positions absent from valid_sites: "
            f"{sorted(extra_m)[:8]}"
        )
    if extra_h:
        raise ValueError(
            f"Read {read_id!r}: mods_h contains positions absent from valid_sites: "
            f"{sorted(extra_h)[:8]}"
        )
    if overlap:
        raise ValueError(
            f"Read {read_id!r}: positions cannot be both 5mC and 5hmC: {sorted(overlap)[:8]}"
        )


def _bin_counts(
    positions: Iterable[int],
    *,
    xmin: int,
    xmax: int,
    bin_size: int,
) -> Counter[int]:
    return Counter(
        (position // bin_size) * bin_size
        for position in positions
        if xmin <= position < xmax
    )


def summarise_sites(
    table: pd.DataFrame,
    config: MethylPlotter,
    *,
    source: str | Path | None = None,
) -> pd.DataFrame:
    config.validate()
    xmin, xmax = config.resolved_xlim()
    require_schema(table, REQUIRED_COLUMNS, source=source)
    if table.empty:
        raise ValueError("Input methylation table is empty.")

    identity = table[["read_id", "chrom", "chrom_phase", "arm"]]
    missing_identity = identity.isna() | identity.astype(str).apply(
        lambda column: column.str.strip().eq("")
    )
    if missing_identity.any(axis=None):
        raise ValueError("read_id, chrom, chrom_phase, and arm must be non-empty for every row.")

    duplicate_reads = (
        table.loc[table["read_id"].duplicated(keep=False), "read_id"]
        .astype(str)
        .drop_duplicates()
        .head(8)
        .tolist()
    )
    if duplicate_reads:
        raise ValueError(
            "Duplicate read_id values found in methylation table; expected one row "
            f"per read: {', '.join(duplicate_reads)}"
        )

    counts: defaultdict[tuple[str, str, str, int], BinCounts] = defaultdict(BinCounts)
    arm_read_counts: defaultdict[tuple[str, str, str], int] = defaultdict(int)
    for row in table.itertuples(index=False):
        read_id = str(row.read_id)
        chrom = str(row.chrom)
        phase = normalize_chrom_phase(row.chrom_phase)
        arm = str(row.arm)
        arm_read_counts[(arm, chrom, phase)] += 1

        valid_sites = set(parse_positions(row.valid_sites, column="valid_sites", read_id=read_id))
        mods_m = set(parse_positions(row.mods_m, column="mods_m", read_id=read_id))
        mods_h = set(parse_positions(row.mods_h, column="mods_h", read_id=read_id))
        validate_call_sets(
            read_id=read_id,
            valid_sites=valid_sites,
            mods_m=mods_m,
            mods_h=mods_h,
        )

        valid_by_bin = _bin_counts(valid_sites, xmin=xmin, xmax=xmax, bin_size=config.bin_size)
        methylated_by_bin = _bin_counts(mods_m, xmin=xmin, xmax=xmax, bin_size=config.bin_size)
        for bin_start, valid_n in valid_by_bin.items():
            key = (arm, chrom, phase, int(bin_start))
            cell = counts[key]
            cell.valid_calls += int(valid_n)
            cell.methylated_calls += int(methylated_by_bin.get(bin_start, 0))
            cell.contributing_reads += 1

    if not counts:
        raise ValueError(f"No valid cytosine calls overlap [{xmin}, {xmax}) bp.")

    records: list[dict[str, object]] = []
    for (arm, chrom, phase, bin_start), cell in counts.items():
        records.append(
            {
                "arm": arm,
                "chrom": chrom,
                "chrom_phase": phase,
                "bin_start_bp": bin_start,
                "bin_end_bp": bin_start + config.bin_size,
                "bin_center_bp": bin_start + config.bin_size / 2,
                "valid_calls": cell.valid_calls,
                "methylated_calls": cell.methylated_calls,
                "contributing_reads": cell.contributing_reads,
                "arm_read_count": arm_read_counts[(arm, chrom, phase)],
                "methylation_frequency": cell.methylated_calls / cell.valid_calls,
            }
        )

    summary = pd.DataFrame.from_records(records)
    summary["_arm_order"] = summary["arm"].map(arm_sort_key)
    summary["_chrom_order"] = summary["chrom"].astype(str).map(NATURAL_SORT_KEY)
    summary["_phase_order"] = summary["chrom_phase"].astype(str).map(NATURAL_SORT_KEY)
    summary.sort_values(
        ["_arm_order", "_chrom_order", "_phase_order", "bin_start_bp"],
        inplace=True,
    )
    summary.drop(columns=["_arm_order", "_chrom_order", "_phase_order"], inplace=True)
    summary.reset_index(drop=True, inplace=True)

    return summary


def heatmap_rows(
    arm_summary: pd.DataFrame,
    config: MethylPlotter,
) -> list[tuple[str, str]]:
    eligible = arm_summary.loc[
        arm_summary["arm_read_count"].ge(config.min_bin_reads)
        & arm_summary["methylated_calls"].gt(0),
        ["chrom", "chrom_phase"],
    ].drop_duplicates()
    return natsorted(
        eligible.itertuples(index=False, name=None)
    )


def phase_presentations(phases: Sequence[object]) -> dict[str, PhasePresentation]:
    phase_values = tuple(dict.fromkeys(str(phase) for phase in phases))
    phased_values = tuple(
        phase for phase in phase_values if chrom_phase_key(phase)
    )
    phase_colors = categorical_color_map(phased_values, column="chrom_phase")
    return {
        phase: PhasePresentation(
            label=(phase_label if (phase_label := chrom_phase_label(phase)) else None),
            curve_color=phase_colors.get(phase, UNPHASED_CURVE_COLOR),
            row_label_color=phase_colors.get(phase, INK),
        )
        for phase in phase_values
    }


def methyl_row_label(chrom: object, phase: object, arm: object) -> str:
    chrom_arm = f"{display_chrom(chrom)}{arm}"
    phase_label = chrom_phase_label(phase)
    return f"{phase_label}-{chrom_arm}" if phase_label else chrom_arm


def gaussian_kernel(config: MethylPlotter) -> np.ndarray:
    half_window = config.smooth_window / 2.0
    radius_bins = max(1, int(np.floor(half_window / config.bin_size)))
    offsets_bp = np.arange(-radius_bins, radius_bins + 1, dtype=float) * config.bin_size
    sigma_bp = config.smooth_window / 4.0
    weights = np.exp(-0.5 * (offsets_bp / sigma_bp) ** 2)
    return weights / weights.sum()


def convolve_same(values: np.ndarray, kernel: np.ndarray) -> np.ndarray:
    full = np.convolve(values, kernel, mode="full")
    start = (len(kernel) - 1) // 2
    return full[start : start + len(values)]


def smooth_supported(values: np.ndarray, kernel: np.ndarray, support_mask: np.ndarray) -> np.ndarray:
    finite = np.isfinite(values)
    numerator = convolve_same(np.where(finite, values, 0.0), kernel)
    denominator = convolve_same(finite.astype(float), kernel)
    smoothed = np.divide(
        numerator,
        denominator,
        out=np.full_like(numerator, np.nan, dtype=float),
        where=denominator > 0,
    )
    smoothed[~support_mask] = np.nan
    return smoothed


def phase_frequency_matrix(
    arm_summary: pd.DataFrame,
    phase: str,
    config: MethylPlotter,
) -> tuple[np.ndarray, np.ndarray, list[str]]:
    xmin, xmax = config.resolved_xlim()
    bin_starts = np.arange(xmin, xmax, config.bin_size, dtype=int)
    phase_summary = arm_summary.loc[arm_summary["chrom_phase"] == phase]
    chromosomes = natsorted(phase_summary["chrom"].unique().tolist())
    eligible = phase_summary.loc[
        phase_summary["contributing_reads"] >= config.min_bin_reads,
        ["chrom", "bin_start_bp", "methylation_frequency"],
    ]
    matrix = (
        eligible.pivot(
            index="chrom",
            columns="bin_start_bp",
            values="methylation_frequency",
        )
        .reindex(index=chromosomes, columns=bin_starts)
        .to_numpy(dtype=float)
    )

    return matrix, bin_starts, chromosomes


def aggregate_phase_curve(
    arm_summary: pd.DataFrame,
    phase: str,
    config: MethylPlotter,
) -> pd.DataFrame:
    matrix, bin_starts, _chromosomes = phase_frequency_matrix(arm_summary, phase, config)
    n_arms = np.isfinite(matrix).sum(axis=0)
    supported = n_arms >= config.min_bin_arms

    raw_mean = np.full(len(bin_starts), np.nan, dtype=float)
    if matrix.size:
        sums = np.nansum(matrix, axis=0)
        raw_mean = np.divide(
            sums,
            n_arms,
            out=np.full(len(bin_starts), np.nan, dtype=float),
            where=n_arms > 0,
        )
    raw_mean[~supported] = np.nan

    raw_sd = np.full(len(bin_starts), np.nan, dtype=float)
    sd_supported = n_arms >= max(2, config.min_bin_arms)
    if config.show_sd and np.any(sd_supported):
        raw_sd[sd_supported] = np.nanstd(
            matrix[:, sd_supported],
            axis=0,
            ddof=1,
        )

    kernel = gaussian_kernel(config)
    mean = smooth_supported(raw_mean, kernel, supported)
    sd = smooth_supported(raw_sd, kernel, supported & np.isfinite(raw_sd)) if config.show_sd else raw_sd

    return pd.DataFrame(
        {
            "bin_start_bp": bin_starts,
            "bin_center_bp": bin_starts + config.bin_size / 2,
            "n_arms": n_arms,
            "raw_mean_frequency": raw_mean,
            "methylation_frequency": mean,
            "raw_chromosome_sd": raw_sd,
            "chromosome_sd": sd,
        }
    )


def _heatmap_cells(
    arm_summary: pd.DataFrame,
    rows: list[tuple[str, str]],
    config: MethylPlotter,
) -> pd.DataFrame:
    cmap = mpl.colormaps[config.heatmap_cmap]
    row_index = {pair: index for index, pair in enumerate(rows)}
    visible = arm_summary.loc[arm_summary["contributing_reads"] >= config.min_bin_reads].copy()
    visible["_row"] = list(
        visible[["chrom", "chrom_phase"]].itertuples(index=False, name=None)
    )
    visible = visible.loc[visible["_row"].isin(row_index)].copy()
    if visible.empty:
        return pd.DataFrame(columns=["start", "end", "y", "color"])

    visible["y"] = visible["_row"].map(row_index)
    visible["start"] = visible["bin_start_bp"].astype(float)
    visible["end"] = visible["bin_end_bp"].astype(float)
    visible["color"] = [mpl.colors.to_hex(cmap(value), keep_alpha=False) for value in visible["methylation_frequency"]]
    return visible[["start", "end", "y", "color"]]


def _add_vector_colorbar(ax: mpl.axes.Axes, cmap_name: str) -> None:

    cmap = mpl.colormaps[cmap_name]
    cax = ax.inset_axes([1.012, 0.18, COLORBAR_WIDTH, 0.64])
    edges = np.linspace(0.0, 1.0, 65)
    for lower, upper in zip(edges[:-1], edges[1:]):
        cax.axhspan(
            lower,
            upper,
            xmin=0.0,
            xmax=1.0,
            facecolor=cmap((lower + upper) / 2.0),
            edgecolor="none",
            linewidth=0,
        )
    cax.set_xlim(0, 1)
    cax.set_ylim(0, 1)
    cax.set_xticks([])
    cax.set_yticks([0, 0.5, 1])
    cax.yaxis.tick_right()
    cax.yaxis.set_label_position("right")
    cax.set_ylabel("5mC frequency", rotation=270, va="bottom", labelpad=8)
    for spine in cax.spines.values():
        spine.set_visible(False)


def _despine_panel(ax: mpl.axes.Axes) -> None:
    ax.spines["top"].set_visible(False)
    ax.spines["right"].set_visible(False)


def add_region_labels(ax: mpl.axes.Axes, config: MethylPlotter) -> None:
    xmin, xmax = config.resolved_xlim()
    if xmin < 0:
        ax.text(0.02, 0.96, "Subtelomere",transform=ax.transAxes, ha="left", va="top", color="0.3",)
    if xmax > 0:
        ax.text(0.98, 0.96, "Telomere",transform=ax.transAxes, ha="right", va="top", color="0.3",)


def figure_height(
    config: MethylPlotter,
    row_count: int,
    *,
    show_heatmap: bool,
) -> float:
    """Resolve the final page height for one methylation arm panel."""
    if config.height is not None:
        return config.height
    if show_heatmap:
        return max(4.2, 2.0 + 0.115 * row_count)
    return 2.3


def plot_arm(
    summary: pd.DataFrame,
    *,
    arm: str,
    output_stem: Path,
    config: MethylPlotter,
) -> list[Path]:

    arm_summary = summary.loc[summary["arm"] == arm].copy()
    if arm_summary.empty:
        return []

    xmin, xmax = config.resolved_xlim()
    rows = heatmap_rows(arm_summary, config) if config.show_heatmap else []
    show_heatmap = config.show_heatmap and bool(rows)
    phases = natsorted(arm_summary["chrom_phase"].unique().tolist())
    presentations = phase_presentations(phases)

    fig_height = figure_height(config, len(rows), show_heatmap=show_heatmap)

    with nature_rc():
        if show_heatmap:
            fig = plt.figure(figsize=(config.width, fig_height), layout="constrained")
            grid = fig.add_gridspec(nrows=2, ncols=1, height_ratios=(1.35, 5.0), hspace=0.04)
            curve_ax = fig.add_subplot(grid[0])
            heatmap_ax = fig.add_subplot(grid[1], sharex=curve_ax)
        else:
            fig, curve_ax = plt.subplots(figsize=(config.width, fig_height), layout="constrained")
            heatmap_ax = None

        for phase in phases:
            curve = aggregate_phase_curve(arm_summary, phase, config)
            x = curve["bin_center_bp"].to_numpy(dtype=float)
            mean = curve["methylation_frequency"].to_numpy(dtype=float)
            presentation = presentations[phase]
            curve_ax.plot(
                x, mean,
                color=presentation.curve_color,
                linewidth=STYLE.data_linewidth,
                label=presentation.label,
            )

            if config.show_sd:
                sd = curve["chromosome_sd"].to_numpy(dtype=float)
                finite = np.isfinite(mean) & np.isfinite(sd)
                curve_ax.fill_between(
                    x,
                    np.clip(mean - sd, 0.0, 1.0),
                    np.clip(mean + sd, 0.0, 1.0),
                    where=finite,
                    color=presentation.curve_color,
                    alpha=0.16,
                    linewidth=0,
                    interpolate=False,
                    zorder=1,
                )

        if xmin < 0 < xmax:
            curve_ax.axvline(0, color="0.25", linestyle="--", linewidth=0.75, zorder=0)
        curve_ax.set_xlim(xmin, xmax)
        curve_ax.set_ylim(0, 1)
        curve_ax.set_ylabel("5mC frequency")
        _despine_panel(curve_ax)
        add_region_labels(curve_ax, config)
        if any(presentation.label for presentation in presentations.values()):
            curve_ax.legend(loc="upper left",bbox_to_anchor=(1.01, 1.0),ncol=1,handlelength=2.2,borderaxespad=0.0,)

        if heatmap_ax is not None:
            curve_ax.tick_params(axis="x", which="both", bottom=False, labelbottom=False)
            heatmap_ax.set_facecolor(MISSING_COLOR)
            cells = _heatmap_cells(arm_summary, rows, config)
            plot_intervals(
                heatmap_ax,
                cells,
                start="start",
                end="end",
                y="y",
                color="color",
                thickness=1.0,
                rasterized=False,
                zorder=1,
            )

            if len(rows) > 1:
                heatmap_ax.hlines(
                    np.arange(len(rows) - 1, dtype=float) + 0.5,
                    xmin,
                    xmax,
                    colors="white",
                    linewidth=0.25,
                    zorder=2,
                )
            if xmin < 0 < xmax:
                heatmap_ax.axvline(0, color="0.2", linestyle="--", linewidth=0.75, zorder=3)

            heatmap_ax.set_xlim(xmin, xmax)
            heatmap_ax.set_ylim(len(rows) - 0.5, -0.5)
            heatmap_ax.set_xlabel("Distance to telomere-subtelomere boundary (kb)")
            heatmap_ax.set_yticks(np.arange(len(rows)))
            heatmap_ax.set_yticklabels(
                [
                    methyl_row_label(chrom, phase, arm)
                    for chrom, phase in rows
                ],
                fontsize=STYLE.row_label_size,
            )
            for label, (_chrom, phase) in zip(heatmap_ax.get_yticklabels(), rows):
                label.set_color(presentations[phase].row_label_color)

            format_kb_axis(heatmap_ax)
            _despine_panel(heatmap_ax)
            _add_vector_colorbar(heatmap_ax, config.heatmap_cmap)
        else:
            curve_ax.set_xlabel("Distance to telomere-subtelomere boundary (kb)")
            format_kb_axis(curve_ax)

        try:
            path = save_figure(fig, output_stem, output_format=config.output_format)
        finally:
            plt.close(fig)

    return [path]


def _read_methyl_table(path: Path) -> pd.DataFrame:
    if not path.exists():
        raise ValueError(f"Methylation plot-data file not found: {path}")
    return pd.read_csv(
        path,
        sep="\t",
        dtype={"read_id": str, "chrom": str, "chrom_phase": str, "arm": str},
        keep_default_na=False,
    )


def plot_arm_methyl(
    mods_file: str | Path,
    outdir: str | Path = ".",
    prefix: str = "sample",
    config: MethylPlotter | None = None,
    logger=None,
) -> list[Path]:

    logger = logger or log_utils.get_logger()
    config = config or MethylPlotter()
    config.validate()

    output_dir = Path(outdir)
    output_dir.mkdir(parents=True, exist_ok=True)
    mods_path = Path(mods_file)
    logger.info(f"Loading sample data: {prefix}")
    table = _read_methyl_table(mods_path)
    summary = summarise_sites(table, config, source=mods_path)

    written: list[Path] = []
    arms = sorted(summary["arm"].unique().tolist(), key=arm_sort_key)
    for arm in arms:
        logger.info(f"Plotting 5mC boundary profile: {prefix}, {arm}-arm")
        figure_paths = plot_arm(
            summary,
            arm=arm,
            output_stem=output_dir / f"{prefix}.{arm}-arm.methyl_plot",
            config=config,
        )
        written.extend(figure_paths)

    return written
