# Copyright (C) 2025 Huihui Li <hhui.li@outlook.com>. Licensed under GNU GPL v3.0.
from pathlib import Path
from dataclasses import dataclass
from typing import Optional, Union

import numpy as np
import pandas as pd
from natsort import natsorted

from . import telox_utils
from . import logger as log_utils
from .viz_core import (
    STYLE,
    chrom_phase_key,
    chrom_phase_label,
    display_chrom,
    nature_rc,
    normalize_chrom_phase,
    require_schema,
    resolve_requested_values,
    save_figure,
)
import matplotlib.pyplot as plt
import matplotlib.ticker as ticker
import matplotlib.patches as patches
from matplotlib.lines import Line2D

LENGTH_RASTERIZE_THRESHOLD = 20_000


@dataclass
class LengthPlotter:
    chroms: Optional[str] = None
    arm_colors: Optional[str] = None
    rasterize: Optional[bool] = None
    output_format: str = "pdf"
    width: float = 8.0
    height: float = 3.0

    @classmethod
    def from_params(cls, **kwargs):
        return cls(**{k: v for k, v in kwargs.items() if k in cls.__dataclass_fields__ and v is not None})


@dataclass(frozen=True, slots=True)
class _LengthData:
    plot_levels: tuple[object, ...]
    arm_labels: tuple[str, ...]
    row_order: tuple[str | None, ...]
    grouped_lengths: dict[tuple[object, ...], np.ndarray]
    row_medians: dict[str | None, float]
    n_reads: int


def _rasterize_points(config: LengthPlotter, n_reads: int) -> bool:
    if config.rasterize is not None:
        return config.rasterize
    return n_reads > LENGTH_RASTERIZE_THRESHOLD


def _prepare_length_data(tel_length: pd.DataFrame, config: LengthPlotter) -> _LengthData | None:
    require_schema(tel_length)
    if "tel_length" not in tel_length.columns:
        raise ValueError("Missing required column: tel_length.")
    if tel_length.empty:
        return None

    tidy = tel_length.copy()
    tidy["tel_length"] = pd.to_numeric(tidy["tel_length"], errors="raise")

    mitochondrial = (
        tidy["chrom"]
        .map(display_chrom)
        .str.casefold()
        .isin({"m", "mt", "mito"})
    )
    tidy = tidy.loc[~mitochondrial].copy()
    if tidy.empty:
        return None

    tidy["_chrom_phase"] = tidy["chrom_phase"].map(normalize_chrom_phase)
    is_haplotype_data = tidy["_chrom_phase"].map(chrom_phase_key).ne("").any()

    observed = tuple(dict.fromkeys(tidy["chrom"].dropna().tolist()))
    selected_text = resolve_requested_values(
        config.chroms,
        observed,
        label="chromosome",
    )
    observed_by_text = {str(value): value for value in observed}
    selected = tuple(observed_by_text[value] for value in selected_text)
    if not selected:
        return None

    tidy = tidy.loc[tidy["chrom"].isin(selected)].copy()
    if tidy.empty:
        return None

    if config.chroms:
        initial_chrom_levels = selected
    elif is_haplotype_data:
        initial_chrom_levels = tuple(natsorted(selected))
    else:
        sorter = telox_utils.create_chromosome_sorter(selected)
        initial_chrom_levels = tuple(sorted(selected, key=sorter))

    tidy["chrom"] = pd.Categorical(
        tidy["chrom"], categories=initial_chrom_levels, ordered=True
    )
    tidy.dropna(subset=["chrom"], inplace=True)
    if tidy.empty:
        return None

    observed_arms = {str(arm) for arm in tidy["arm"].dropna().unique()}
    if not observed_arms:
        return None
    if observed_arms <= {"p", "q"}:
        arm_labels = tuple(arm for arm in ("p", "q") if arm in observed_arms)
    elif observed_arms <= {"L", "R"}:
        arm_labels = tuple(arm for arm in ("L", "R") if arm in observed_arms)
    else:
        observed_labels = ", ".join(natsorted(observed_arms))
        raise ValueError(
            "Length plotting requires one arm convention (p/q or L/R); "
            f"observed mixed or unsupported arm values: {observed_labels}."
        )

    row_var: str | None = None
    if is_haplotype_data:
        tidy["_phase_group"] = tidy["_chrom_phase"].map(chrom_phase_key)
        if tidy["_phase_group"].nunique() > 1:
            row_var = "_phase_group"
        plot_levels = tuple(natsorted(tidy["chrom"].unique()))
        tidy["chrom"] = pd.Categorical(tidy["chrom"], categories=plot_levels, ordered=True)
    else:
        plot_levels = initial_chrom_levels

    if row_var:
        row_order: tuple[str | None, ...] = tuple(natsorted(tidy[row_var].dropna().unique()))
        group_columns = [row_var, "chrom", "arm"]
        row_medians = {
            str(phase): float(median)
            for phase, median in tidy.groupby(row_var, observed=True, sort=False)["tel_length"].median().items()
        }
    else:
        row_order = (None,)
        group_columns = ["chrom", "arm"]
        row_medians = {None: float(tidy["tel_length"].median())}

    grouped_lengths = {
        key if isinstance(key, tuple) else (key,): group["tel_length"]
        .dropna()
        .to_numpy()
        for key, group in tidy.groupby(
            group_columns,
            observed=True,
            sort=False,
        )
    }
    n_reads = int(tidy.loc[tidy["arm"].isin(arm_labels), "tel_length"].notna().sum())
    return _LengthData(
        plot_levels=plot_levels,
        arm_labels=arm_labels,
        row_order=row_order,
        grouped_lengths=grouped_lengths,
        row_medians=row_medians,
        n_reads=n_reads,
    )


def _draw_length_panel(
    ax,
    *,
    row_name: str | None,
    data: _LengthData,
    palette: dict[str, str],
    rng: np.random.Generator,
    rasterized: bool,
) -> None:
    for index in range(len(data.plot_levels) - 1):
        ax.axvline(
            x=index + 0.5,
            linestyle="--",
            color="black",
            alpha=0.2,
            linewidth=STYLE.separator_linewidth,
        )

    median_val = data.row_medians.get(row_name, np.nan)
    if pd.notna(median_val):
        ax.text(
            0.02,
            0.95,
            f"Genome-wide median: {median_val:.0f} bp",
            transform=ax.transAxes,
            color="black",
            fontsize=8,
            ha="left",
            va="top",
        )

    for chrom_index, chrom in enumerate(data.plot_levels):
        for arm_index, arm in enumerate(data.arm_labels):
            key = (
                (row_name, chrom, arm)
                if row_name is not None
                else (chrom, arm)
            )
            arm_data = data.grouped_lengths.get(key)
            if arm_data is None or len(arm_data) == 0:
                continue

            color = palette[arm]
            position = chrom_index + (-0.2 if arm_index == 0 else 0.2)
            n_points = len(arm_data)
            dynamic_alpha = max(0.1, min(0.6, 300.0 / n_points))
            jitter = rng.uniform(-0.05, 0.05, size=n_points)
            ax.scatter(
                position + jitter,
                arm_data,
                color=color,
                s=1.5,
                alpha=dynamic_alpha,
                linewidths=0,
                zorder=1,
                rasterized=rasterized,
            )
            ax.boxplot(
                arm_data,
                positions=[position],
                widths=0.3,
                showfliers=False,
                patch_artist=True,
                boxprops=dict(facecolor="none", edgecolor=color, linewidth=0.7, zorder=2),
                whiskerprops=dict(color=color, linewidth=0.7, zorder=2),
                capprops=dict(color=color, linewidth=0.7, zorder=2),
                medianprops=dict(color=color, linewidth=0.7, zorder=3),
            )


def plot_length(
        length_input: Union[str, Path],
        outdir: Union[str, Path],
        prefix: str,
        config: LengthPlotter,
        logger = None
    ):
    logger = logger or log_utils.get_logger()
    outdir = Path(outdir)

    tel_length = pd.read_csv(length_input, sep='\t')
    data = _prepare_length_data(tel_length, config)
    if data is None:
        return

    default_colors = ['#45732b','#f97644']
    user_colors = [c.strip() for c in config.arm_colors.split(",")] if config.arm_colors else []
    palette = dict(zip(
        data.arm_labels,
        user_colors if len(user_colors) >= len(data.arm_labels) else default_colors,
    ))

    n_rows = len(data.row_order)
    n_cols = len(data.plot_levels)
    rasterized = _rasterize_points(config, data.n_reads)

    with nature_rc():
        fig, axes = plt.subplots(
            nrows=n_rows, ncols=1,
            figsize=(config.width, config.height),
            sharex=True, sharey=True, squeeze=False,
            constrained_layout=True
        )
        fig.suppressComposite = rasterized

        rng = np.random.default_rng(42)

        for row_index, row_name in enumerate(data.row_order):
            ax = axes[row_index, 0]
            _draw_length_panel(
                ax,
                row_name=row_name,
                data=data,
                palette=palette,
                rng=rng,
                rasterized=rasterized,
            )

            ax.yaxis.set_major_formatter(ticker.FuncFormatter(lambda x, p: f"{x/1000:g}"))
            for spine in ax.spines.values():
                spine.set_visible(True)
                spine.set_edgecolor('black')

            ax.tick_params(axis='x', length=3)
            ax.tick_params(axis='y', length=3)

            if row_name:
                rect = patches.Rectangle(
                    (1.0, 0), 0.03, 1.0,
                    transform=ax.transAxes,
                    clip_on=False,
                    facecolor='#D9D9D9',
                    edgecolor='black',
                    linewidth=STYLE.axis_linewidth,
                    zorder=10
                )
                ax.add_patch(rect)

                ax.text(
                    1.015, 0.5,
                    chrom_phase_label(row_name),
                    transform=ax.transAxes,
                    rotation=-90,
                    va='center',
                    ha='center',
                    fontsize=8,
                    color='black',
                    clip_on=False,
                    zorder=100
                )

        bottom_ax = axes[-1, 0]
        bottom_ax.set_xticks(range(n_cols))
        bottom_ax.set_xticklabels([display_chrom(chrom) for chrom in data.plot_levels])
        bottom_ax.set_xlabel("Chromosome", fontsize=8)

        if n_rows % 2 == 1:
            target_ax = axes[n_rows // 2, 0]
            y_pos = 0.5
        else:
            target_ax = axes[n_rows // 2, 0]
            y_pos = 1.0

        target_ax.set_ylabel("Telomere length (kb)", y=y_pos, fontsize=8)

        legend_elements = [
            Line2D([0], [0], color=palette[arm], lw=1.5, label=f"{arm}-arm")
            for arm in data.arm_labels
        ]
        axes[0, 0].legend(
            handles=legend_elements,
            loc="upper right",
            bbox_to_anchor=(1.0, 1.0),
            title=None,
            frameon=True,
            facecolor='white',
            framealpha=0.5,
            ncol=2,
            columnspacing=1.0,
            handletextpad=0.4,
            handlelength=1.0,
            fontsize=8
        )

        output = save_figure(fig, outdir / f"{prefix}.length_plot.pdf", output_format=config.output_format,)
        plt.close(fig)
    logger.info(f"Plot saved to: {output}")
