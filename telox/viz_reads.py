# Copyright (C) 2025 Huihui Li <hhui.li@outlook.com>. Licensed under GNU GPL v3.0.
from __future__ import annotations

from collections import Counter
from dataclasses import dataclass
from pathlib import Path
from typing import Mapping, Sequence

import matplotlib.pyplot as plt
import matplotlib.colors as mcolors
import numpy as np
import pandas as pd
from matplotlib.collections import LineCollection
from matplotlib.lines import Line2D
from natsort import natsorted

from . import logger as log_utils
from .viz_core import (
    STYLE,
    arm_group_key,
    chrom_phase_key,
    chrom_phase_label,
    display_chrom,
    format_kb_axis,
    load_motif_color_overrides,
    motif_color_map,
    nature_rc,
    normalize_figure_format,
    parse_motif_blocks,
    parse_sample_sheet,
    require_schema,
    resolve_requested_values,
    row_major,
    safe_name,
    save_figure,
    split_csv,
    write_motif_color_map,
)


TVR_REQUIRED = {"read_id", "chrom", "chrom_phase", "arm", "tvr_hap", "tel_length", "tvrs",}
MODS_REQUIRED = {"read_id", "chrom", "chrom_phase", "arm", "read_length", "valid_sites", "mods_m"}
CONS_REQUIRED = {
    "chrom", "chrom_phase", "arm", "tvr_hap", "read_support", "hap_length", "motif_blocks",
}

READ_BG = "#E7EAED"
OTHER_COLOR = "#CCCCCC"
BOUNDARY = "#000000"
METHYL_COLOR = "#D62728"
UNMETHYL_COLOR = "#C8CDD3"
MOD_READ_BG = "#E7EAED"


@dataclass(slots=True)
class TVRReadPlotter:
    canonical_motif: str = "TTAGGG"
    chroms: str | Sequence[str] | None = "all"
    pool_chroms: bool = False
    xlim: str | Sequence[int | float] = (0, 2000)
    top_tvr: int = 10
    tvr_colors: Path | None = None
    read_thickness: float = 0.8
    consensus_thickness: float = 1.0
    show_outliers: bool = False
    show_unmethylated_c: bool = False
    legend_ncol: int = 4
    sample_label: str | None = None
    width: float = STYLE.inches(STYLE.double_column_mm)
    height: float | None = None
    output_format: str = "pdf"

    @classmethod
    def from_params(cls, **kwargs):
        return cls(**{k: v for k, v in kwargs.items() if k in cls.__dataclass_fields__ and v is not None})


@dataclass(frozen=True, slots=True)
class Panel:
    arm_key: str
    arm: str
    pooled: bool
    stem: str
    title: str
    tvr: pd.DataFrame
    mods: pd.DataFrame | None
    consensus: pd.DataFrame | None


@dataclass(frozen=True, slots=True)
class HaplotypeGroup:
    label: str
    tvr: pd.DataFrame
    consensus: pd.DataFrame | None
    mods: pd.DataFrame | None


@dataclass(slots=True)
class PreparedPanel:
    tvr: pd.DataFrame
    mods: pd.DataFrame | None
    cons_bg: pd.DataFrame
    cons_fg: pd.DataFrame
    labels: list[tuple[str, float]]
    total_y: float
    cons_visible: set[str]


def _xlim(value: str | Sequence[int | float]) -> tuple[float, float] | None:
    if isinstance(value, str) and value.strip().casefold() == "none":
        return None
    parts = split_csv(value) if isinstance(value, str) else tuple(map(str, value))
    if len(parts) != 2:
        raise ValueError("xlim must contain exactly MIN_BP,MAX_BP.")
    xmin, xmax = map(float, parts)
    if not np.isfinite([xmin, xmax]).all() or xmin >= xmax:
        raise ValueError("xlim requires two finite values with MIN_BP < MAX_BP.")
    return xmin, xmax


def _panel_xlim(panel: Panel) -> tuple[float, float]:
    xmin = 0.0
    xmax = 0.0

    span_starts = panel.tvr["span_start"].to_numpy(dtype=float)
    span_ends = panel.tvr["span_end"].to_numpy(dtype=float)
    finite_starts = span_starts[np.isfinite(span_starts)]
    finite_ends = span_ends[np.isfinite(span_ends)]
    if finite_ends.size:
        xmax = max(xmax, float(finite_ends.max()))

    if panel.mods is not None:
        xmin = -xmax
    elif finite_starts.size:
        xmin = min(xmin, float(finite_starts.min()))

    if panel.consensus is not None and not panel.consensus.empty:
        values = panel.consensus["hap_length"].to_numpy(dtype=float)
        finite = values[np.isfinite(values)]
        if finite.size:
            xmax = max(xmax, float(finite.max()))

    if panel.mods is not None and not panel.mods.empty:
        keys = ["read_id", "chrom", "_phase", "_arm"]
        identities = panel.tvr[keys].drop_duplicates()
        matched_mods = panel.mods.merge(
            identities,
            on=keys,
            how="inner",
            validate="many_to_one",
        )
        for column in ("valid_pos", "methyl_pos"):
            for positions in matched_mods[column]:
                finite = np.asarray(positions, dtype=float)
                finite = finite[np.isfinite(finite)]
                if finite.size:
                    xmax = max(xmax, float(finite.max()))

    if np.isclose(xmin, xmax):
        xmax = xmin + 1.0
    return xmin, xmax


def _positions(value: object) -> np.ndarray:
    if pd.isna(value) or str(value).strip() in {"", "."}:
        return np.empty(0, dtype=float)
    return np.fromstring(str(value), sep=",", dtype=float)


def _tvr_calls(value: object) -> tuple[tuple[float, float, str], ...]:
    if pd.isna(value) or str(value).strip() in {"", "."}:
        return ()
    calls = []
    for token in str(value).split(";"):
        fields = token.split(",", 2)
        if len(fields) != 3:
            raise ValueError(f"Invalid TVR call: {token!r}")
        calls.append((float(fields[0]), float(fields[1]), fields[2]))
    return tuple(calls)


def _phase_group_key(value: object) -> str:
    return chrom_phase_key(value).casefold()


def _arm_label(values: pd.Series, key: str) -> str:
    observed = {str(v).lower() for v in values}
    if key == "left":
        return "p" if "p" in observed else "L"
    if key == "right":
        return "q" if "q" in observed else "R"
    return key.removeprefix("other:")


def _group_label(chrom: object, phase: object, arm: str) -> str:
    chrom_arm = f"{display_chrom(chrom)}{arm}"
    parent = chrom_phase_label(phase)
    return f"{parent}-{chrom_arm}" if parent else chrom_arm


def _plot_title(panel: Panel, title_prefix: str) -> str:
    phase = chrom_phase_label(panel.tvr["_phase"].iloc[0])
    phase_suffix = f" ({phase})" if phase else ""
    chromosome_label = (
        f"{panel.arm}-arm"
        if panel.pooled
        else f"{display_chrom(panel.tvr['chrom'].iloc[0])}{panel.arm}"
    )
    return f"{title_prefix} chromosome {chromosome_label}{phase_suffix}"


def _load_tvr(path: Path, sample: str, *, show_outliers: bool = False) -> pd.DataFrame:
    df = pd.read_csv(path, sep="\t", na_values=".")
    require_schema(df, TVR_REQUIRED, source=path)
    if not show_outliers:
        df = df[df["tvr_hap"].fillna("").astype(str).str.casefold().ne("outlier")].copy()
    df["tel_length"] = pd.to_numeric(df["tel_length"], errors="raise")
    cache = {raw: _tvr_calls(raw) for raw in pd.unique(df["tvrs"].fillna("."))}
    df["tvr_calls"] = df["tvrs"].fillna(".").map(cache)
    df["sample_id"] = sample
    df["_phase"] = df["chrom_phase"].map(_phase_group_key)
    df["_arm"] = df["arm"].map(arm_group_key)
    df["span_start"] = 0.0
    df["span_end"] = df["tel_length"].astype(float)
    return df


def _load_mods(path: Path, *, include_unmethylated: bool) -> pd.DataFrame:
    df = pd.read_csv(path, sep="\t", na_values=".")
    require_schema(df, MODS_REQUIRED, source=path)
    df["read_length"] = pd.to_numeric(df["read_length"], errors="raise")
    df["_phase"] = df["chrom_phase"].map(_phase_group_key)
    df["_arm"] = df["arm"].map(arm_group_key)
    for source, target in (("valid_sites", "valid_pos"), ("mods_m", "methyl_pos")):
        raw = df[source].fillna(".")
        cache = {value: _positions(value) for value in pd.unique(raw)}
        df[target] = raw.map(cache)
    if include_unmethylated:
        df["unmethyl_pos"] = [
            valid if len(methyl) == 0 else valid[~np.isin(valid, methyl)]
            for valid, methyl in zip(df["valid_pos"], df["methyl_pos"], strict=True)
        ]
    return df


def _load_consensus(path: Path) -> pd.DataFrame:
    df = pd.read_csv(path, sep="\t", na_values=".")
    require_schema(df, CONS_REQUIRED, source=path)
    df["read_support"] = pd.to_numeric(df["read_support"], errors="raise")
    df["hap_length"] = pd.to_numeric(df["hap_length"], errors="raise")
    raw = df["motif_blocks"].astype(str)
    cache = {value: parse_motif_blocks(value) for value in pd.unique(raw)}
    df["tokens"] = raw.map(cache)
    df["_phase"] = df["chrom_phase"].map(_phase_group_key)
    df["_arm"] = df["arm"].map(arm_group_key)
    return df


def _add_read_lengths(tvr: pd.DataFrame, mods: pd.DataFrame | None) -> pd.DataFrame:
    if mods is None:
        return tvr
    keys = ["read_id", "chrom", "_phase", "_arm"]
    lengths = mods[keys + ["read_length"]].drop_duplicates(keys)
    tvr = tvr.merge(lengths, on=keys, how="left", validate="many_to_one")
    tvr["span_start"] = np.where(
        tvr["read_length"].notna(),
        tvr["tel_length"].astype(float) - tvr["read_length"].astype(float),
        0.0,
    )
    return tvr


def _top_tvrs(
    df: pd.DataFrame,
    samples: Sequence[str],
    top_n: int,
    xlim: tuple[float, float] | None,
    color_file: Path | None,
    canonical_motif: str,
) -> tuple[dict[tuple[str, str], tuple[str, ...]], dict[str, str]]:
    xmin, xmax = xlim if xlim is not None else (-np.inf, np.inf)
    by_sample_arm: dict[tuple[str, str], tuple[str, ...]] = {}
    candidates: list[str] = []
    seen: set[str] = set()
    global_counts: Counter[str] = Counter()
    arms = natsorted(df["_arm"].unique())
    sample_groups = {
        str(sample): rows
        for sample, rows in df.groupby("sample_id", sort=False, observed=True)
    }

    for sample in samples:
        sample_df = sample_groups.get(str(sample))
        arm_groups = (
            {}
            if sample_df is None
            else {
                str(arm): rows
                for arm, rows in sample_df.groupby("_arm", sort=False, observed=True)
            }
        )
        for arm in arms:
            counts: Counter[str] = Counter()
            arm_rows = arm_groups.get(str(arm))
            if arm_rows is not None:
                for calls in arm_rows["tvr_calls"]:
                    counts.update(
                        motif
                        for start, end, motif in calls
                        if end > xmin and start < xmax
                    )
            global_counts.update(counts)
            top = tuple(m for m, _ in counts.most_common(top_n))
            by_sample_arm[(sample, arm)] = top
            for motif in top:
                if motif not in seen:
                    seen.add(motif)
                    candidates.append(motif)

    abundance = Counter({motif: global_counts[motif] for motif in candidates})
    colors, _ = motif_color_map(
        abundance,
        len(abundance),
        load_motif_color_overrides(color_file),
        priority_motifs=(canonical_motif,),
    )
    return by_sample_arm, colors


def _panel_motifs(
    by_sample_arm: Mapping[tuple[str, str], tuple[str, ...]],
    samples: Sequence[str],
    arm: str,
) -> tuple[str, ...]:
    seen: set[str] = set()
    result = []
    for sample in samples:
        for motif in by_sample_arm.get((sample, arm), ()):
            if motif not in seen:
                seen.add(motif)
                result.append(motif)
    return tuple(result)


def _select_chroms(df: pd.DataFrame, value: str | Sequence[str] | None) -> tuple[str, ...]:
    observed = natsorted(df["chrom"].astype(str).unique())
    return resolve_requested_values(value, observed, label="chromosome")


def _panels(
    tvr: pd.DataFrame,
    mods: pd.DataFrame | None,
    consensus: pd.DataFrame | None,
    pool_chroms: bool,
) -> list[Panel]:
    group_columns = ["_phase", "_arm"] if pool_chroms else ["chrom", "_phase", "_arm"]
    selected_chroms = set(tvr["chrom"])

    def related_groups(
        frame: pd.DataFrame | None,
    ) -> dict[tuple[object, ...], pd.DataFrame]:
        if frame is None:
            return {}
        selected = frame.loc[frame["chrom"].isin(selected_chroms)] if pool_chroms else frame
        return {
            key if isinstance(key, tuple) else (key,): rows.copy()
            for key, rows in selected.groupby(
                group_columns,
                sort=False,
                observed=True,
            )
        }

    mods_groups = related_groups(mods)
    consensus_groups = related_groups(consensus)
    result: list[Panel] = []
    if pool_chroms:
        groups = list(tvr.groupby(["_phase", "_arm"], sort=False, observed=True))
        groups = natsorted(groups, key=lambda item: item[0])
        for (phase, arm_key), rows in groups:
            arm = _arm_label(rows["arm"], arm_key)
            phase_tag = safe_name(phase) if phase else "unphased"
            phase_title = f" ({chrom_phase_label(phase)})" if phase else ""
            key = (phase, arm_key)
            result.append(Panel(
                arm_key, arm, True, f"all_chroms.{phase_tag}.{arm}-arm",
                f"pooled chromosomes {arm}-arm{phase_title}", rows.copy(),
                mods_groups.get(key), consensus_groups.get(key),
            ))
        return result

    groups = list(tvr.groupby(["chrom", "_phase", "_arm"], sort=False, observed=True))
    groups = natsorted(groups, key=lambda item: item[0])
    for (_chrom, phase, arm_key), rows in groups:
        chrom = str(rows["chrom"].iloc[0])
        arm = _arm_label(rows["arm"], arm_key)
        phase_suffix = f"_{safe_name(phase)}" if phase else ""
        key = (_chrom, phase, arm_key)
        result.append(Panel(
            arm_key, arm, False, f"{safe_name(chrom)}{phase_suffix}.{arm}-arm",
            _group_label(chrom, phase, arm), rows.copy(),
            mods_groups.get(key), consensus_groups.get(key),
        ))
    return result


def _haplotype_groups(panel: Panel, multi_sample: bool) -> list[HaplotypeGroup]:
    df = panel.tvr
    if multi_sample:
        raw_groups = list(df.groupby(["sample_id", "tvr_hap"], sort=False, observed=True))
        label = lambda key, rows: (
            f"{key[0]} {display_chrom(rows['chrom'].iloc[0])}{panel.arm} {key[1]}"
        )
    elif panel.pooled:
        raw_groups = list(df.groupby(["chrom", "_phase", "tvr_hap"], sort=False, observed=True))
        label = lambda key, rows: (
            f"{display_chrom(rows['chrom'].iloc[0])}{panel.arm} {key[2]}"
        )
    else:
        raw_groups = list(df.groupby("tvr_hap", sort=False, observed=True))
        label = lambda key, _rows: str(key)

    def match(
        source: pd.DataFrame | None,
        reads: pd.DataFrame,
        keys: list[str],
    ) -> pd.DataFrame | None:
        if source is None:
            return None
        identity = reads[keys].drop_duplicates()
        return source.merge(identity, on=keys, how="inner", validate="many_to_one")

    result: list[HaplotypeGroup] = []
    for key, reads in natsorted(raw_groups, key=lambda item: item[0]):
        consensus = match(panel.consensus, reads, ["chrom", "_phase", "_arm", "tvr_hap"])
        mods = match(panel.mods, reads, ["read_id", "chrom", "_phase", "_arm"])
        result.append(HaplotypeGroup(label(key, reads), reads.copy(), consensus, mods))
    return result


def _assign_read_y(df: pd.DataFrame) -> tuple[pd.DataFrame, float, float]:
    df = df.sort_values(["tel_length", "read_id"], ascending=[False, True], kind="stable").copy()
    df["y"] = np.arange(len(df), dtype=float)
    return df, float(df["y"].mean()), float(len(df))


def _consensus_intervals(
    df: pd.DataFrame | None,
    selected: Sequence[str],
    colors: Mapping[str, str],
    canonical: str,
    xlim: tuple[float, float],
) -> tuple[pd.DataFrame, pd.DataFrame, float, set[str]]:
    cols = ["start", "end", "y", "color"]
    if df is None or df.empty:
        empty = pd.DataFrame(columns=cols)
        return empty, empty.copy(), 0.0, set()

    df = df.sort_values(
        ["chrom", "_phase", "read_support", "hap_length"],
        ascending=[True, True, False, False], kind="stable",
    )
    xmin, xmax = xlim
    selected = set(selected)
    bg, fg, visible = [], [], set()

    for y, row in enumerate(df.itertuples(index=False)):
        start, end = max(0.0, xmin), min(float(row.hap_length), xmax)
        if end > start:
            bg.append((start, end, float(y), READ_BG))
        pos = 0.0
        for motif, count in row.tokens:
            end = pos + len(motif) * int(count)
            a, b = max(pos, xmin), min(end, xmax)
            if motif != canonical and b > a:
                visible.add(motif)
                fg.append((a, b, float(y), colors.get(motif, OTHER_COLOR) if motif in selected else OTHER_COLOR))
            pos = end
            if pos >= xmax:
                break

    return (
        pd.DataFrame(bg, columns=cols), pd.DataFrame(fg, columns=cols),
        float(len(df)), visible,
    )


def _align_mods(mods: pd.DataFrame | None, tvr: pd.DataFrame) -> pd.DataFrame | None:
    if mods is None or mods.empty:
        return None
    keys = ["read_id", "chrom", "_phase", "_arm"]
    placement = tvr[keys + ["y", "span_start", "span_end"]].drop_duplicates(keys)
    aligned = mods.merge(placement, on=keys, how="inner", validate="many_to_one")
    return aligned if not aligned.empty else None


def _prepare_panel(
    panel: Panel,
    motifs: Sequence[str],
    colors: Mapping[str, str],
    config: TVRReadPlotter,
    xlim: tuple[float, float],
    multi_sample: bool,
) -> PreparedPanel:
    groups = _haplotype_groups(panel, multi_sample)
    show_group = len(groups) > 1 or panel.pooled or multi_sample
    clearance = max(0.0, (config.consensus_thickness - 1.0) / 2.0)
    read_parts: list[pd.DataFrame] = []
    mod_parts: list[pd.DataFrame] = []
    cons_bg_parts: list[pd.DataFrame] = []
    cons_fg_parts: list[pd.DataFrame] = []
    labels: list[tuple[str, float]] = []
    visible: set[str] = set()
    cursor = 0.0

    def track_label(track: str, group: HaplotypeGroup) -> str:
        if not show_group:
            return track
        return "-".join((track, *str(group.label).split()))

    for index, group in enumerate(groups):
        if index:
            cursor += 1.0

        reads, read_center, read_height = _assign_read_y(group.tvr)
        mods = _align_mods(group.mods, reads)
        reads["y"] = reads["y"].to_numpy(float) + cursor
        read_parts.append(reads)
        labels.append((track_label("TVR", group), read_center + cursor))
        cursor += read_height

        cons_bg, cons_fg, cons_height, cons_motifs = _consensus_intervals(
            group.consensus, motifs, colors, config.canonical_motif, xlim
        )
        if cons_height:
            cons_offset = cursor + 0.75 + clearance
            cons_bg["y"] = cons_bg["y"].to_numpy(float) + cons_offset
            cons_fg["y"] = cons_fg["y"].to_numpy(float) + cons_offset
            cons_bg_parts.append(cons_bg)
            cons_fg_parts.append(cons_fg)
            labels.append((track_label("Consensus", group), cons_offset + (cons_height - 1.0) / 2.0))
            cursor = cons_offset + cons_height + clearance
            visible.update(cons_motifs)

        if mods is not None:
            mods_offset = cursor + 0.75
            mods["y"] = mods["y"].to_numpy(float) + mods_offset
            mod_parts.append(mods)
            labels.append((track_label("5mC", group), read_center + mods_offset))
            cursor = mods_offset + read_height

    interval_cols = ["start", "end", "y", "color"]
    empty_intervals = pd.DataFrame(columns=interval_cols)
    return PreparedPanel(
        tvr=pd.concat(read_parts, ignore_index=True),
        mods=pd.concat(mod_parts, ignore_index=True) if mod_parts else None,
        cons_bg=pd.concat(cons_bg_parts, ignore_index=True) if cons_bg_parts else empty_intervals,
        cons_fg=pd.concat(cons_fg_parts, ignore_index=True) if cons_fg_parts else empty_intervals.copy(),
        labels=labels,
        total_y=cursor,
        cons_visible=visible,
    )


def _linewidth(ax: plt.Axes, thickness: float) -> float:
    y = ax.transData.transform([(0, 0), (0, thickness)])[:, 1]
    return abs(y[1] - y[0]) * 72.0 / ax.figure.dpi


def _xaxis(ax: plt.Axes, xlim: tuple[float, float]) -> None:
    ax.set_xlim(*xlim)
    format_kb_axis(ax)


def _decorate(
    ax: plt.Axes,
    xlim: tuple[float, float],
    total_y: float,
    labels: Sequence[tuple[str, float]],
) -> None:
    ax.set(xlim=xlim, ylim=(max(0.5, total_y - 0.5), -0.5))
    ax.set_yticks([y for _, y in labels], [label for label, _ in labels])
    # Hidden ticks still contribute their nominal length to the tick-label
    # transform.  Collapse that unused geometry so labels retain only
    # Matplotlib's normal typographic padding from the y-axis spine.
    ax.tick_params(axis="y", which="both", left=False, length=0)
    ax.spines[["top", "right"]].set_visible(False)
    # TVR segments can start exactly at the left boundary, and a thick
    # consensus can reach the bottom boundary.  Keep the visible frame above
    # those data artists so neither axis line is hidden.
    for name in ("left", "bottom"):
        ax.spines[name].set_zorder(5)
    _xaxis(ax, xlim)
    ax.set_xlabel("Distance to telomere-subtelomere boundary (kb)", labelpad=3)


def _line_collection(
    ax: plt.Axes,
    starts: np.ndarray,
    ends: np.ndarray,
    ys: np.ndarray,
    colors,
    linewidth: float,
    capstyle: str,
    zorder: int,
) -> None:
    if len(starts) == 0:
        return
    segments = np.stack((np.column_stack((starts, ys)), np.column_stack((ends, ys))), axis=1)
    ax.add_collection(LineCollection(
        segments, colors=colors, linewidths=linewidth, capstyle=capstyle,
        rasterized=False, zorder=zorder,
    ))


def _plot_vertical_calls(
    ax: plt.Axes,
    frame: pd.DataFrame,
    *,
    x: str,
    y: str,
    probability: str | None = None,
    color: str = METHYL_COLOR,
    height: float = 0.7,
    linewidth: float = 0.5,
    rasterized: bool = False,
    zorder: int = 3,
) -> LineCollection | None:
    if frame.empty:
        return None
    xs, ys = frame[x].to_numpy(float), frame[y].to_numpy(float)
    segments = np.stack((
        np.column_stack((xs, ys - height / 2)),
        np.column_stack((xs, ys + height / 2)),
    ), axis=1)
    if color in frame.columns:
        rgba = np.asarray([mcolors.to_rgba(value) for value in frame[color]], dtype=float)
    else:
        rgba = np.tile(mcolors.to_rgba(color), (len(frame), 1))
    if probability is not None:
        rgba[:, 3] = np.clip(frame[probability].to_numpy(float), 0.1, 1.0)
    collection = LineCollection(
        segments,
        colors=rgba,
        linewidths=linewidth,
        rasterized=rasterized,
        zorder=zorder,
    )
    ax.add_collection(collection)
    return collection


def _draw_reads(
    ax: plt.Axes,
    df: pd.DataFrame,
    selected: Sequence[str],
    colors: Mapping[str, str],
    xlim: tuple[float, float],
    linewidth: float,
) -> None:
    xmin, xmax = xlim
    bg, fg, fg_colors = [], [], []
    selected = set(selected)
    for row in df.itertuples(index=False):
        a, b = max(float(row.span_start), xmin), min(float(row.span_end), xmax)
        if b > a:
            bg.append((a, b, float(row.y)))
        for start, end, motif in row.tvr_calls:
            a, b = max(start, xmin), min(end, xmax)
            if b > a:
                fg.append((a, b, float(row.y)))
                fg_colors.append(colors.get(motif, OTHER_COLOR) if motif in selected else OTHER_COLOR)
    if bg:
        arr = np.asarray(bg)
        _line_collection(ax, arr[:, 0], arr[:, 1], arr[:, 2], READ_BG, linewidth, "round", 1)
    if fg:
        arr = np.asarray(fg)
        _line_collection(ax, arr[:, 0], arr[:, 1], arr[:, 2], fg_colors, linewidth, "butt", 3)


def _draw_consensus(
    ax: plt.Axes,
    bg: pd.DataFrame,
    fg: pd.DataFrame,
    linewidth: float,
    *,
    left_flat: bool,
) -> None:
    if not bg.empty:
        _line_collection(
            ax,
            bg["start"].to_numpy(float), bg["end"].to_numpy(float),
            bg["y"].to_numpy(float), bg["color"].tolist(),
            linewidth, "butt" if left_flat else "round", 1,
        )
        if left_flat:
            xs = bg["end"].to_numpy(float)
            ys = bg["y"].to_numpy(float)
            ax.scatter(
                xs, ys, s=(linewidth * 0.95) ** 2, c=bg["color"].tolist(),
                edgecolors="none", rasterized=False, zorder=2,
            )
    if not fg.empty:
        _line_collection(
            ax,
            fg["start"].to_numpy(float), fg["end"].to_numpy(float),
            fg["y"].to_numpy(float), fg["color"].tolist(),
            linewidth, "butt", 3,
        )


def _collect_sites(
    df: pd.DataFrame,
    xlim: tuple[float, float],
    *,
    methylated: bool,
) -> pd.DataFrame:
    xmin, xmax = xlim
    xs: list[float] = []
    ys: list[float] = []
    color = METHYL_COLOR if methylated else UNMETHYL_COLOR
    for row in df.itertuples(index=False):
        y = float(row.y)
        positions = row.methyl_pos if methylated else row.unmethyl_pos
        if len(positions):
            visible = positions[(positions >= xmin) & (positions <= xmax)]
            xs.extend(visible)
            ys.extend(np.full(len(visible), y))
    return pd.DataFrame({"x": xs, "y": ys, "color": color})


def _draw_5mc(
    ax: plt.Axes,
    df: pd.DataFrame,
    xlim: tuple[float, float],
    *,
    show_unmethylated: bool,
) -> None:
    xmin, xmax = xlim
    bg = []
    for row in df.itertuples(index=False):
        y = float(row.y)
        a, b = max(float(row.span_start), xmin), min(float(row.span_end), xmax)
        if b > a:
            bg.append((a, b, y))
    if bg:
        arr = np.asarray(bg)
        _line_collection(
            ax, arr[:, 0], arr[:, 1], arr[:, 2], MOD_READ_BG,
            _linewidth(ax, 0.14), "round", 1,
        )
    if show_unmethylated:
        unmethylated = _collect_sites(df, xlim, methylated=False)
        _plot_vertical_calls(
            ax, unmethylated, x="x", y="y", color="color", probability=None,
            height=0.72, linewidth=0.4,
            rasterized=False, zorder=3,
        )
    methylated = _collect_sites(df, xlim, methylated=True)
    _plot_vertical_calls(
        ax, methylated, x="x", y="y", color="color", probability=None,
        height=0.72, linewidth=0.4,
        rasterized=False, zorder=4,
    )


def _legend_columns(labels: Sequence[str], width: float, requested: int) -> int:
    if not labels:
        return 1
    if requested:
        return min(requested, len(labels))
    usable_width = max(1.0, width - 1.25)
    entry_width = max(0.85, 0.063 * max(map(len, labels)) + 0.48)
    return max(1, min(len(labels), int(usable_width // entry_width)))


def _legend_height(item_count: int, columns: int) -> float:
    if item_count == 0:
        return 0.0
    rows = int(np.ceil(item_count / columns))
    line_height = STYLE.legend_size * 1.55 / 72.0
    return 0.22 + rows * line_height


def _render(
    panel: Panel,
    motifs: Sequence[str],
    colors: Mapping[str, str],
    config: TVRReadPlotter,
    xlim: tuple[float, float],
    multi_sample: bool,
    title_prefix: str,
) -> plt.Figure:
    prepared = _prepare_panel(panel, motifs, colors, config, xlim, multi_sample)

    xmin, xmax = xlim
    visible = {
        motif
        for calls in prepared.tvr["tvr_calls"]
        for start, end, motif in calls
        if end > xmin and start < xmax
    } | prepared.cons_visible
    handles = [
        Line2D([0], [0], color=colors[m], lw=2.8, label=m)
        for m in motifs if m in visible and m in colors
    ]
    if visible.difference(motifs):
        handles.append(Line2D([0], [0], color=OTHER_COLOR, lw=2.8, label="Others"))

    legend_labels = [handle.get_label() for handle in handles]
    ncol = _legend_columns(legend_labels, config.width, config.legend_ncol)
    legend_h = _legend_height(len(handles), ncol)

    main_h = max(0.90, 0.065 * prepared.total_y)
    fig_height = config.height or main_h + legend_h + 0.81

    fig = plt.figure(figsize=(config.width, fig_height), layout="constrained")
    if handles:
        grid = fig.add_gridspec(2, 1, height_ratios=[main_h, legend_h], hspace=0.02)
        ax = fig.add_subplot(grid[0, 0])
        legend_ax = fig.add_subplot(grid[1, 0])
        legend_ax.set_axis_off()
    else:
        ax = fig.add_subplot(111)
        legend_ax = None

    _decorate(ax, xlim, prepared.total_y, prepared.labels)
    ax.set_title(_plot_title(panel, title_prefix), loc="left", fontweight="bold", pad=4)

    if legend_ax is not None:
        legend_ax.legend(
            handles=row_major(handles, ncol),
            title="TVR motif",
            ncol=ncol,
            loc="center",
            frameon=False,
            handlelength=1.55,
            handletextpad=0.55,
            columnspacing=1.15,
            labelspacing=0.35,
            borderaxespad=0,
            fontsize=STYLE.legend_size,
            title_fontsize=STYLE.legend_size,
        )

    # Finalize constrained geometry before converting y-units to physical linewidths.
    fig.canvas.draw()

    _draw_reads(ax, prepared.tvr, motifs, colors, xlim, _linewidth(ax, config.read_thickness))
    if not prepared.cons_bg.empty:
        _draw_consensus(
            ax, prepared.cons_bg, prepared.cons_fg,
            _linewidth(ax, config.consensus_thickness), left_flat=True,
        )
    if prepared.mods is not None:
        _draw_5mc(ax, prepared.mods, xlim, show_unmethylated=config.show_unmethylated_c)
    if xmin < 0 < xmax:
        ax.axvline(0.0, color=BOUNDARY, linestyle="--", lw=STYLE.axis_linewidth, zorder=10)

    return fig


def plot_tvr_reads(
    tvr_file: str | Path | None = None,
    sample_sheet: str | Path | None = None,
    mods_file: str | Path | None = None,
    consensus_file: str | Path | None = None,
    outdir: str | Path = ".",
    prefix: str = "sample",
    config: TVRReadPlotter | None = None,
    logger=None,
) -> None:
    config = config or TVRReadPlotter()
    logger = logger or log_utils.get_logger()
    outdir = Path(outdir)
    outdir.mkdir(parents=True, exist_ok=True)

    if config.width <= 0 or (config.height is not None and config.height <= 0):
        raise ValueError("Figure width/height must be positive.")
    normalize_figure_format(config.output_format)
    if config.top_tvr < 0:
        raise ValueError("top_tvr must be non-negative.")
    if not np.isfinite(config.read_thickness) or not 0 < config.read_thickness <= 1:
        raise ValueError("read_thickness must be in (0, 1].")
    if not np.isfinite(config.consensus_thickness) or config.consensus_thickness <= 0:
        raise ValueError("consensus_thickness must be positive and finite.")
    if config.legend_ncol < 0:
        raise ValueError("legend_ncol must be non-negative.")
    if sample_sheet is not None and (mods_file is not None or consensus_file is not None):
        raise ValueError("--sample-sheet mode supports TVR files only; --mods/--consensus are not supported.")

    manifest = parse_sample_sheet(primary_file=tvr_file, sample_sheet=sample_sheet, prefix=prefix)
    if len(manifest) > 1 and config.pool_chroms:
        raise ValueError("--pool-chroms is not supported in multi-sample mode.")
    multi_sample = len(manifest) > 1
    requested_xlim = _xlim(config.xlim)

    frames = []
    for sample, path in manifest.items():
        logger.info(f"Loading TVR reads: {sample}")
        frame = _load_tvr(
            Path(path),
            str(sample),
            show_outliers=config.show_outliers,
        )
        if not frame.empty:
            frames.append(frame)
    if not frames:
        logger.warning("No valid TVR reads to plot.")
        return

    tvr = pd.concat(frames, ignore_index=True)
    mods = (
        _load_mods(
            Path(mods_file),
            include_unmethylated=config.show_unmethylated_c,
        )
        if mods_file is not None
        else None
    )
    consensus = _load_consensus(Path(consensus_file)) if consensus_file is not None else None
    tvr = _add_read_lengths(tvr, mods)

    top_by_sample_arm, colors = _top_tvrs(tvr,tuple(manifest),config.top_tvr,requested_xlim,config.tvr_colors,config.canonical_motif,)
    chroms = _select_chroms(tvr, config.chroms)
    tvr = tvr[tvr["chrom"].astype(str).isin(chroms)].copy()
    if tvr.empty:
        logger.warning("No reads remain after --chroms filtering.")
        return
    selected_keys = set(chroms)
    if mods is not None:
        mods = mods[mods["chrom"].isin(selected_keys)].copy()
    if consensus is not None:
        consensus = consensus[consensus["chrom"].isin(selected_keys)].copy()

    color_map_path = outdir / f"{prefix}.tvr_colors.tsv"
    write_motif_color_map(
        color_map_path,
        colors,
        load_motif_color_overrides(config.tvr_colors),
    )
    logger.info(f"Writing TVR color map: {color_map_path}")
    if len(colors) > 60:
        logger.warning(
            f"The figure uses {len(colors)} distinct TVR colors. Colors remain "
            "unique, but a legend of this size may be difficult to distinguish."
        )

    title_prefix = config.sample_label or prefix
    with nature_rc():
        for panel in _panels(tvr, mods, consensus, config.pool_chroms):
            logger.info(f"Plotting {panel.title}")
            motifs = _panel_motifs(top_by_sample_arm, tuple(manifest), panel.arm_key)
            xlim = requested_xlim or _panel_xlim(panel)
            fig = _render(panel, motifs, colors, config, xlim, multi_sample, title_prefix)
            try:
                save_figure(
                    fig,
                    outdir / f"{prefix}.{panel.stem}.reads_plot",
                    output_format=config.output_format,
                )
            finally:
                plt.close(fig)
