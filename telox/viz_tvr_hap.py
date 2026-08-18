# Copyright (C) 2025 Huihui Li <hhui.li@outlook.com>. Licensed under GNU GPL v3.0.
from __future__ import annotations

from collections import Counter
from dataclasses import dataclass
from functools import lru_cache
import math
from pathlib import Path
from typing import Mapping, Sequence
import warnings

import numpy as np
import pandas as pd
from natsort import natsort_keygen, natsorted
from scipy.cluster import hierarchy
from scipy.spatial.distance import squareform

from telox import motiflev
from . import logger as log_utils
from .viz_core import (
    BACKGROUND_COLOR,
    INK,
    MISSING_COLOR,
    OTHER_COLOR,
    STYLE,
    arm_group_key,
    arm_sort_key,
    categorical_color_map,
    chrom_phase_key,
    chrom_phase_label,
    format_kb_axis,
    load_color_overrides,
    load_motif_color_overrides,
    motif_color,
    motif_color_map,
    nature_rc,
    normalize_figure_format,
    parse_motif_blocks,
    parse_sample_sheet,
    plot_intervals,
    require_schema,
    resolve_requested_values,
    row_major,
    safe_name,
    display_chrom,
    save_figure,
    split_csv,
    write_motif_color_map,
)
import matplotlib as mpl
import matplotlib.pyplot as plt
from matplotlib.collections import PolyCollection, QuadMesh
from matplotlib.lines import Line2D
from matplotlib.patches import Patch


PLOT_REQUIRED_COLUMNS = frozenset({
    "chrom", "chrom_phase", "arm", "tvr_hap", "read_support",
    "hap_frequency", "hap_length", "motif_blocks",
})
ROW_LABEL_JOIN_RULES = {("chrom", "arm"): ""}
ANNOTATION_LABELS = {
    "sample_id": "Sample", "chrom": "Chrom", "chrom_phase": "Phase",
    "arm": "Arm", "tvr_hap": "Haplotype",
}
ANNOTATION_MARKERS = ("o", "s", "D", "^", "v", "P", "X")
CANONICAL_TRACK_BG = "#F4F4F4"
HEATMAP_RASTERIZE_THRESHOLD = 1000


@dataclass(frozen=True, slots=True)
class HapLayout:
    outer_margin: float = 0.10
    tree_width: float = 0.36
    tree_gap: float = 0.025
    label_width_min: float = 0.14
    label_width_max: float = 0.60
    annotation_width: float = 0.105
    annotation_gap: float = 0.10
    panel_gap: float = 0.08
    heatmap_width: float = 1.03
    min_main_height: float = 3.35
    default_main_aspect: float = 0.54
    footer_gap: float = 0.52
    footer_column_gap: float = 0.22
    similarity_legend_width: float = 0.80
    dense_row_threshold_pt: float = 3.6
    title_pad: float = 4.0
    dendrogram_lw: float = 0.50
    cbar_height: float = 0.10


LAYOUT = HapLayout()
NATURAL_SORT_KEY = natsort_keygen()


@dataclass(frozen=True, slots=True)
class HaplotypeFilter:
    max_haplotypes: int
    single_min_frequency: float
    multiple_min_frequency: float

    def __post_init__(self) -> None:
        if (
            isinstance(self.max_haplotypes, bool)
            or not isinstance(self.max_haplotypes, int)
            or self.max_haplotypes < 1
        ):
            raise ValueError("MAX must be a positive integer.")
        frequencies = (self.single_min_frequency, self.multiple_min_frequency)
        if not all(math.isfinite(value) and 0 <= value <= 1 for value in frequencies):
            raise ValueError("SINGLE and MULTI must be finite frequencies in [0, 1].")

    @classmethod
    def parse(cls, value: object) -> HaplotypeFilter:
        if isinstance(value, cls):
            return value
        if isinstance(value, str):
            text = value.strip()
            if text.casefold() == "hprc2":
                return cls(2, 0.4, 0.2)
            parts = [part.strip() for part in text.split(",")]
        elif isinstance(value, Sequence):
            parts = list(value)
        else:
            raise TypeError("Expected HPRC2 or MAX,SINGLE,MULTI.")
        if len(parts) != 3:
            raise ValueError("Expected HPRC2 or MAX,SINGLE,MULTI (for example 2,0.4,0.2).")
        try:
            maximum = int(parts[0])
            if str(parts[0]).strip() != str(maximum):
                raise ValueError
            single = float(parts[1])
            multiple = float(parts[2])
        except (TypeError, ValueError) as exc:
            raise ValueError("Expected a positive integer MAX and numeric SINGLE,MULTI frequencies.") from exc
        return cls(maximum, single, multiple)


@dataclass(slots=True)
class TVRHapPlotter:
    canonical_motif: str = "TTAGGG"
    width: float = STYLE.inches(STYLE.double_column_mm)
    height: float | None = None
    chroms: str | None = "all"
    pool_arms: bool = False
    pool_chroms: bool = False
    top_tvr: int = 10
    tvr_colors: Path | None = None
    xlim: float | str | None = 2000.0
    track_thickness: float = 0.78
    cluster_rows: bool = True
    cluster_proximal_bp: int = 250
    no_heatmap: bool = False
    heatmap_cmap: str = "Blues"
    row_label_columns: str | Sequence[str] = "auto"
    strip_chr_prefix: bool = True
    annotation_columns: str | Sequence[str] = ("chrom")
    annotation_file: Path | None = None
    annotation_colors: Path | None = None
    legend_ncol: int = 0
    rasterize: bool | None = None
    hap_filter: HaplotypeFilter | str | Sequence[float] | None = None
    sample_label: str | None = None
    output_format: str = "pdf"

    def __post_init__(self) -> None:
        if self.hap_filter is not None:
            self.hap_filter = HaplotypeFilter.parse(self.hap_filter)
        if isinstance(self.xlim, str):
            value = self.xlim.strip()
            if value.casefold() == "none":
                self.xlim = None
            else:
                try:
                    self.xlim = float(value)
                except ValueError as exc:
                    raise ValueError("xlim must be a positive number or 'none'.") from exc

    @classmethod
    def from_params(cls, **kwargs):
        return cls(**{
            key: value
            for key, value in kwargs.items()
            if key in cls.__dataclass_fields__ and value is not None
        })


@dataclass(frozen=True, slots=True)
class MotifStyle:
    canonical: str
    selected: tuple[str, ...]
    colors: Mapping[str, str]
    background: str
    other_label: str = "Others"


@dataclass(frozen=True, slots=True)
class AnnotationTrack:
    label: str
    marker: str
    values: tuple[str, ...]
    colors: tuple[str, ...]
    palette: Mapping[str, str]
    show_legend: bool


@dataclass(slots=True)
class PreparedHaplotypes:
    output_stem: str
    title_suffix: str
    row_labels: list[str]
    row_ids: list[str]
    backgrounds: pd.DataFrame
    intervals: pd.DataFrame
    motif_legend: tuple[str, ...]
    motif_colors: Mapping[str, str]
    similarity: np.ndarray | None
    linkage: np.ndarray | None
    annotations: tuple[AnnotationTrack, ...]
    xlim: tuple[float, float]


@dataclass(frozen=True, slots=True)
class FooterBlock:
    key: str
    title: str
    labels: tuple[str, ...] = ()
    kind: str = "annotation"
    columns: int | None = None
    fixed_width: float | None = None
    fixed_height: float | None = None


@dataclass(frozen=True, slots=True)
class FooterPlacement:
    key: str
    x: float
    width: float
    height: float
    columns: int = 0


@dataclass(frozen=True, slots=True)
class FooterPlan:
    placements: tuple[FooterPlacement, ...]
    width: float
    height: float


@dataclass(frozen=True, slots=True)
class Geometry:
    figure_height: float
    main_y: float
    main_height: float
    tree: tuple[float, float]
    row_labels: tuple[float, float]
    annotation_tracks: tuple[tuple[float, float], ...]
    tracks: tuple[float, float]
    heatmap: tuple[float, float] | None
    footer_blocks: tuple[FooterPlacement, ...]
    footer_y: float
    footer_height: float
    row_pitch: float
    dense_rows: bool
    vertical_scale: float


def _require_text(frame: pd.DataFrame, columns: Sequence[str], source: Path) -> None:
    for column in columns:
        bad = frame[column].isna() | frame[column].astype(str).str.strip().eq("")
        if bad.any():
            rows = ", ".join(str(index + 2) for index in frame.index[bad][:8])
            raise ValueError(f"{source}: column {column!r} is empty at input row(s): {rows}")


def _stable_ids(frame: pd.DataFrame, columns: Sequence[str], sep: str = "-") -> list[str]:
    base = frame[list(columns)].astype(str).agg(sep.join, axis=1)
    occurrence = base.groupby(base, sort=False).cumcount()
    duplicated = base.duplicated(keep=False)
    return base.where(~duplicated, base + sep + occurrence.astype(str)).tolist()


def _read_annotation_table(path: str | Path | None) -> pd.DataFrame | None:
    if path is None:
        return None
    path = Path(path)
    if not path.exists():
        raise ValueError(f"Annotation file not found: {path}")
    table = pd.read_csv(path, sep=None, engine="python", dtype=str)
    if table.shape[1] < 2:
        raise ValueError("Annotation file must contain an ID column and at least one annotation column.")
    table = table.set_index(table.columns[0])
    if table.index.has_duplicates:
        duplicates = table.index[table.index.duplicated(keep=False)].unique().tolist()[:8]
        raise ValueError("Duplicate annotation ID(s): " + ", ".join(map(str, duplicates)))
    return table


def _parse_tokens(frame: pd.DataFrame, source: Path) -> pd.Series:
    values = frame["motif_blocks"].astype(str)
    cache: dict[str, tuple[tuple[str, int], ...]] = {}
    for raw in pd.unique(values):
        try:
            cache[raw] = parse_motif_blocks(raw)
        except ValueError as exc:
            raise ValueError(f"{source}: invalid motif_blocks value {raw!r}: {exc}") from exc
    return values.map(cache)


def _validate_frame(frame: pd.DataFrame, source: Path) -> pd.DataFrame:
    require_schema(frame, PLOT_REQUIRED_COLUMNS, source=source)
    if frame.empty:
        return frame.copy()

    result = frame.copy()
    text_columns = ("chrom", "chrom_phase", "arm", "tvr_hap", "motif_blocks")
    _require_text(result, text_columns, source)
    for column in text_columns:
        result[column] = result[column].astype(str).str.strip()

    for column in ("read_support", "hap_frequency", "hap_length"):
        result[column] = pd.to_numeric(result[column], errors="raise")
        if result[column].isna().any():
            rows = ", ".join(str(index + 2) for index in result.index[result[column].isna()][:8])
            raise ValueError(f"{source}: column {column!r} is NA at input row(s): {rows}")

    if (result["read_support"] < 0).any():
        raise ValueError(f"{source}: read_support must be non-negative.")
    if (result["hap_length"] <= 0).any():
        raise ValueError(f"{source}: hap_length must be positive.")
    if not result["hap_frequency"].between(0, 1).all():
        raise ValueError(f"{source}: hap_frequency must be in [0, 1].")

    result["tokens"] = _parse_tokens(result, source)
    return result


def load_haplotypes(
    consensus_file: str | Path | None,
    sample_sheet: str | Path | None,
    prefix: str,
) -> pd.DataFrame:

    manifest = parse_sample_sheet(primary_file=consensus_file,sample_sheet=sample_sheet,prefix=prefix,)
    frames: list[pd.DataFrame] = []
    for sample, path in manifest.items():
        path = Path(path)
        frame = pd.read_csv(path, sep="\t", na_values=".")
        frame = _validate_frame(frame, path)
        frame["sample_id"] = str(sample)
        frames.append(frame)
    return pd.concat(frames, ignore_index=True) if frames else pd.DataFrame()


def _validate_config(config: TVRHapPlotter) -> None:
    if config.width <= 0 or (config.height is not None and config.height <= 0):
        raise ValueError("Figure width/height must be positive.")
    normalize_figure_format(config.output_format)
    if config.xlim is not None and config.xlim <= 0:
        raise ValueError("xlim must be positive or None.")
    if config.top_tvr < 0:
        raise ValueError("top_tvr must be non-negative.")
    if not 0 < config.track_thickness <= 1:
        raise ValueError("track_thickness must be in (0, 1].")
    if config.cluster_proximal_bp <= 0:
        raise ValueError("cluster_proximal_bp must be positive.")
    if config.legend_ncol < 0:
        raise ValueError("legend_ncol must be non-negative.")
    if config.heatmap_cmap not in mpl.colormaps:
        raise ValueError(f"Unknown Matplotlib colormap: {config.heatmap_cmap!r}")


def _filter_haplotypes(frame: pd.DataFrame, config: TVRHapPlotter) -> pd.DataFrame:
    policy = config.hap_filter
    if policy is None:
        return frame.copy()
    selected = frame.copy()
    selected["_phase_group"] = selected["chrom_phase"].map(chrom_phase_key)
    group = ["sample_id", "chrom", "_phase_group", "arm"]
    distinct_haplotypes = selected.groupby(group, observed=True)["tvr_hap"].transform("nunique")
    selected["_min_hap_frequency"] = np.where(
        distinct_haplotypes.eq(1),
        policy.single_min_frequency,
        policy.multiple_min_frequency,
    )
    selected = (
        selected.sort_values(
            group + ["hap_frequency", "read_support", "tvr_hap"],
            ascending=[True] * len(group) + [False, False, True],
        )
        .groupby(group, observed=True)
        .head(policy.max_haplotypes)
        .copy()
    )
    selected = selected[
        selected["hap_frequency"].astype(float) >= selected["_min_hap_frequency"]
    ].copy()
    return selected.drop(columns=["_phase_group", "_min_hap_frequency"])


def _filter_requested_rows(frame: pd.DataFrame, config: TVRHapPlotter) -> pd.DataFrame:
    observed = natsorted(frame["chrom"].astype(str).unique())
    resolved = resolve_requested_values(config.chroms, observed, label="chrom value")
    return frame[frame["chrom"].astype(str).isin(resolved)].copy()

@dataclass(frozen=True, slots=True)
class ArmGroup:
    key: str
    label: str
    values: tuple[str, ...]


@dataclass(frozen=True, slots=True)
class PanelSpec:
    arm_key: str
    output_stem: str
    title_suffix: str
    rows: pd.DataFrame


@dataclass(frozen=True, slots=True)
class TopTVRStats:
    by_sample_arm: Mapping[tuple[str, str], tuple[str, ...]]
    arm_counts: Mapping[str, Counter[str]]
    candidates: tuple[str, ...]
    candidate_set: frozenset[str]
    global_counts: Counter[str]


def _arm_groups(frame: pd.DataFrame) -> tuple[ArmGroup, ...]:

    observed = sorted(frame["arm"].astype(str).unique(), key=arm_sort_key)
    lowered = {value.lower() for value in observed}
    is_pq = bool(lowered & {"p", "q"})

    if is_pq:
        definitions = (("p", "p-arm", "left"),("q", "q-arm", "right"),)
    else:
        definitions = (("L", "L-arm", "left"),("R", "R-arm", "right"),)

    groups: list[ArmGroup] = []
    consumed: set[str] = set()
    for key, label, side in definitions:
        values = tuple(value for value in observed if arm_group_key(value) == side)
        if values:
            groups.append(ArmGroup(key=key, label=label, values=values))
            consumed.update(values)

    for value in observed:
        if value not in consumed:
            groups.append(ArmGroup(key=value, label=f"{value}-arm", values=(value,)))
    return tuple(groups)


def _pooled_arm_group(groups: Sequence[ArmGroup]) -> ArmGroup:
    keys = [group.key for group in groups]
    values = tuple(value for group in groups for value in group.values)
    if keys == ["p", "q"]:
        return ArmGroup("pq", "pq-arms", values)
    if keys == ["L", "R"]:
        return ArmGroup("LR", "LR-arms", values)
    return ArmGroup("pooled", "pooled-arms", values)


def _trim_tokens(tokens: tuple[tuple[str, int], ...], bp: int) -> tuple[tuple[str, int], ...]:
    if bp <= 0:
        return tokens

    kept: list[tuple[str, int]] = []
    used = 0
    for motif, count in tokens:
        motif_length = len(motif)
        token_length = motif_length * int(count)
        if used + token_length <= bp:
            kept.append((motif, int(count)))
            used += token_length
            continue

        remaining = bp - used
        allowed_count = remaining // motif_length
        if allowed_count > 0:
            kept.append((motif, allowed_count))
        break

    return tuple(kept)


def _count_tvr_copies(frame: pd.DataFrame,canonical: str,xlim: float | None,) -> Counter[str]:
    counts: Counter[str] = Counter()
    for tokens in frame["tokens"]:
        ranking_tokens = (_trim_tokens(tokens, int(math.floor(xlim))) if xlim is not None else tokens)
        for motif, count in ranking_tokens:
            if motif != canonical:
                counts[motif] += int(count)
    return counts


def _rank_counts(counts: Counter[str], limit: int | None = None) -> tuple[str, ...]:
    if not counts or limit == 0:
        return ()
    ordered = sorted(counts, key=counts.get, reverse=True)
    return tuple(ordered if limit is None else ordered[:limit])


def _unique_union(groups: Sequence[Sequence[str]]) -> tuple[str, ...]:
    seen: set[str] = set()
    result: list[str] = []
    for group in groups:
        for value in group:
            if value not in seen:
                seen.add(value)
                result.append(value)
    return tuple(result)


def _top_tvr_stats(
    frame: pd.DataFrame,
    config: TVRHapPlotter,
    arm_groups: Sequence[ArmGroup],
) -> TopTVRStats:

    by_sample_arm: dict[tuple[str, str], tuple[str, ...]] = {}
    arm_counts: dict[str, Counter[str]] = {group.key: Counter() for group in arm_groups}
    global_counts: Counter[str] = Counter()
    candidates: list[str] = []
    candidate_seen: set[str] = set()
    samples = natsorted(frame["sample_id"].astype(str).unique())
    sample_groups = {
        str(sample): rows
        for sample, rows in frame.groupby(
            frame["sample_id"].astype(str),
            sort=False,
            observed=True,
        )
    }
    arm_group_by_value = {
        str(value): group.key
        for group in arm_groups
        for value in group.values
    }

    for sample in samples:
        sample_rows = sample_groups[str(sample)]
        group_keys = sample_rows["arm"].astype(str).map(arm_group_by_value)
        grouped_rows = {
            str(key): rows
            for key, rows in sample_rows.groupby(
                group_keys,
                sort=False,
                observed=True,
            )
        }
        for group in arm_groups:
            arm_rows = grouped_rows.get(group.key, sample_rows.iloc[0:0])
            counts = _count_tvr_copies(arm_rows, config.canonical_motif, config.xlim)
            arm_counts[group.key].update(counts)
            global_counts.update(counts)
            local_top = _rank_counts(counts, config.top_tvr)
            by_sample_arm[(sample, group.key)] = local_top
            for motif in local_top:
                if motif not in candidate_seen:
                    candidate_seen.add(motif)
                    candidates.append(motif)

    if config.pool_arms and arm_groups:
        pooled = _pooled_arm_group(arm_groups)
        pooled_counts: Counter[str] = Counter()
        for group in arm_groups:
            pooled_counts.update(arm_counts[group.key])
        arm_counts[pooled.key] = pooled_counts
        for sample in samples:
            by_sample_arm[(sample, pooled.key)] = _unique_union(
                [by_sample_arm.get((sample, group.key), ()) for group in arm_groups]
            )

    return TopTVRStats(
        by_sample_arm=by_sample_arm,
        arm_counts=arm_counts,
        candidates=tuple(candidates),
        candidate_set=frozenset(candidates),
        global_counts=global_counts,
    )


def _shared_motif_colors(stats: TopTVRStats, config: TVRHapPlotter) -> dict[str, str]:
    overrides = load_motif_color_overrides(config.tvr_colors)
    canonical_color = motif_color(config.canonical_motif, overrides)
    colors: dict[str, str] = {}
    if stats.candidates:
        candidate_counts = Counter(
            {motif: stats.global_counts.get(motif, 0) for motif in stats.candidates}
        )
        mapped, _ = motif_color_map(
            candidate_counts,
            len(stats.candidates),
            overrides,
            exclude=(config.canonical_motif,),
            reserved_colors=(canonical_color,),
        )
        colors.update(mapped)
    colors[config.canonical_motif] = canonical_color
    colors["Other TVR motif"] = OTHER_COLOR
    return colors


def _panel_motif_style(
    spec: PanelSpec,
    stats: TopTVRStats,
    colors: Mapping[str, str],
    config: TVRHapPlotter,
    *,
    multi_sample: bool,
) -> MotifStyle:
    if multi_sample:
        selected = stats.candidates
        ranked_here = _rank_counts(stats.arm_counts.get(spec.arm_key, Counter()))
        ordered = tuple(motif for motif in ranked_here if motif in stats.candidate_set)
        ordered = _unique_union((ordered, selected))
    else:
        sample = str(spec.rows["sample_id"].iloc[0])
        selected = stats.by_sample_arm.get((sample, spec.arm_key), ())
        ordered = selected

    return MotifStyle(
        canonical=config.canonical_motif,
        selected=ordered,
        colors=colors,
        background=CANONICAL_TRACK_BG,
    )


def _distance_matrix(
    sequences: Sequence[tuple[tuple[str, int], ...]], proximal_bp: int
) -> np.ndarray | None:
    if len(sequences) < 2:
        return None
    trimmed = [_trim_tokens(sequence, proximal_bp) for sequence in sequences]
    unique = list(dict.fromkeys(trimmed))
    lookup = {sequence: index for index, sequence in enumerate(unique)}
    inverse = np.fromiter((lookup[sequence] for sequence in trimmed), int, len(trimmed))

    aligner = motiflev.Aligner()
    distances = np.zeros((len(unique), len(unique)))
    for i, first in enumerate(unique):
        for j in range(i + 1, len(unique)):
            distances[i, j] = distances[j, i] = float(aligner.distance(first, unique[j]))
    return distances[np.ix_(inverse, inverse)]


def _row_label_columns(value: str | Sequence[str], *, multi_sample: bool) -> tuple[str, ...]:
    if isinstance(value, str):
        text = value.strip()
        if text.lower() == "none":
            return ()
        if text.lower() == "auto":
            return ("sample_id", "chrom", "arm") if multi_sample else ("chrom", "arm")
    columns = split_csv(value)
    duplicates = [name for name, count in Counter(columns).items() if count > 1]
    if duplicates:
        raise ValueError("Duplicate row-label column(s): " + ", ".join(duplicates))
    return columns


def _annotation_columns(value: str | Sequence[str]) -> tuple[str, ...]:
    if isinstance(value, str) and value.strip().casefold() == "none":
        return ()
    return split_csv(value)


def _row_labels(rows: pd.DataFrame, columns: Sequence[str], strip_chr: bool) -> list[str]:
    labels: list[str] = []
    for values in rows.loc[:, list(columns)].itertuples(index=False, name=None):
        parts: list[str] = []
        previous_column: str | None = None
        for column, value in zip(columns, values, strict=True):
            text = str(value).strip()
            if column == "chrom_phase":
                text = chrom_phase_label(value)
            if column == "chrom" and strip_chr:
                text = display_chrom(text)
            if not text:
                continue
            if previous_column is not None:
                parts.append(ROW_LABEL_JOIN_RULES.get((previous_column, column), "-"))
            parts.append(text)
            previous_column = column
        labels.append("".join(parts))
    return labels

def _row_ids(rows: pd.DataFrame) -> list[str]:
    identity = pd.DataFrame({
        "sample_id": rows["sample_id"].astype(str),
        "chrom_arm": (rows["chrom"].astype(str) + rows["arm"].astype(str).str.lower()),
        "chrom_phase": rows["chrom_phase"].astype(str),
        "tvr_hap": rows["tvr_hap"].astype(str),
    })

    return _stable_ids(identity,("sample_id", "chrom_arm", "chrom_phase", "tvr_hap"),)

def _external_annotations(
    rows: pd.DataFrame,
    row_ids: Sequence[str],
    source: pd.DataFrame | None,
) -> pd.DataFrame | None:
    if source is None:
        return None
    sample_ids = rows["sample_id"].astype(str).tolist()
    if any(row_id in source.index for row_id in row_ids):
        keys = list(row_ids)
    elif any(sample in source.index for sample in sample_ids):
        keys = sample_ids
    else:
        keys = list(row_ids)
    return source.reindex(keys).fillna("NA").astype(str).reset_index(drop=True)


def _annotation_palettes(
    frame: pd.DataFrame,
    columns: Sequence[str],
    external: pd.DataFrame | None,
    overrides: Mapping[str, str],
) -> tuple[dict[str, dict[str, str]], dict[str, dict[str, str]]]:
    builtin = {
        column: categorical_color_map(
            frame[column].fillna("NA").astype(str).tolist(),
            column=column,
            overrides=overrides,
        )
        for column in columns
    }
    custom = {} if external is None else {
        str(column): categorical_color_map(
            external[column].fillna("NA").astype(str).tolist(),
            column=str(column),
            overrides=overrides,
        )
        for column in external.columns
    }
    return builtin, custom


def _unique_label(label: str, used: set[str]) -> str:
    if label not in used:
        used.add(label)
        return label
    index = 2
    while f"{label} ({index})" in used:
        index += 1
    result = f"{label} ({index})"
    used.add(result)
    return result


def _show_annotation_legend(
    column: str,
    values: pd.Series,
    *,
    multi_sample: bool,
    builtin: bool,
) -> bool:
    observed = values[values.ne("NA")].nunique(dropna=False)
    if observed <= 1:
        return False
    return multi_sample if builtin and column == "chrom" else True


def _prepare_annotations(
    rows: pd.DataFrame,
    row_ids: Sequence[str],
    columns: Sequence[str],
    external_source: pd.DataFrame | None,
    builtin_palettes: Mapping[str, Mapping[str, str]],
    external_palettes: Mapping[str, Mapping[str, str]],
    *,
    multi_sample: bool,
) -> tuple[AnnotationTrack, ...]:
    items: list[tuple[str, str, pd.Series, Mapping[str, str], bool]] = []
    for column in columns:
        items.append((
            column,
            ANNOTATION_LABELS.get(column, column.replace("_", " ").title()),
            rows[column].fillna("NA").astype(str).reset_index(drop=True),
            builtin_palettes[column],
            True,
        ))

    external = _external_annotations(rows, row_ids, external_source)
    if external is not None:
        for column in external.columns:
            items.append((
                str(column),
                str(column),
                external[column].fillna("NA").astype(str).reset_index(drop=True),
                external_palettes[str(column)],
                False,
            ))

    tracks: list[AnnotationTrack] = []
    used_labels: set[str] = set()
    for index, (column, label, values, palette, builtin) in enumerate(items):
        label = _unique_label(label, used_labels)
        colors = tuple(palette.get(value, MISSING_COLOR) for value in values)
        tracks.append(AnnotationTrack(
            label=label,
            marker=ANNOTATION_MARKERS[index % len(ANNOTATION_MARKERS)],
            values=tuple(values),
            colors=colors,
            palette=palette,
            show_legend=_show_annotation_legend(
                column, values, multi_sample=multi_sample, builtin=builtin
            ),
        ))
    return tuple(tracks)


def _prepare_intervals(
    rows: pd.DataFrame,
    style: MotifStyle,
    xlim: float | None,
) -> tuple[pd.DataFrame, pd.DataFrame, tuple[str, ...], tuple[float, float]]:
    backgrounds: list[tuple[float, float, float, str]] = []
    intervals: list[tuple[float, float, float, str, str]] = []
    seen: set[str] = set()
    max_length = 0.0
    selected = set(style.selected)

    for y, row in enumerate(rows.itertuples(index=False)):
        full_length = float(row.hap_length)
        visible_length = min(full_length, xlim) if xlim is not None else full_length
        backgrounds.append((0.0, visible_length, float(y), style.background))

        position = 0.0
        for motif, count in row.tokens:
            end = position + len(motif) * count
            display_end = min(end, xlim) if xlim is not None else end
            if motif != style.canonical and display_end > position:
                label = motif if motif in selected else style.other_label
                intervals.append((position, display_end, float(y), style.colors.get(label, OTHER_COLOR), label))
                seen.add(label)
            position = end
            if xlim is not None and position >= xlim:
                break
        max_length = max(max_length, full_length)

    legend = tuple(
        [motif for motif in style.selected if motif in seen]
        + ([style.other_label] if style.other_label in seen else [])
    )
    display_max = float(xlim if xlim is not None else max_length)
    return (
        pd.DataFrame(backgrounds, columns=["start", "end", "y", "color"]),
        pd.DataFrame(intervals, columns=["start", "end", "y", "color", "motif"]),
        legend,
        (0.0, display_max),
    )


def _build_panel_specs(
    selected: pd.DataFrame,
    config: TVRHapPlotter,
    arm_groups: Sequence[ArmGroup],
    *,
    multi_sample: bool,
) -> tuple[PanelSpec, ...]:
    groups = (
        (_pooled_arm_group(arm_groups),)
        if config.pool_arms and arm_groups
        else tuple(arm_groups)
    )
    specs: list[PanelSpec] = []

    if multi_sample and not config.pool_chroms:
        for chrom in natsorted(selected["chrom"].astype(str).unique()):
            chrom_rows = selected[selected["chrom"].astype(str).eq(chrom)]
            for group in groups:
                rows = chrom_rows[chrom_rows["arm"].astype(str).isin(group.values)].copy()
                if rows.empty:
                    continue
                safe_chrom = safe_name(chrom)
                specs.append(PanelSpec(
                    arm_key=group.key,
                    output_stem=f"{safe_chrom}.{group.label}",
                    title_suffix=f"{chrom} {group.label}",
                    rows=rows,
                ))
        return tuple(specs)

    for group in groups:
        rows = selected[selected["arm"].astype(str).isin(group.values)].copy()
        if rows.empty:
            continue
        if multi_sample and config.pool_chroms:
            output_stem = f"pooled_chroms.{group.label}"
            title_suffix = f"pooled chromosomes {group.label}"
        else:
            output_stem = group.label
            title_suffix = group.label
        specs.append(PanelSpec(
            arm_key=group.key,
            output_stem=output_stem,
            title_suffix=title_suffix,
            rows=rows,
        ))
    return tuple(specs)


def prepare_haplotypes(
    frame: pd.DataFrame,
    config: TVRHapPlotter,
    *,
    write_similarity_score: bool = False,
    logger=None,
) -> dict[str, PreparedHaplotypes]:

    _validate_config(config)
    if frame.empty:
        return {}

    multi_sample = frame["sample_id"].nunique(dropna=False) > 1
    filtered = _filter_haplotypes(frame, config)
    if config.hap_filter is not None and logger is not None:
        policy = config.hap_filter
        logger.info(
            "Haplotype filter "
            f"(max={policy.max_haplotypes}, single>={policy.single_min_frequency:g}, "
            f"multi>={policy.multiple_min_frequency:g}): "
            f"{len(frame)} -> {len(filtered)} haplotypes"
        )
    if filtered.empty:
        return {}

    arm_groups = _arm_groups(filtered)
    if not arm_groups:
        return {}

    top_stats = _top_tvr_stats(filtered, config, arm_groups)
    shared_colors = _shared_motif_colors(top_stats, config)

    selected = _filter_requested_rows(filtered, config)
    if selected.empty:
        return {}

    label_columns = _row_label_columns(config.row_label_columns, multi_sample=multi_sample)
    annotation_columns = _annotation_columns(config.annotation_columns)
    unknown_labels = [column for column in label_columns if column not in selected.columns]
    unknown_annotations = [column for column in annotation_columns if column not in selected.columns]
    if unknown_labels:
        raise ValueError("Unknown row-label column(s): " + ", ".join(unknown_labels))
    if unknown_annotations:
        raise ValueError("Unknown annotation column(s): " + ", ".join(unknown_annotations))

    external_source = _read_annotation_table(config.annotation_file)
    overrides = load_color_overrides(config.annotation_colors, key_name="VALUE")
    builtin_palettes, external_palettes = _annotation_palettes(
        selected, annotation_columns, external_source, overrides
    )

    outputs: dict[str, PreparedHaplotypes] = {}
    for spec in _build_panel_specs(
        selected, config, arm_groups, multi_sample=multi_sample
    ):
        natural_columns = (
            ["sample_id", "chrom", "chrom_phase", "tvr_hap"]
            if multi_sample
            else ["chrom", "chrom_phase", "tvr_hap"]
        )
        rows = (
            spec.rows.sort_values(
                natural_columns,
                key=lambda column: column.astype(str).map(NATURAL_SORT_KEY),
                kind="stable",
            )
            .reset_index(drop=True)
        )
        if rows.empty:
            continue

        need_distance = config.cluster_rows or not config.no_heatmap or write_similarity_score
        distance = (
            _distance_matrix(rows["tokens"].tolist(), config.cluster_proximal_bp)
            if need_distance else None
        )
        linkage_matrix = None
        if config.cluster_rows and distance is not None and np.any(distance):
            linkage_matrix = hierarchy.linkage(squareform(distance), method="average")
            order = hierarchy.dendrogram(linkage_matrix, no_plot=True)["leaves"]
            rows = rows.iloc[order].reset_index(drop=True)
            distance = distance[np.ix_(order, order)]

        similarity = (
            np.clip(1.0 - distance, 0.0, 1.0).astype(np.float32, copy=False)
            if distance is not None else None
        )
        labels = _row_labels(rows, label_columns, config.strip_chr_prefix) if label_columns else []
        row_ids = _row_ids(rows)
        annotations = _prepare_annotations(
            rows, row_ids, annotation_columns,
            external_source, builtin_palettes, external_palettes,
            multi_sample=multi_sample,
        )
        style = _panel_motif_style(
            spec, top_stats, shared_colors, config, multi_sample=multi_sample
        )
        backgrounds, intervals, motif_legend, panel_xlim = _prepare_intervals(
            rows, style, config.xlim
        )
        outputs[spec.output_stem] = PreparedHaplotypes(
            output_stem=spec.output_stem,
            title_suffix=spec.title_suffix,
            row_labels=labels,
            row_ids=row_ids,
            backgrounds=backgrounds,
            intervals=intervals,
            motif_legend=motif_legend,
            motif_colors=style.colors,
            similarity=similarity,
            linkage=linkage_matrix,
            annotations=annotations,
            xlim=panel_xlim,
        )
    return outputs


def _label_width(labels: Sequence[str]) -> float:
    if not labels:
        return 0.0
    estimate = max(map(len, labels)) * STYLE.row_label_size * 0.56 / 72.0 + 0.045
    return float(np.clip(estimate, LAYOUT.label_width_min, LAYOUT.label_width_max))


def _text_line_height(fontsize: float, *, linespacing: float = 1.22) -> float:
    return fontsize * linespacing / 72.0


def _top_reserve() -> float:
    title_size = mpl.font_manager.FontProperties(
        size=mpl.rcParams["axes.titlesize"]
    ).get_size_in_points()
    return _text_line_height(title_size) + LAYOUT.title_pad / 72.0 + LAYOUT.outer_margin


def _similarity_legend_metrics() -> tuple[float, float, float, float, float]:
    title_h = _text_line_height(STYLE.legend_size)
    tick_h = _text_line_height(STYLE.annotation_size)
    title_gap = max(0.030, title_h * 0.35)
    tick_gap = max(0.020, tick_h * 0.12)
    bar_h = max(LAYOUT.cbar_height, tick_h * 0.72)
    return title_h, tick_h, title_gap, tick_gap, bar_h


def _similarity_legend_height() -> float:
    title_h, tick_h, title_gap, tick_gap, bar_h = _similarity_legend_metrics()
    return title_h + title_gap + bar_h + tick_gap + tick_h


def _measurement_handles(labels: tuple[str, ...], *, motif: bool) -> list[mpl.artist.Artist]:
    if motif:
        return [Patch(facecolor="0.5", edgecolor="none", label=label) for label in labels]
    return [
        Line2D(
            [0], [0], marker="o", linestyle="None", color="none",
            markerfacecolor="0.5", markeredgecolor="none", markersize=4.5,
            label=label,
        )
        for label in labels
    ]


def _discrete_legend_kwargs(
    handles: Sequence[mpl.artist.Artist],
    *,
    title: str,
    columns: int,
    motif: bool,
) -> dict[str, object]:
    kwargs: dict[str, object] = dict(
        handles=row_major(handles, columns),
        title=title,
        loc="upper left",
        bbox_to_anchor=(0.0, 1.0),
        ncols=columns,
        frameon=False,
        fontsize=STYLE.annotation_size,
        handlelength=1.0 if motif else 0.9,
        handletextpad=0.38 if motif else 0.34,
        columnspacing=0.82 if motif else 0.78,
        labelspacing=0.34,
        borderpad=0.0,
        borderaxespad=0.0,
        title_fontsize=STYLE.legend_size,
        alignment="left",
    )
    if motif:
        kwargs["handleheight"] = 0.75
    return kwargs


@lru_cache(maxsize=512)
def _measure_discrete_legend(
    title: str,
    labels: tuple[str, ...],
    columns: int,
    motif: bool,
) -> tuple[float, float]:
    figure = plt.figure(figsize=(1.0, 1.0), dpi=72)
    try:
        ax = figure.add_axes((0.0, 0.0, 1.0, 1.0))
        ax.set_axis_off()
        handles = _measurement_handles(labels, motif=motif)
        legend = ax.legend(**_discrete_legend_kwargs(
            handles, title=title, columns=columns, motif=motif
        ))
        _style_discrete_legend(legend)
        figure.canvas.draw()
        bbox = legend.get_window_extent(figure.canvas.get_renderer())
        return bbox.width / figure.dpi, bbox.height / figure.dpi
    finally:
        plt.close(figure)


def _footer_block_size(
    block: FooterBlock,
    target_rows: int,
) -> tuple[float, float, int]:
    if block.labels:
        if block.columns is not None:
            columns = min(block.columns, len(block.labels))
        else:
            columns = min(len(block.labels),max(1, math.ceil(len(block.labels) / target_rows)),)
        width, height = _measure_discrete_legend(
            block.title,
            block.labels,
            columns,
            block.kind == "motif",
        )
        return width, height, columns

    if block.fixed_width is None or block.fixed_height is None:
        raise ValueError(f"Footer block {block.key!r} has no measurable content.")
    return block.fixed_width, block.fixed_height, 0

def _plan_footer_blocks(
    blocks: Sequence[FooterBlock],
    available_width: float,
) -> FooterPlan:
    blocks = tuple(blocks)
    if not blocks:
        return FooterPlan((), 0.0, 0.0)

    gap = LAYOUT.footer_column_gap
    flexible = [block for block in blocks if block.labels and block.columns is None]
    max_rows = max((len(block.labels) for block in flexible), default=1)
    choice: list[tuple[float, float, int]] | None = None

    for target_rows in range(1, max_rows + 1):
        sizes = [_footer_block_size(block, target_rows) for block in blocks]
        total_width = sum(width for width, _height, _cols in sizes)
        total_width += gap * (len(sizes) - 1)
        if total_width <= available_width + 1e-9:
            choice = sizes
            break

    if choice is None:
        choice = [_footer_block_size(block, max_rows) for block in blocks]
        total_width = sum(width for width, _height, _cols in choice)
        total_width += gap * (len(choice) - 1)
        if total_width > available_width + 1e-9:
            warnings.warn(
                "Bottom legends exceed the available figure width; consider "
                "increasing --width or shortening annotation labels.",
                RuntimeWarning,
                stacklevel=2,
            )
    else:
        total_width = sum(width for width, _height, _cols in choice)
        total_width += gap * (len(choice) - 1)

    placements: list[FooterPlacement] = []
    x = 0.0
    for block, (width, height, columns) in zip(blocks, choice, strict=True):
        placements.append(FooterPlacement(
            key=block.key,
            x=x,
            width=width,
            height=height,
            columns=columns,
        ))
        x += width + gap

    return FooterPlan(
        placements=tuple(placements),
        width=total_width,
        height=max((placement.height for placement in placements), default=0.0),
    )


def _annotation_legend_entries(track: AnnotationTrack,) -> tuple[tuple[str, str], ...]:
    if not track.show_legend:
        return ()
    observed = set(track.values)
    return tuple(
        (value, color)
        for value, color in track.palette.items()
        if value != "NA" and value in observed
    )


def _annotation_legend_tracks(panel: PreparedHaplotypes) -> tuple[AnnotationTrack, ...]:
    return tuple(track for track in panel.annotations if _annotation_legend_entries(track))


def _footer_blocks(
    panel: PreparedHaplotypes,
    config: TVRHapPlotter,
    *,
    heatmap_visible: bool,
) -> tuple[FooterBlock, ...]:
    blocks: list[FooterBlock] = []
    if panel.motif_legend:
        blocks.append(FooterBlock(
            key="motif",
            title="TVR motif",
            labels=tuple(panel.motif_legend),
            kind="motif",
            columns=config.legend_ncol or None,
        ))

    for index, track in enumerate(_annotation_legend_tracks(panel)):
        blocks.append(FooterBlock(
            key=f"annotation:{index}",
            title=track.label,
            labels=tuple(value for value, _color in _annotation_legend_entries(track)),
        ))

    if heatmap_visible:
        blocks.append(FooterBlock(
            key="similarity",
            title="TVR Similarity",
            kind="continuous",
            fixed_width=LAYOUT.similarity_legend_width,
            fixed_height=_similarity_legend_height(),
        ))
    return tuple(blocks)


def figure_geometry(panel: PreparedHaplotypes, config: TVRHapPlotter) -> Geometry:
    tree_width = LAYOUT.tree_width if panel.linkage is not None else 0.03
    tree_gap = LAYOUT.tree_gap if panel.linkage is not None else 0.0
    label_width = _label_width(panel.row_labels)
    annotation_width = len(panel.annotations) * LAYOUT.annotation_width
    heatmap_visible = not config.no_heatmap and panel.similarity is not None
    heatmap_width = LAYOUT.heatmap_width if heatmap_visible else 0.0
    heatmap_gap = LAYOUT.panel_gap if heatmap_visible else 0.0

    occupied = (
        2 * LAYOUT.outer_margin + tree_width + tree_gap + label_width + annotation_width
        + LAYOUT.annotation_gap + heatmap_gap + heatmap_width
    )
    track_width = config.width - occupied
    x = LAYOUT.outer_margin
    tree = (x, tree_width)
    x += tree_width + tree_gap
    row_labels = (x, label_width)
    x += label_width
    annotation_tracks = tuple(
        (x + index * LAYOUT.annotation_width, LAYOUT.annotation_width)
        for index in range(len(panel.annotations))
    )
    x += annotation_width + LAYOUT.annotation_gap
    tracks = (x, track_width)
    x += track_width
    heatmap = None
    if heatmap_visible:
        x += LAYOUT.panel_gap
        heatmap = (x, LAYOUT.heatmap_width)

    footer_left = LAYOUT.outer_margin
    footer_right = config.width - LAYOUT.outer_margin
    footer_width = max(0.0, footer_right - footer_left)
    footer_blocks = _footer_blocks(panel, config, heatmap_visible=heatmap_visible)
    footer_plan = _plan_footer_blocks(footer_blocks, footer_width)

    footer_group_left = footer_left + max(0.0, (footer_width - footer_plan.width) / 2)
    footer_placements = tuple(
        FooterPlacement(
            key=placement.key,
            x=footer_group_left + placement.x,
            width=placement.width,
            height=placement.height,
            columns=placement.columns,
        )
        for placement in footer_plan.placements
    )
    footer_height = footer_plan.height

    preferred_main = max(LAYOUT.min_main_height, config.width * LAYOUT.default_main_aspect)
    top_reserve = _top_reserve()
    bottom_margin = LAYOUT.outer_margin
    preferred_fixed = bottom_margin + footer_height + LAYOUT.footer_gap + top_reserve
    if config.height is None:
        figure_height = preferred_fixed + preferred_main
        main_height = preferred_main
        footer_gap = LAYOUT.footer_gap
        vertical_scale = 1.0
    elif config.height > preferred_fixed:
        figure_height = config.height
        main_height = figure_height - preferred_fixed
        footer_gap = LAYOUT.footer_gap
        vertical_scale = 1.0
    else:
        figure_height = config.height
        template = preferred_fixed + preferred_main
        vertical_scale = figure_height / template
        main_height = preferred_main * vertical_scale
        footer_gap = LAYOUT.footer_gap * vertical_scale
        footer_height *= vertical_scale
        bottom_margin *= vertical_scale

    footer_y = bottom_margin
    main_y = footer_y + footer_height + footer_gap
    row_pitch = main_height / max(1, len(panel.row_ids))
    dense_rows = row_pitch * 72.0 < LAYOUT.dense_row_threshold_pt
    return Geometry(
        figure_height=figure_height,
        main_y=main_y,
        main_height=main_height,
        tree=tree,
        row_labels=row_labels,
        annotation_tracks=annotation_tracks,
        tracks=tracks,
        heatmap=heatmap,
        footer_blocks=footer_placements,
        footer_y=footer_y,
        footer_height=footer_height,
        row_pitch=row_pitch,
        dense_rows=dense_rows,
        vertical_scale=vertical_scale,
    )


def _add_axis(figure: plt.Figure,rectangle: tuple[float, float, float, float],) -> plt.Axes:
    x, y, w, h = rectangle
    width, height = figure.get_size_inches()
    return figure.add_axes((x / width, y / height, w / width, h / height))


def _add_dendrogram(
    ax: plt.Axes,
    linkage_matrix: np.ndarray | None,
    n_rows: int,
    rasterize: bool,
) -> None:
    if linkage_matrix is None or n_rows < 2:
        ax.axis("off")
        return
    hierarchy.dendrogram(
        linkage_matrix,
        orientation="left",
        no_labels=True,
        color_threshold=0,
        above_threshold_color=INK,
        link_color_func=lambda _: INK,
        ax=ax,
    )
    ax.set_ylim(10 * n_rows, 0)
    ax.margins(0)
    for collection in ax.collections:
        collection.set_linewidth(LAYOUT.dendrogram_lw)
        collection.set_rasterized(rasterize)
    ax.axis("off")


def _add_row_labels(ax: plt.Axes, labels: Sequence[str], n_rows: int) -> None:
    ax.set(xlim=(0, 1), ylim=(n_rows - 0.5, -0.5))
    ax.axis("off")
    for y, label in enumerate(labels):
        ax.text(0.96, y, label, ha="right", va="center", fontsize=STYLE.row_label_size)


def _add_annotation(ax: plt.Axes,track: AnnotationTrack,marker_size: float,dense_rows: bool,rasterized: bool = False,) -> None:
    n = len(track.values)
    marker_x = 0.0
    if dense_rows:
        verts = np.asarray([
            [
                (-0.5, y - 0.5),
                (0.5, y - 0.5),
                (0.5, y + 0.5),
                (-0.5, y + 0.5),
            ]
            for y in range(n)
        ], dtype=float)
        ax.add_collection(PolyCollection(
            verts,
            facecolors=track.colors,
            edgecolors="none",
            linewidths=0,
            antialiaseds=False,
            rasterized=rasterized,
        ))
    else:
        ax.scatter(
            np.full(n, marker_x), np.arange(n), s=marker_size, marker=track.marker,
            c=track.colors, clip_on=False, rasterized=rasterized, zorder=3,
        )
    ax.set(xlim=(-0.5, 0.5), ylim=(n - 0.5, -0.5), xticks=[], yticks=[])
    for spine in ax.spines.values():
        spine.set_visible(False)
    ax.text(
        marker_x, -0.018, track.label, transform=ax.get_xaxis_transform(), rotation=45,
        ha="right", va="top",fontsize=STYLE.annotation_size, clip_on=False,
    )


def _add_tracks(ax: plt.Axes,panel: PreparedHaplotypes,config: TVRHapPlotter,rasterized: bool,) -> None:
    plot_intervals(
        ax, panel.backgrounds,
        start="start", end="end", y="y", color="color",
        thickness=config.track_thickness, rasterized=rasterized, zorder=1,
    )
    plot_intervals(
        ax, panel.intervals,
        start="start", end="end", y="y", color="color",
        thickness=config.track_thickness, rasterized=rasterized, zorder=2,
    )
    ax.set(xlim=panel.xlim, ylim=(len(panel.row_ids) - 0.5, -0.5), yticks=[])
    format_kb_axis(ax, nbins=5, steps=(1, 2, 5, 10), min_n_ticks=3)
    ax.tick_params(axis="y", left=False)
    ax.set_xlabel("Distance to telomere-subtelomere boundary (kb)", labelpad=4)
    ax.spines[["top", "right"]].set_visible(False)
    ax.spines[["left", "bottom"]].set_color(INK)
    ax.spines[["left", "bottom"]].set_linewidth(STYLE.axis_linewidth)


def _add_heatmap(
    ax: plt.Axes,
    similarity: np.ndarray,
    cmap: str,
    rasterized: bool = False,
) -> QuadMesh:
    n = len(similarity)
    edges = np.arange(n + 1, dtype=float) - 0.5
    mesh = ax.pcolormesh(
        edges,
        edges,
        similarity,
        cmap=cmap,
        vmin=0,
        vmax=1,
        shading="flat",
        edgecolors="none",
        linewidth=0,
        antialiased=False,
        rasterized=rasterized,
    )
    ax.set(
        xlim=(-0.5, n - 0.5),
        ylim=(n - 0.5, -0.5),
        xticks=[],
        yticks=[],
    )
    for spine in ax.spines.values():
        spine.set_visible(False)
    return mesh


def _rasterize_heatmap(config: TVRHapPlotter, n_haplotypes: int) -> bool:
    """Resolve explicit rasterization or apply the large-heatmap default."""
    if config.rasterize is not None:
        return config.rasterize
    return n_haplotypes >= HEATMAP_RASTERIZE_THRESHOLD


def _annotation_handles(track: AnnotationTrack) -> list[Line2D]:
    return [
        Line2D(
            [0], [0], marker=track.marker, linestyle="None", color="none",
            markerfacecolor=color, markeredgecolor="none", markersize=4.5,
            label=value,
        )
        for value, color in _annotation_legend_entries(track)
    ]


def _motif_handles(panel: PreparedHaplotypes) -> list[Patch]:
    return [
        Patch(facecolor=panel.motif_colors.get(motif, OTHER_COLOR), edgecolor="none", label=motif)
        for motif in panel.motif_legend
    ]


def _footer_block_axis(
    figure: plt.Figure,
    geometry: Geometry,
    placement: FooterPlacement,
) -> plt.Axes:
    height = placement.height * geometry.vertical_scale
    y = geometry.footer_y + geometry.footer_height - height
    ax = _add_axis(figure, (placement.x, y, placement.width, height))
    ax.axis("off")
    return ax


def _style_discrete_legend(legend: mpl.legend.Legend) -> None:
    legend.get_title().set_fontweight("bold")
    legend.get_title().set_ha("left")


def _add_discrete_footer_block(
    ax: plt.Axes,
    handles: Sequence[mpl.artist.Artist],
    *,
    title: str,
    columns: int,
    motif: bool = False,
) -> None:
    legend = ax.legend(**_discrete_legend_kwargs(
        handles, title=title, columns=columns, motif=motif
    ))
    _style_discrete_legend(legend)

def _add_similarity_footer_block(figure: plt.Figure,ax: plt.Axes,heatmap_image: QuadMesh,) -> None:
    ax.text(
        0.0, 1.0, "TVR Similarity",
        transform=ax.transAxes,
        ha="left", va="top",
        fontsize=STYLE.legend_size,
        fontweight="bold",
    )

    bbox = ax.get_position()
    block_h = bbox.height * figure.get_figheight()
    title_h, tick_h, title_gap, tick_gap, bar_h = _similarity_legend_metrics()

    bar_top = block_h - title_h - title_gap
    bar_bottom = tick_h + tick_gap
    usable_bar_h = max(0.01, min(bar_h, bar_top - bar_bottom))
    cax = ax.inset_axes((0.0,bar_bottom / block_h,1.0,usable_bar_h / block_h,))
    colorbar = figure.colorbar(heatmap_image, cax=cax, orientation="horizontal")
    colorbar.set_ticks([0, 0.5, 1])
    colorbar.ax.tick_params(length=2.0,width=STYLE.tick_linewidth,pad=1.0,labelsize=STYLE.annotation_size,)
    colorbar.outline.set_linewidth(0)


def _add_legends(
    figure: plt.Figure,
    geometry: Geometry,
    panel: PreparedHaplotypes,
    heatmap_image: QuadMesh | None,
) -> None:
    placements = {placement.key: placement for placement in geometry.footer_blocks}

    motif_handles = _motif_handles(panel)
    motif_placement = placements.get("motif")
    if motif_handles and motif_placement is not None:
        ax = _footer_block_axis(figure, geometry, motif_placement)
        _add_discrete_footer_block(
            ax, motif_handles,
            title="TVR motifs",
            columns=motif_placement.columns,
            motif=True,
        )

    for index, track in enumerate(_annotation_legend_tracks(panel)):
        placement = placements.get(f"annotation:{index}")
        if placement is None:
            continue
        ax = _footer_block_axis(figure, geometry, placement)
        _add_discrete_footer_block(
            ax, _annotation_handles(track),
            title=track.label,
            columns=placement.columns,
        )

    similarity_placement = placements.get("similarity")
    if heatmap_image is not None and similarity_placement is not None:
        ax = _footer_block_axis(figure, geometry, similarity_placement)
        _add_similarity_footer_block(figure, ax, heatmap_image)


def _render_haplotypes(
    panel: PreparedHaplotypes,
    config: TVRHapPlotter,
    *,
    title: str,
) -> plt.Figure:
    geometry = figure_geometry(panel, config)
    figure = plt.figure(
        figsize=(config.width, geometry.figure_height),
        facecolor=BACKGROUND_COLOR,
    )
    main_rect = lambda item: (item[0], geometry.main_y, item[1], geometry.main_height)

    tree_ax = _add_axis(figure, main_rect(geometry.tree))
    label_ax = _add_axis(figure, main_rect(geometry.row_labels))
    annotation_axes = [
        _add_axis(figure, main_rect(item))
        for item in geometry.annotation_tracks
    ]
    track_ax = _add_axis(figure, main_rect(geometry.tracks))
    heatmap_ax = (
        _add_axis(figure, main_rect(geometry.heatmap))
        if geometry.heatmap is not None else None
    )

    n_rows = len(panel.row_ids)
    rasterized = config.rasterize is True
    _add_dendrogram(tree_ax, panel.linkage, n_rows, rasterized)
    _add_row_labels(label_ax, panel.row_labels, n_rows)
    marker_diameter = float(np.clip(geometry.row_pitch * 72.0 * 0.70, 1.0, 5.3))
    for axis, annotation in zip(annotation_axes, panel.annotations, strict=True):
        _add_annotation(axis, annotation, marker_diameter**2, geometry.dense_rows, rasterized)

    _add_tracks(track_ax, panel, config, rasterized)
    track_ax.set_title(title, fontweight="bold", pad=LAYOUT.title_pad)

    heatmap_image = None
    if heatmap_ax is not None and panel.similarity is not None:
        heatmap_image = _add_heatmap(
            heatmap_ax,
            panel.similarity,
            config.heatmap_cmap,
            _rasterize_heatmap(config, n_rows),
        )

    _add_legends(figure, geometry, panel, heatmap_image)
    return figure


def render_haplotypes(
    panel: PreparedHaplotypes,
    config: TVRHapPlotter,
    *,
    title: str,
) -> plt.Figure:
    with nature_rc():
        return _render_haplotypes(panel, config, title=title)


def _write_similarity(panel: PreparedHaplotypes, path: Path) -> None:

    if panel.similarity is None:
        return

    rows, cols = np.triu_indices_from(panel.similarity, k=1)

    result = pd.DataFrame({
        "tvr_hap1": [panel.row_ids[index] for index in rows],
        "tvr_hap2": [panel.row_ids[index] for index in cols],
        "similarity_score": panel.similarity[rows, cols],
    })

    result.to_csv(path,sep="\t",index=False,)

def plot_tvr_haps(
    consensus_file: str | Path | None = None,
    sample_sheet: str | Path | None = None,
    outdir: str | Path = ".",
    prefix: str = "sample",
    config: TVRHapPlotter | None = None,
    write_similarity_score: bool = False,
    logger=None,
) -> None:
    config = config or TVRHapPlotter()
    _validate_config(config)
    logger = logger or log_utils.get_logger()
    outdir = Path(outdir)
    outdir.mkdir(parents=True, exist_ok=True)

    frame = load_haplotypes(consensus_file, sample_sheet, prefix)
    panels = prepare_haplotypes(
        frame,
        config,
        write_similarity_score=write_similarity_score,
        logger=logger,
    )
    sample_label = config.sample_label or prefix

    if not frame.empty and frame["sample_id"].nunique(dropna=False) <= 1 and config.pool_chroms:
        logger.warning("--pool-chroms only applies to multi-sample plots; ignoring it for this input.")

    resolved_colors: dict[str, str] = {}
    for panel in panels.values():
        for motif, colour in panel.motif_colors.items():
            previous = resolved_colors.setdefault(motif, colour)
            if previous != colour:
                raise RuntimeError(f"Inconsistent TVR color for motif {motif!r} across panels.")
    if resolved_colors:
        color_map_path = outdir / f"{prefix}.tvr_colors.tsv"
        write_motif_color_map(
            color_map_path,
            resolved_colors,
            load_motif_color_overrides(config.tvr_colors),
        )
        logger.info(f"Writing TVR color map: {color_map_path}")
        displayed_color_count = sum(
            motif not in {"", "Other TVR motif"} for motif in resolved_colors
        )
        if displayed_color_count > 60:
            logger.warning(
                f"The figure uses {displayed_color_count} distinct TVR colors. "
                "Colors remain unique, but a legend of this size may be difficult "
                "to distinguish."
            )

    for output_stem, panel in panels.items():
        logger.info(f"Plotting {panel.title_suffix}")
        if write_similarity_score and panel.similarity is not None:
            _write_similarity(panel, outdir / f"{prefix}.{output_stem}.similarity.tsv")

        figure = render_haplotypes(
            panel,
            config,
            title=f"{sample_label} {panel.title_suffix} TVR haplotypes",
        )
        try:
            save_figure(
                figure,
                outdir / f"{prefix}.{output_stem}.hap_plot",
                output_format=config.output_format,
            )
        finally:
            plt.close(figure)
