# Copyright (C) 2025 Huihui Li <hhui.li@outlook.com>. Licensed under GNU GPL v3.0.
from __future__ import annotations

from contextlib import contextmanager
from dataclasses import dataclass
from pathlib import Path
from typing import Iterator, Mapping, Sequence, TypeVar
import re

import matplotlib
matplotlib.use("Agg")

import matplotlib.colors as mcolors
import matplotlib.pyplot as plt
import matplotlib.ticker as mticker
import numpy as np
import pandas as pd
from matplotlib.collections import PolyCollection
from natsort import natsort_keygen, natsorted


REQUIRED_COLUMNS = frozenset({"chrom", "chrom_phase", "arm"})
UNPHASED_VALUES = frozenset({"", "unphased", "all", "none", "nan", "na"})

BACKGROUND_COLOR = "#FFFFFF"
INK = "#202020"
OTHER_COLOR = "#CCCCCC"
MISSING_COLOR = "#E6E6E6"

PREFERRED_TVR_COLORS: dict[str, str] = {
    "TTAGGG": "#E5E5E5", "TTGGGG": "#4DBBD5", "TTAGGGG": "#00A087",
    "TGAGGG": "#1F77B4", "TCAGGG": "#FF9896", "TTCGGG": "#E64B35",
    "TAGGG": "#FFA90E", "TTAAGGG": "#9467BD", "CTAGGG": "#AEC7E8",
    "CGAGGG": "#F7B6D2", "TTTAGGG": "#2CA02C", "GTAGGG": "#FFBB78",
    "TTGGG": "#98DF8A", "TCGGGG": "#3C5488", "TTAGG": "#E7298A",
    "TGGGGG": "#BCBD22", "TTTGGG": "#91D1C2", "TAAGGG": "#9EA9C4",
    "TTCGGGG": "#E377C2", "TTAGAGGG": "#8C564B", "TTAGGGGG": "#92DCE5",
    "TGAGGGG": "#D4AF37", "TGAGG": "#56B4E9", "TTTTAGGG": "#B09C85",
    "TCAGGGG": "#C5B0D5", "TTACGG": "#A96B59", "CTGGGG": "#D8D0C2",
    "TGGGG": "#DBDB8D", "TTAGGGTTG": "#005E59", "TTAGGGTGG": "#5BA396",
    "TTAGGGTG": "#B4D8D1", "TTGGGGG": "#5B6827", "TTAAA": "#CD5C5C",
    "TTAAAA": "#B03060", "TTAAAAA": "#7A2143", "TTAAAAAA": "#619CFF",
    "GAAGAA": "#6A3D9A", "AAGAAG": "#6A3D9A",
    "TG": "#4DBBD5", "TGG": "#92DCE5", "TGGG": "#1F77B4",
    "GG": "#E64B35", "AGGG": "#9467BD", "AG": "#FFA90E",
    "TGGA": "#2CA02C", "TT": "#E377C2", "TTG": "#98DF8A",
    "TGAA": "#FF9896", "AAG": "#AEC7E8", "TAG": "#BCBD22",
    "TA": "#E7298A",
}

TVR_COLOR_POOL = (
    "#4DBBD5", "#00A087", "#1F77B4", "#FF9896", "#E64B35",
    "#FFA90E", "#9467BD", "#AEC7E8", "#F7B6D2", "#2CA02C", "#FFBB78",
    "#98DF8A", "#3C5488", "#E7298A", "#BCBD22", "#91D1C2", "#9EA9C4",
    "#E377C2", "#8C564B", "#92DCE5", "#D4AF37", "#56B4E9", "#B09C85",
    "#C5B0D5", "#A96B59", "#D8D0C2", "#DBDB8D", "#005E59", "#5BA396",
    "#B4D8D1", "#5B6827", "#CD5C5C", "#B03060", "#7A2143", "#619CFF",
    "#6A3D9A", "#9E0000", "#AE00AE", "#6551FF", "#A27900", "#DB39FF",
    "#35D214", "#AE92FF", "#714500", "#865179", "#E37504", "#86AA6D",
    "#007D49", "#AE4900", "#0431FF", "#00C682", "#00F7BA", "#AAEB00",
    "#417D86", "#9A00FF", "#F7A2FF", "#7982AE", "#8A9A00", "#AE7192",
    "#146100", "#FF658A", "#CE0035", "#00FB6D", "#E300C6", "#6108E3",
    "#DB9265", "#7D8649", "#82351C", "#C27DDB", "#7D2882", "#866531",
    "#FBCE00", "#6596AE", "#2D79FF", "#B659B2", "#9269FF", "#696192",
    "#2D5D39", "#497D00", "#9A0071", "#00C6B6", "#FF754D", "#006582",
    "#599265", "#FF61FF", "#CE868E", "#AA9A4D", "#AA3D39", "#9EB2FF",
    "#00F3E7", "#F7285D", "#82BE49", "#9214BE", "#1496D7", "#A28EC2",
    "#9A0031", "#BA793D", "#C2108E", "#B6BA75", "#5D6DBE", "#C62800",
    "#FF8EBA", "#B24DFB", "#00DFFF", "#FF45BE", "#7555C2", "#EBC6FF",
    "#C692BE", "#CA8E18", "#A25569", "#D2DF51", "#009AA2", "#82C296",
    "#007969", "#D26192", "#CEAAFF", "#245DC2", "#797500", "#BA00E7",
    "#654571", "#FF923D", "#868AE3", "#7D414D", "#8A4D9E", "#CA6131",
    "#24BE4D", "#556D92",
)
ANNOTATION_PALETTE = (
    "#4DBBD5", "#00A087", "#3C5488", "#F39B7F", "#E64B35", "#F4A300",
    "#8491B4", "#91D1C2", "#DC0000", "#7E6148", "#B09C85", "#BCBD22",
    "#2CA02C", "#98DF8A", "#9467BD", "#C5B0D5", "#1F77B4", "#AEC7E8",
    "#E377C2", "#F7B6D2", "#A96B59", "#B03060", "#005E59", "#9EA9C4",
)
SPECIAL_ANNOTATION_COLORS: dict[str, dict[str, str]] = {
    "chrom_phase": {
        "mat": "#2C6DB2", "maternal": "#2C6DB2", "hap1": "#2C6DB2",
        "pat": "#F06449", "paternal": "#F06449", "hap2": "#F06449",
    },
    "arm": {"p": "#B8A7D1", "l": "#B8A7D1", "q": "#8FC7B5", "r": "#8FC7B5"},
}
_CHROM_PREFIX_RE = re.compile(
    r"^(?:chromosome|chromsome|chrom|chr)"
    r"(?=$|[\s_.:-]|\d|x|y|m|un|[ivxlcdm]+(?:$|[\s_.:-]))[\s_.:-]*",
    re.IGNORECASE,
)
T = TypeVar("T")
_NATURAL_SORT_KEY = natsort_keygen()
_ARM_ORDER = {"p": 0, "q": 1, "l": 2, "r": 3}


@dataclass(frozen=True, slots=True)
class NatureStyle:
    single_column_mm: float = 89.0
    double_column_mm: float = 183.0
    font_family: tuple[str, ...] = ("Arial", "Helvetica", "Liberation Sans", "DejaVu Sans", "sans-serif")
    label_size: float = 7.0
    title_size: float = 7.0
    tick_size: float = 6.0
    legend_size: float = 6.0
    annotation_size: float = 6.0
    row_label_size: float = 5.8
    panel_label_size: float = 8.0
    axis_linewidth: float = 0.6
    tick_linewidth: float = 0.6
    major_tick_length: float = 3.2
    minor_tick_length: float = 1.8
    data_linewidth: float = 1.0
    separator_linewidth: float = 0.4
    raster_dpi: int = 600

    @staticmethod
    def inches(mm: float) -> float:
        return mm / 25.4


STYLE = NatureStyle()
FIGURE_FORMATS = ("pdf", "png", "svg")


@contextmanager
def nature_rc(style: NatureStyle = STYLE) -> Iterator[None]:
    with plt.rc_context({
        "font.family": "sans-serif",
        "font.sans-serif": list(style.font_family),
        "font.size": style.label_size,
        "text.color": INK,
        "axes.labelsize": style.label_size,
        "axes.labelcolor": INK,
        "axes.titlesize": style.title_size,
        "axes.edgecolor": INK,
        "axes.linewidth": style.axis_linewidth,
        "xtick.labelsize": style.tick_size,
        "ytick.labelsize": style.tick_size,
        "xtick.color": INK,
        "ytick.color": INK,
        "xtick.major.size": style.major_tick_length,
        "ytick.major.size": style.major_tick_length,
        "xtick.major.width": style.tick_linewidth,
        "ytick.major.width": style.tick_linewidth,
        "xtick.minor.size": style.minor_tick_length,
        "ytick.minor.size": style.minor_tick_length,
        "xtick.minor.width": style.tick_linewidth,
        "ytick.minor.width": style.tick_linewidth,
        "xtick.direction": "out",
        "ytick.direction": "out",
        "legend.fontsize": style.legend_size,
        "legend.frameon": False,
        "axes.grid": False,
        "figure.facecolor": BACKGROUND_COLOR,
        "axes.facecolor": BACKGROUND_COLOR,
        "savefig.facecolor": BACKGROUND_COLOR,
        "pdf.fonttype": 42,
        "ps.fonttype": 42,
        "svg.fonttype": "none",
    }):
        yield


def split_csv(value: str | Sequence[str] | None) -> tuple[str, ...]:
    if value is None:
        return ()
    if isinstance(value, str):
        return tuple(part.strip() for part in value.split(",") if part.strip())
    return tuple(str(part).strip() for part in value if str(part).strip())


def safe_name(value: object) -> str:
    return re.sub(r"[^A-Za-z0-9._-]+", "_", str(value).strip()).strip("._") or "value"


def display_chrom(value: object) -> str:
    return _CHROM_PREFIX_RE.sub("", str(value).strip(), count=1)


def resolve_requested_values(
    value: str | Sequence[str] | None,
    observed: Sequence[object],
    *,
    label: str = "value",
) -> tuple[str, ...]:
    available = tuple(dict.fromkeys(str(item) for item in observed))
    wanted = split_csv(value)
    if not wanted or (len(wanted) == 1 and wanted[0].lower() == "all"):
        return available

    available_set = set(available)
    selected = tuple(dict.fromkeys(item for item in wanted if item in available_set))
    missing = tuple(item for item in wanted if item not in available_set)
    if missing:
        raise ValueError(
            f"Unknown {label}(s): {', '.join(missing)}; "
            f"observed: {', '.join(available)}"
        )
    return selected


def require_schema(
    frame: pd.DataFrame,
    required: Sequence[str] = REQUIRED_COLUMNS,
    *,
    source: str | Path | None = None,
) -> None:
    missing = sorted(set(required).difference(frame.columns))
    if missing:
        if source is not None:
            raise ValueError(f"{source}: missing required column(s): {', '.join(missing)}")
        raise ValueError(
            f"Missing required column: {missing[0]}. "
            "Regenerate this file using the current Telox pipeline."
        )


def normalize_chrom_phase(value: object) -> str:
    if pd.isna(value):
        return "unphased"
    phase = str(value).strip()
    return "unphased" if phase.lower() in UNPHASED_VALUES else phase


def chrom_phase_key(value: object) -> str:
    phase = normalize_chrom_phase(value)
    return "" if phase == "unphased" else phase


def chrom_phase_label(value: object) -> str:
    return chrom_phase_key(value).upper()


def arm_group_key(value: object) -> str:
    text = str(value).strip()
    if text.lower() in {"p", "l"}:
        return "left"
    if text.lower() in {"q", "r"}:
        return "right"
    return f"other:{text}"


def arm_sort_key(value: object) -> tuple[int, object]:
    label = str(value).strip()
    return _ARM_ORDER.get(label.casefold(), len(_ARM_ORDER)), _NATURAL_SORT_KEY(label)


def row_major(handles: Sequence[T], columns: int) -> list[T]:
    if columns <= 1:
        return list(handles)
    rows = (len(handles) + columns - 1) // columns
    return [
        handles[index]
        for column in range(columns)
        for row in range(rows)
        if (index := row * columns + column) < len(handles)
    ]


def plot_intervals(
    ax,
    frame: pd.DataFrame,
    *,
    start: str,
    end: str,
    y: str,
    color: str,
    thickness: float | str = 1.0,
    rasterized: bool = False,
    zorder: int = 1,
    alpha: float = 1.0,
) -> PolyCollection | None:
    if frame.empty:
        return None
    x0, x1, ys = (frame[column].to_numpy(float) for column in (start, end, y))
    heights = (
        frame[thickness].to_numpy(float)
        if isinstance(thickness, str)
        else np.full(len(frame), thickness)
    )
    verts = np.stack((
        np.column_stack((x0, ys - heights / 2)),
        np.column_stack((x1, ys - heights / 2)),
        np.column_stack((x1, ys + heights / 2)),
        np.column_stack((x0, ys + heights / 2)),
    ), axis=1)
    collection = PolyCollection(
        verts,
        facecolors=frame[color].fillna(MISSING_COLOR),
        edgecolors="none",
        linewidths=0,
        rasterized=rasterized,
        zorder=zorder,
        alpha=alpha,
    )
    ax.add_collection(collection)
    return collection


def normalize_figure_format(value: object) -> str:
    figure_format = str(value).strip().lower().lstrip(".")
    if figure_format not in FIGURE_FORMATS:
        supported = ", ".join(FIGURE_FORMATS)
        raise ValueError(
            f"Unsupported figure format: {value!r}. Choose one of: {supported}."
        )
    return figure_format


def save_figure(
    fig,
    path: str | Path,
    *,
    output_format: str | None = None,
    style: NatureStyle = STYLE,
) -> Path:
    path = Path(path)
    selected_format = output_format if output_format is not None else (path.suffix or "pdf")
    figure_format = normalize_figure_format(selected_format)
    existing_format = path.suffix.lower().lstrip(".")
    if existing_format in FIGURE_FORMATS:
        path = path.with_suffix(f".{figure_format}")
    else:
        path = path.with_name(f"{path.name}.{figure_format}")
    with nature_rc(style):
        fig.savefig(
            path,
            format=figure_format,
            dpi=style.raster_dpi,
            metadata={"Creator": "Teloxplorer", "Title": path.stem},
        )
    return path


def kb_formatter(value: float, _position: int | None = None) -> str:
    value /= 1000.0
    return f"{value:.0f}" if np.isclose(value, round(value)) else f"{value:.1f}"


def format_kb_axis(
    ax,
    *,
    nbins: int = 6,
    steps: Sequence[float] = (1, 2, 2.5, 5, 10),
    min_n_ticks: int | None = None,
    minor_intervals: int = 2,
    major_pad: float = 2,
) -> None:
    locator_kwargs = {"nbins": nbins, "steps": steps}
    if min_n_ticks is not None:
        locator_kwargs["min_n_ticks"] = min_n_ticks
    ax.xaxis.set_major_locator(mticker.MaxNLocator(**locator_kwargs))
    ax.xaxis.set_minor_locator(mticker.AutoMinorLocator(minor_intervals))
    ax.xaxis.set_major_formatter(mticker.FuncFormatter(kb_formatter))
    ax.tick_params(axis="x", which="major", pad=major_pad)


def parse_motif_blocks(blocks: object) -> tuple[tuple[str, int], ...]:
    if pd.isna(blocks) or str(blocks).strip() in {"", "."}:
        return ()
    parsed: list[tuple[str, int]] = []
    for token in str(blocks).split(","):
        motif, separator, count_text = token.strip().rpartition(":")
        if not separator or not motif:
            raise ValueError(f"Invalid motif block: {token!r}")
        try:
            count = int(count_text)
        except ValueError as exc:
            raise ValueError(f"Invalid motif count in block: {token!r}") from exc
        if count <= 0:
            raise ValueError(f"Motif block count must be positive: {token!r}")
        parsed.append((motif, count))
    return tuple(parsed)


def load_color_overrides(
    path: str | Path | None,
    *,
    key_name: str = "VALUE",
) -> dict[str, str]:
    if path is None:
        return {}
    overrides: dict[str, str] = {}
    for line_number, line in enumerate(Path(path).read_text().splitlines(), start=1):
        if not line.strip() or line.lstrip().startswith("#"):
            continue
        fields = [field.strip().strip("\"'") for field in line.split("\t")]
        if len(fields) < 2 or not fields[0] or not mcolors.is_color_like(fields[1]):
            raise ValueError(f"Invalid colour on line {line_number}: expected {key_name}<TAB>COLOR.")
        overrides[fields[0]] = fields[1]
    return overrides


def load_motif_color_overrides(path: str | Path | None) -> dict[str, str]:
    overrides = load_color_overrides(path, key_name="MOTIF")
    return _normalize_unique_colors(overrides, label="TVR color override")


def _normalize_color(value: str) -> str:
    rgba = mcolors.to_rgba(value)
    return mcolors.to_hex(rgba, keep_alpha=not np.isclose(rgba[3], 1.0)).upper()


def _normalize_unique_colors(
    colors: Mapping[str, str],
    *,
    label: str,
) -> dict[str, str]:
    normalized = {str(key): _normalize_color(value) for key, value in colors.items()}
    by_color: dict[str, list[str]] = {}
    for key, colour in normalized.items():
        by_color.setdefault(colour, []).append(key)
    duplicates = {colour: keys for colour, keys in by_color.items() if len(keys) > 1}
    if duplicates:
        detail = "; ".join(
            f"{colour}: {', '.join(keys)}" for colour, keys in duplicates.items()
        )
        raise ValueError(f"{label}s must be unique; duplicate assignments: {detail}.")
    return normalized


def motif_color_map(
    abundance: Mapping[str, float] | pd.Series,
    top_n: int,
    overrides: Mapping[str, str] | None = None,
    *,
    exclude: Sequence[str] = (),
    reserved_colors: Sequence[str] = (),
    priority_motifs: Sequence[str] = (),
) -> tuple[dict[str, str], tuple[str, ...]]:
    if top_n < 0:
        raise ValueError("top_n must be non-negative.")
    excluded = set(exclude)
    totals = {str(motif): float(value) for motif, value in dict(abundance).items()}
    ordered = tuple(
        sorted(
            (motif for motif in totals if motif not in excluded),
            key=lambda motif: (-totals[motif], motif),
        )[:top_n]
    )

    overrides = _normalize_unique_colors(overrides or {}, label="TVR color override")
    used = {
        _normalize_color(colour)
        for colour in (*reserved_colors, OTHER_COLOR, MISSING_COLOR, BACKGROUND_COLOR)
    }
    assigned: dict[str, str] = {}
    prioritized = tuple(
        motif for motif in dict.fromkeys((*priority_motifs, *ordered))
        if motif in ordered
    )

    for motif in prioritized:
        if motif not in overrides:
            continue
        colour = overrides[motif]
        if colour in used:
            raise ValueError(
                f"TVR color override for {motif!r} uses reserved color {colour}. "
                "Assign a unique color for every displayed motif."
            )
        assigned[motif] = colour
        used.add(colour)

    unresolved: list[str] = []
    for motif in prioritized:
        if motif in assigned:
            continue
        preferred = PREFERRED_TVR_COLORS.get(motif)
        if preferred is not None and _normalize_color(preferred) not in used:
            colour = _normalize_color(preferred)
            assigned[motif] = colour
            used.add(colour)
        else:
            unresolved.append(motif)

    available = (colour for colour in TVR_COLOR_POOL if colour not in used)
    for motif in unresolved:
        try:
            colour = next(available)
        except StopIteration as exc:
            raise ValueError(
                f"Need {len(ordered)} unique TVR colors, but the built-in pool "
                f"supports at most {len(TVR_COLOR_POOL)}. Provide additional "
                "unique colors with --tvr-colors or reduce --top-tvr."
            ) from exc
        assigned[motif] = colour
        used.add(colour)

    colors = {motif: assigned[motif] for motif in ordered}
    if len(set(colors.values())) != len(colors):
        raise RuntimeError("Internal error: displayed TVR motif colors are not unique.")
    return colors, ordered


def motif_color(motif: str, overrides: Mapping[str, str] | None = None) -> str:
    overrides = overrides or {}
    colour = overrides.get(motif) or PREFERRED_TVR_COLORS.get(motif) or MISSING_COLOR
    return _normalize_color(colour)


def write_motif_color_map(
    path: str | Path,
    colors: Mapping[str, str],
    overrides: Mapping[str, str] | None = None,
) -> None:
    overrides = _normalize_unique_colors(
        overrides or {}, label="TVR color override"
    )
    base_colors = frozenset(PREFERRED_TVR_COLORS.values())
    lines = ["#MOTIF\tCOLOR\tSOURCE"]
    for motif, value in colors.items():
        if not motif or motif == "Other TVR motif":
            continue
        if any(separator in motif for separator in ("\t", "\n", "\r")):
            raise ValueError(f"Invalid motif for color-map output: {motif!r}.")
        colour = _normalize_color(value)
        if motif in overrides:
            source = "override"
        elif PREFERRED_TVR_COLORS.get(motif) == colour:
            source = "preferred"
        elif colour in base_colors:
            source = "base-pool"
        else:
            source = "extended-pool"
        lines.append(f"{motif}\t{colour}\t{source}")
    Path(path).write_text("\n".join(lines) + "\n")


def categorical_color_map(
    values: Sequence[object],
    *,
    column: str | None = None,
    overrides: Mapping[str, str] | None = None,
) -> dict[str, str]:
    overrides = overrides or {}
    unique = natsorted({str(value) for value in values})
    preferred = SPECIAL_ANNOTATION_COLORS.get(column or "", {})
    unresolved = [
        value for value in unique
        if value != "NA" and value not in overrides and value.lower() not in preferred
    ]
    fallback = {
        value: ANNOTATION_PALETTE[index % len(ANNOTATION_PALETTE)]
        for index, value in enumerate(unresolved)
    }
    palette = {"NA": MISSING_COLOR}
    for value in unique:
        if value != "NA":
            palette[value] = overrides.get(
                value,
                preferred.get(value.lower(), fallback.get(value, MISSING_COLOR)),
            )
    return palette


def parse_sample_sheet(
    *,
    primary_file: str | Path | None,
    sample_sheet: str | Path | None,
    prefix: str,
) -> dict[str, Path]:
    if sample_sheet is None:
        if primary_file is None:
            raise ValueError("Provide an input file or a sample sheet.")
        return {prefix: Path(primary_file)}

    rows: dict[str, Path] = {}
    for line_number, line in enumerate(Path(sample_sheet).read_text().splitlines(), start=1):
        if not line.strip() or line.lstrip().startswith("#"):
            continue
        fields = line.split()
        if len(fields) < 2:
            raise ValueError(f"Invalid sample-sheet line {line_number}: expected SAMPLE<TAB>PATH.")
        name, path = fields[:2]
        if name in rows:
            raise ValueError(f"Duplicate sample name in sample sheet: {name!r}")
        rows[name] = Path(path)
    if not rows:
        raise ValueError("Sample sheet contains no samples.")
    return rows
