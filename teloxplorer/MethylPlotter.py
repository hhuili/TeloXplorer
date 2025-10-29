import argparse
import os
import logging
import pandas as pd
import numpy as np
import matplotlib
matplotlib.use('Agg')
from natsort import natsorted
from plotnine import *
from plotnine.themes import theme
from .utilities import setup_basic_logging,log_parameter_summary


def smooth(x, window_len=8, window='hanning'):
    """
    Smoothing functions adapted from concepts in methylartist (citation needed).
    Original algorithm: sliding window averaging with convolution-based smoothing.
    We reimplemented the core smoothing logic for our specific use case.
    """
    if x.size < window_len:
        return x

    if window_len < 3:
        return x
    
    if window_len % 2 != 0:
        window_len += 1

    window_functions = {
        'flat': lambda size: np.ones(size, 'd'),
        'hanning': np.hanning,
        'hamming': np.hamming,
        'bartlett': np.bartlett,
        'blackman': np.blackman
    }

    s = np.pad(x, (window_len-1, window_len-1), mode='reflect')

    w = window_functions.get(window, np.hanning)(window_len)
    w = w / w.sum()

    y = np.convolve(w, s, mode='valid')
    
    half_window = int(window_len/2)

    return y[half_window-1:-(half_window)]

def slide_window(positions, values, width=1000, slide=100):

    if len(positions) == 0:
        return {}, {}
    
    midpt_min = np.min(positions)
    midpt_max = np.max(positions)
    
    win_start = int(midpt_min - width/2)
    win_end = win_start + width
    
    methyl_frac = {}
    
    while int((win_start+win_end)/2) < midpt_max:
        win_start += slide
        win_end += slide

        window_indices = (positions > win_start) & (positions < win_end)
        window_values = values[window_indices]
        
        midpt = int((win_start+win_end)/2)
        
        if len(window_values) > 0:
            methyl_frac[midpt] = np.mean(window_values)
    
    return methyl_frac


def calculate_smoothed_data(data2plot, window_size=1000, step_size=100, smooth_window=20):

    smoothed_data_list = []
    
    for (chr_arm, id_group), group_data in data2plot.groupby(['chr_arm', 'id']):

        clean_data = group_data.dropna(subset=['start_position', 'percent_modified'])
        
        if len(clean_data) == 0:
            continue

        positions = clean_data['start_position'].values
        values = clean_data['percent_modified'].values

        sort_idx = np.argsort(positions)
        positions_sorted = positions[sort_idx]
        values_sorted = values[sort_idx]

        windowed_data = slide_window(positions_sorted, values_sorted, width=window_size, slide=step_size)
        
        if len(windowed_data) > 0:
            window_positions = np.array(list(windowed_data.keys()))
            window_values = np.array(list(windowed_data.values()))

            if len(window_values) >= smooth_window:
                smoothed_values = smooth(window_values, window_len=smooth_window, window='hanning')
            else:
                smoothed_values = window_values

            for pos, val in zip(window_positions, smoothed_values):
                smoothed_data_list.append({
                    'chr_arm': chr_arm,
                    'id': id_group,
                    'start_position': pos,
                    'percent_modified': val,
                    'data_type': 'smoothed'
                })
    
    return pd.DataFrame(smoothed_data_list)


def plot_ggline_wo_dots(data2plot, smoothed_data, region_size, hide_raw=False):

    unique_ids = data2plot['id'].nunique()

    p = (
        ggplot(data2plot) +
        aes(x="start_position", y="percent_modified")
    )

    if not hide_raw:
        p += geom_line(aes(color="id"), size=0.2, alpha=0.3)

    p += geom_line(
        aes(x="start_position", y="percent_modified", color="id"),
        data=smoothed_data,
        size=0.8
    )
    
    p += scale_color_manual(values=methyline_colors)

    p += facet_wrap(
        "~chr_arm", 
        ncol=6
    )
    
    p += ylim(0, 100)
    p += scale_x_continuous(
        labels=lambda x: [f'{val/1000:.0f}' for val in x],
        limits=(-region_size, region_size)
    )
    p += labs(
        x=f"Coordinate (Unit: kb)",
        y="DNA 5mC methylation level"
    )
    p += theme_matplotlib()
    p += theme(
        panel_border=element_rect(color='black'),
        panel_grid=element_blank(),
        axis_text_x=element_text(size=8, angle=45, vjust=1, hjust=0.5),
        axis_text_y=element_text(size=8),
        strip_text=element_text(size=8),
        legend_title=element_blank(),
        legend_position="bottom",
        aspect_ratio=0.5
    )

    if unique_ids <= 1:
        p += theme(legend_position="none")

    return p

def run_methyl_plot(bedmethyl,boundary,outdir,out_prefix,args):

    #### Setup colors
    global methyline_colors

    set1_colors = ['#E41A1C', '#377EB8', '#4DAF4A', '#984EA3', '#FF7F00', '#FFFF33', '#A65628', '#F781BF', '#999999']
    set2_colors = ['#66C2A5', '#FC8D62', '#8DA0CB', '#E78AC3', '#A6D854', '#FFD92F', '#E5C494', '#B3B3B3']
    dark2_colors = ['#1B9E77', '#D95F02', '#7570B3', '#E7298A', '#66A61E', '#E6AB02', '#A6761D', '#666666']

    methyline_colors = set1_colors[:8] + set2_colors[:6] + dark2_colors[:4]

    methybed_header = [
        "chrom", "start_position", "end_position", "modified_base_code_and_motif",
        "score", "strand", "start_position1", "end_position1", "color",
        "N_valid_cov", "percent_modified", "N_mod", "N_canonical",
        "N_other_mod", "N_delete", "N_fail", "N_diff", "N_nocall"
    ]

    #### LOAD
    dt_bedmethyl = pd.read_csv(bedmethyl, sep="\t", header=None, names=methybed_header)

    exlude_chrs = args.skip_chr.split(',')
    
    # Get unique chromosomes and sort them naturally
    chr_level = dt_bedmethyl["chrom"].str.replace("_MATERNAL|_PATERNAL", "", regex=True).unique()
    
    # Create chr_arm levels with natural sorting
    chr_arm_level = [f"{c}{arm}" for c in chr_level for arm in ["p", "q"]]
    chr_arm_level = [c for c in chr_arm_level if c not in exlude_chrs]

    df_boundary = pd.read_csv(boundary, sep="\t")
    df_boundary["arm"] = np.where(df_boundary["tel_arm"] == "L", "p", "q")
    df_boundary["chrom"] = df_boundary["chrom"].replace({"MT": "M", "Mito": "M"})
    df_boundary["chr_arm"] = df_boundary["chrom"] + df_boundary["arm"]
    df_boundary = df_boundary[["chr_arm", "chrom", "arm", "tel_length", "tel_boundary"]].drop_duplicates()

    df_boundary_wider = df_boundary.pivot(index="chrom", columns="arm", values="tel_boundary").reset_index()
    df_boundary_wider.columns.name = None

    #### TIDY
    dt_bedmethyl_lite = dt_bedmethyl[dt_bedmethyl["score"] >= args.min_valid_cov][["chrom", "start_position", "percent_modified"]]

    dt_bedmethyl_separm = pd.merge(dt_bedmethyl_lite, df_boundary_wider, on="chrom", how="left")

    # Add 'chr' prefix if not present
    dt_bedmethyl_separm["chrom"] = dt_bedmethyl_separm["chrom"].apply(lambda x: f"chr{x}" if not x.startswith("chr") else x)
    dt_bedmethyl_separm["chrom"] = dt_bedmethyl_separm["chrom"].replace({"chrMito": "chrM"})

    # Process phasing
    chrom_split = dt_bedmethyl_separm["chrom"].str.partition("_")
    dt_bedmethyl_separm["chrom"] = chrom_split[0]
    dt_bedmethyl_separm["haplotype"] = chrom_split[2]

    dt_bedmethyl_separm["id"] = "Methyl"
    dt_bedmethyl_separm["id"] = dt_bedmethyl_separm.apply(lambda row: f"{row['id']} {row['haplotype']}" if pd.notna(row['haplotype']) else row['id'], axis=1)

    # Allocate chrom arm
    conditions = [
        dt_bedmethyl_separm["start_position"] <= dt_bedmethyl_separm["p"] + args.plot_flank,
        dt_bedmethyl_separm["start_position"] >= dt_bedmethyl_separm["q"] - args.plot_flank
    ]
    choices = ["p", "q"]
    dt_bedmethyl_separm["arm"] = np.select(conditions, choices, default=None)

    dt_bedmethyl_separm = dt_bedmethyl_separm.dropna(subset=["arm"])
    dt_bedmethyl_separm["tel_boundary"] = np.where(dt_bedmethyl_separm["arm"] == "p", dt_bedmethyl_separm["p"], dt_bedmethyl_separm["q"])
    dt_bedmethyl_separm["chr_arm"] = dt_bedmethyl_separm["chrom"] + dt_bedmethyl_separm["arm"]

    # Convert to relative coordinates
    dt_bedmethyl_relative = dt_bedmethyl_separm.copy()
    dt_bedmethyl_relative["start_position"] = dt_bedmethyl_relative["start_position"] - dt_bedmethyl_relative["tel_boundary"]

    # Adjust chromosomal arm tags
    dt_bedmethyl_tmp = dt_bedmethyl_relative.copy()
    if args.arm_label == "LR":
        dt_bedmethyl_tmp["arm"] = dt_bedmethyl_tmp["arm"].replace({"p": "L", "q": "R"})
        dt_bedmethyl_tmp["chr_arm"] = dt_bedmethyl_tmp["chrom"] + "_" + dt_bedmethyl_tmp["arm"]
        chr_level = [c.replace("p", "_L").replace("q", "_R") for c in chr_arm_level]
    else:
        chr_level = chr_arm_level

    # Sort chromosomes with natural order
    chr_level_sorted = natsorted(chr_level)
    dt_bedmethyl_tmp["chr_arm"] = pd.Categorical(dt_bedmethyl_tmp["chr_arm"], categories=chr_level_sorted, ordered=True)

    # Complete data for plotting
    all_combs = pd.MultiIndex.from_product([dt_bedmethyl_tmp['id'].unique(), dt_bedmethyl_tmp['chr_arm'].cat.categories], names=['id', 'chr_arm']).to_frame(index=False)
    dt_bedmethyl_2plot = pd.merge(all_combs, dt_bedmethyl_tmp, on=['id', 'chr_arm'], how='left')

    # Calculate smoothed data using methylartist method
    smoothed_data = calculate_smoothed_data(
        dt_bedmethyl_2plot, 
        window_size=args.window_size,
        step_size=args.step_size, 
        smooth_window=args.smooth_window
    )

    # Set chromosomes factor level
    dt_bedmethyl_2plot['chr_arm'] = pd.Categorical(dt_bedmethyl_2plot['chr_arm'], categories=chr_level_sorted, ordered=True)
    smoothed_data['chr_arm'] = pd.Categorical(smoothed_data['chr_arm'],categories=chr_level_sorted,ordered=True)

    #### output data
    raw_data_output= os.path.join(outdir, f"{out_prefix}.raw_data_for_plot.tsv")
    smoothed_data_output = os.path.join(outdir, f"{out_prefix}.smoothed_data_for_plot.tsv")

    #### writing plot data
    logging.info(f"> Writing raw data for plotting to: {raw_data_output}")
    dt_bedmethyl_2plot.to_csv(raw_data_output, sep='\t', index=False)

    logging.info(f"> Writing smoothed data for plotting to: {smoothed_data_output}")
    smoothed_data.to_csv(smoothed_data_output, sep='\t', index=False)

    #### PLOT
    matplotlib.set_loglevel("error")
    matplotlib.rcParams['pdf.fonttype'] = 42
    matplotlib.rcParams['font.family'] = 'Arial'
    matplotlib.rcParams['font.size'] = 10
    logging.getLogger('fontTools').setLevel(logging.WARNING)
    logging.getLogger('matplotlib').setLevel(logging.WARNING)

    p = plot_ggline_wo_dots(
        dt_bedmethyl_2plot,
        smoothed_data,
        args.plot_flank,
        set_theme="bw",
        hide_raw=args.hide_raw
    )
    p += geom_vline(xintercept=0, linetype="dashed",)

    output = os.path.join(outdir, f"{out_prefix}.tel_methyl.pdf")
    p.save(output, width=args.width, height=args.height, dpi=300, verbose=False)

def setup_common_args(parser):

    group = parser.add_argument_group("Input/Output Options")
    group.add_argument("--bed-methyl",
                       required=True,
                       help="bedMethyl file from ModKit")
    group.add_argument("--boundary",
                       required=True,
                       help="Telomere-subtelomere boundary file from telox-asm")
    group.add_argument(f"--out-prefix",
                       default="plot",
                       help="Output prefix [default: %(default)s].")
    group.add_argument(f"--outdir",
                       default=".",
                       help="Output directory [default: %(default)s].")

    return parser

def setup_plotting_args(parser):

    group_name = 'Plotting Options'
    group = parser.add_argument_group(group_name, 'DNA methylation plotting parameters')

    group.add_argument("--min-valid-cov",
                       type=int,
                       default=0,
                       help="Minimum valid coverage (Nmod + Nother_mod + Ncanonical), also used as the score in the bedMethyl [default: %(default)s]")
    group.add_argument("--plot-flank",
                       type=int,
                       default=20000,
                       help="Flanking size around boundary for plotting (bp) [default: %(default)s]")
    group.add_argument("--window-size",
                       type=int,
                       default=1000,
                       help="Sliding window size (bp) [default: %(default)s]")
    group.add_argument("--step-size",
                       type=int,
                       default=100,
                       help="Sliding step size (bp) [default: %(default)s]")
    group.add_argument("--smooth-window",
                       type=int,
                       default=20, 
                       help="Smoothing window size [default: %(default)s]")
    group.add_argument("--arm-label",
                       choices=["pq", "LR"],
                       default="pq",
                       help="Chromosome arm labeling scheme [default: %(default)s]")
    group.add_argument("--skip-chr",
                       default="chrMp,chrMq",
                       help="Chromosomes to exclude [default: %(default)s]")
    group.add_argument("--hide-raw",
                       action="store_true",
                       help="Hide raw methylation data in plot")
    group.add_argument("-W", "--width",
                       type=int,
                       default=8,
                       help="Plot width (inches) [default: %(default)s]")
    group.add_argument("-H", "--height",
                       type=int,
                       default=8,
                       help="Plot height (inches) [default: %(default)s]")

    return parser


def parse_arguments():

    parser = argparse.ArgumentParser(
        description="DNA methylation plotting parameters."
    )

    parser = setup_common_args(parser)
    parser = setup_plotting_args(parser)

    return parser.parse_args()

def methyl_plot_cli():

    module = "methyl-plot"
    args = parse_arguments()

    os.makedirs(args.outdir, exist_ok=True)
    setup_basic_logging()

    #### Log startup information
    log_parameter_summary(args, module)

    #### Plot boundary-flanking methylation levels
    logging.info(f"\n> Plotting telomere DNA methylation levels")
    run_methyl_plot(bedmethyl=args.bed_methyl,
                    boundary=args.boundary,
                    outdir=args.outdir,
                    out_prefix=args.out_prefix,
                    args=args)

    logging.info(f"\n**** Plotting completed! ****")


if __name__ == "__main__":

    methyl_plot_cli()


