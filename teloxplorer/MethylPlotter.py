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
        return {}
    
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


def calculate_smoothed_data(raw_data, window_size=1000, step_size=100, smooth_window=20):

    smoothed_data_list = []

    for (chr_arm, haplotype), group_data in raw_data.groupby(['chr_arm', 'haplotype']):

        clean_data = group_data.dropna(subset=['relative_position', 'percent_modified'])
        
        if len(clean_data) == 0:
            continue

        positions = clean_data['relative_position'].values
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
                    'haplotype': haplotype,
                    'relative_position': pos,
                    'percent_modified': val,
                    'data_type': 'smoothed'
                })
    
    return pd.DataFrame(smoothed_data_list)


def plot_ggline_wo_dots(data2plot, smoothed_data, region_size, hap_colors, hide_raw=False):

    unique_haps = data2plot['haplotype'].nunique()

    p = (
        ggplot(data2plot) +
        aes(x="relative_position", y="percent_modified")
    )

    if not hide_raw:
        p += geom_line(aes(color="haplotype"), size=0.2, alpha=0.3)

    p += geom_line(
        aes(x="relative_position", y="percent_modified", color="haplotype"),
        data=smoothed_data,
        size=0.8
    )
    
    p += scale_color_manual(values=hap_colors)

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
        x=f"Coordinate relative to telo-subtelo boundary (Unit: kb)",
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

    if unique_haps <= 1:
        p += theme(legend_position="none")

    return p

def run_methyl_plot(bedmethyl,boundary,outdir,out_prefix,args):

    # --- 1. Processing bedmethyl file ---
    methybed_header = [
        "chrom", "start_position", "end_position", "modified_base_code_and_motif",
        "score", "strand", "start_position1", "end_position1", "color",
        "N_valid_cov", "percent_modified", "N_mod", "N_canonical",
        "N_other_mod", "N_delete", "N_fail", "N_diff", "N_nocall"
    ]
    bedmethyl_df = pd.read_csv(bedmethyl, sep="\t", header=None, names=methybed_header)

    bedmethyl_df_tidy = bedmethyl_df.copy()

    # Filter out mitochondrial chromosomes
    mito_regex = r'(?i)^(chr|chromosome)?(M|MT|Mito)$'
    is_mito = bedmethyl_df_tidy['chrom'].str.match(mito_regex).fillna(False)
    bedmethyl_df_tidy = bedmethyl_df_tidy[~is_mito]

    # Filter out low cov site
    cols_to_keep = [
        "chrom", 
        "start_position", 
        "percent_modified"
    ]
    row_filter = (bedmethyl_df_tidy["score"] >= args.min_valid_cov)
    bedmethyl_df_tidy = bedmethyl_df_tidy.loc[row_filter, cols_to_keep]

    # --- 2. Processing boundary file ---
    boundary_df = pd.read_csv(boundary, sep="\t")
    boundary_df = boundary_df[["chrom", "tel_arm", "tel_length", "tel_boundary"]].drop_duplicates()

    boundary_df_wider = boundary_df.pivot(index="chrom", columns="tel_arm", values="tel_boundary").reset_index()
    boundary_df_wider.columns.name = None

    # --- 3. Join bedmethyl and boundary ---
    bedmethyl_join = pd.merge(bedmethyl_df_tidy, boundary_df_wider, on="chrom", how="left")

    # Add 'chr' prefix if not present
    bedmethyl_join["chrom"] = bedmethyl_join["chrom"].apply(
        lambda x: f"chr{x}" if not x.startswith("chr") else x
    )

    # Process phasing chromosomes
    chrom_split = bedmethyl_join["chrom"].str.partition("_")
    bedmethyl_join["chrom"] = chrom_split[0]
    bedmethyl_join["haplotype"] = chrom_split[2].replace("", "hap0")

    # Assign chromosome arm for each methyl site
    conditions = [
        bedmethyl_join["start_position"] <= bedmethyl_join["L"] + args.plot_flank,
        bedmethyl_join["start_position"] >= bedmethyl_join["R"] - args.plot_flank
    ]
    choices = ["L", "R"]
    bedmethyl_join["arm"] = np.select(conditions, choices, default=None)

    bedmethyl_join = bedmethyl_join.drop(columns=['L', 'R'])

    bedmethyl_join = bedmethyl_join.dropna(subset=["arm"])
    bedmethyl_join["arm_boundary"] = np.where(bedmethyl_join["arm"] == "L", bedmethyl_join["L"], bedmethyl_join["R"])
    bedmethyl_join["chr_arm"] = bedmethyl_join["chrom"] + "-" + bedmethyl_join["arm"]

    # Convert methyl position to relative coordinates
    bedmethyl_join["relative_position"] = bedmethyl_join["start_position"] - bedmethyl_join["arm_boundary"]

    # Set chromosome order based on user input or natural sort
    if args.plot_chr:
        logging.info(f"Using user-defined chromosome list and order: {args.plot_chr}")
        chrom_level = args.plot_chr.split(',')

        # Filter the DataFrame to only these chromosomes
        bedmethyl_join = bedmethyl_join[bedmethyl_join['chrom'].isin(chrom_level)].copy()

    else:
        chrom_unique = bedmethyl_join['chrom'].unique()
        chrom_level = natsorted(chrom_unique)

    chr_arm_level = [f"{c}-{arm}" for c in chrom_level for arm in ["L", "R"]]

    # Set chromosomal arm tags
    if args.arm_label == "p,q":
        bedmethyl_join["arm"] = bedmethyl_join["arm"].replace({"L": "p", "R": "q"})
        bedmethyl_join["chr_arm"] = bedmethyl_join["chrom"] + bedmethyl_join["arm"]
        chr_arm_level = [c.replace("-L", "p").replace("-R", "q") for c in chr_arm_level]

    # Set chromosome factor level for raw data
    bedmethyl_join["chr_arm"] = pd.Categorical(
        bedmethyl_join["chr_arm"],
        categories=chr_arm_level,
        ordered=True
    )

    # --- 4. methyl plot ---
    # Complete data for plotting
    all_combs = pd.MultiIndex.from_product([bedmethyl_join['haplotype'].unique(), bedmethyl_join['chr_arm'].cat.categories], names=['haplotype', 'chr_arm']).to_frame(index=False)
    bedmethyl_2plot = pd.merge(all_combs, bedmethyl_join, on=['haplotype', 'chr_arm'], how='left')

    # Calculate smoothed data using methylartist method
    smoothed_data = calculate_smoothed_data(
        bedmethyl_2plot, 
        window_size=args.window_size,
        step_size=args.step_size, 
        smooth_window=args.smooth_window
    )

    # Set chromosomes factor level for smoothed data
    bedmethyl_2plot['chr_arm'] = pd.Categorical(
        bedmethyl_2plot['chr_arm'],
        categories=chr_arm_level,
        ordered=True
    )
    smoothed_data['chr_arm'] = pd.Categorical(
        smoothed_data['chr_arm'],
        categories=chr_arm_level,
        ordered=True
    )

    # Writing plot data
    raw_data_output= os.path.join(outdir, f"{out_prefix}.raw_data_for_plot.tsv")
    smoothed_data_output = os.path.join(outdir, f"{out_prefix}.smoothed_data_for_plot.tsv")

    logging.info(f"> Writing raw plotting data to: {raw_data_output}")
    bedmethyl_2plot.to_csv(raw_data_output, sep='\t', index=False)

    logging.info(f"> Writing smoothed plotting data to: {smoothed_data_output}")
    smoothed_data.to_csv(smoothed_data_output, sep='\t', index=False)

    # plotting
    matplotlib.set_loglevel("error")
    matplotlib.rcParams['pdf.fonttype'] = 42
    matplotlib.rcParams['font.family'] = 'Arial'
    matplotlib.rcParams['font.size'] = 10
    logging.getLogger('fontTools').setLevel(logging.WARNING)
    logging.getLogger('matplotlib').setLevel(logging.WARNING)

    if args.hap_colors:
        logging.info(f"Using custom haplotype colors: {args.hap_colors}")
        color_values = args.hap_colors.split(',')
    else:
        color_values = ['#E41A1C', '#377EB8']

    p = plot_ggline_wo_dots(
        bedmethyl_2plot,
        smoothed_data,
        args.plot_flank,
        color_values,
        hide_raw=args.hide_raw
    )

    p += geom_vline(xintercept=0, linetype="dashed",)

    output = os.path.join(outdir, f"{out_prefix}.tel_methyl.pdf")
    p.save(output, width=args.width, height=args.height, dpi=300, verbose=False)

    logging.info(f"Plot completed and saved to: {output}")

def setup_common_args(parser):

    group = parser.add_argument_group("Input/Output Options")
    group.add_argument("--bedmethyl",
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
                       help="Minimum valid coverage (Nmod + Nother_mod + Ncanonical), "
                            "also used as the score in the bedMethyl [default: %(default)s]")
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
    group.add_argument("--plot-flank",
                       type=int,
                       default=20000,
                       help="Flanking size (bp) around telomere-subtelomere boundary for plotting [default: %(default)s]")
    group.add_argument("--arm-label",
                       choices=["p,q", "L,R"],
                       default="L,R",
                       help="Chromosome arm labels [default: %(default)s]")
    group.add_argument("--plot-chr",
                       default=None,
                       help="A comma-separated list of chromosomes to plot in the given order"
                            "(e.g. 'chr1,chrX') [default: %(default)s]")
    group.add_argument("--hap-colors",
                       default=None,
                       help="List of colors for haplotypes (comma-separated). "
                            "If only one haplotype exists, the first color provided is used. "
                            "Defaults to an automatic color palette if not set. [default: %(default)s]")
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
    run_methyl_plot(bedmethyl=args.bedmethyl,
                    boundary=args.boundary,
                    outdir=args.outdir,
                    out_prefix=args.out_prefix,
                    args=args)

    logging.info(f"\n**** Plotting completed! ****")


if __name__ == "__main__":

    methyl_plot_cli()


