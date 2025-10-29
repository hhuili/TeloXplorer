import argparse
import re
import os
import logging
import pandas as pd
import matplotlib
matplotlib.use('Agg')
from natsort import natsorted
from plotnine import *
from plotnine.themes import theme
from mizani.breaks import breaks_extended
from .utilities import setup_basic_logging,log_parameter_summary

def run_length_plot(length_input,outdir,out_prefix,args):

    # setup logging
    matplotlib.set_loglevel("error")
    matplotlib.rcParams['pdf.fonttype'] = 42
    matplotlib.rcParams['font.family'] = 'Arial'
    matplotlib.rcParams['font.size'] = 10
    logging.getLogger('fontTools').setLevel(logging.WARNING)
    logging.getLogger('matplotlib').setLevel(logging.WARNING)

    # output file
    output = os.path.join(outdir, f"{out_prefix}.tel_length.pdf")

    # --- 1. Data Loading and Processing ---
    tel_length = pd.read_csv(length_input, sep='\t')
    tel_length_tidy = tel_length.copy()

    # Remove 'chr' prefix from chromosome names
    tel_length_tidy['chrom'] = tel_length_tidy['chrom'].str.replace(
        r'^(chr|Chr|CHR|chromosome|Chromosome|CHROMOSOME)\s*',
        '',
        regex=True,
        flags=re.IGNORECASE
    )

    # Identify haplotype chromosomes and split chromosome names safely
    tel_length_tidy['is_hap'] = tel_length_tidy['chrom'].apply(
        lambda x: 'yes' if '_' in x else 'no'
    )

    # Safely split chromosome names to handle cases without '_'
    split_chrom = tel_length_tidy['chrom'].str.split('_', expand=True)

    # Handle chromosome splitting safely
    tel_length_tidy['chr'] = split_chrom[0]

    # Create haplotype column - use second part if exists, otherwise empty string
    if split_chrom.shape[1] > 1:
        tel_length_tidy['haplotype'] = split_chrom[1].fillna('')
    else:
        tel_length_tidy['haplotype'] = ''

    # Filter out mitochondrial chromosomes
    tel_length_tidy = tel_length_tidy[~tel_length_tidy['chr'].isin(["M", "Mito", "Mt"])]

    # Set chromosome order based on user input or natural sort
    if args.plot_chr:
        logging.info(f"Using user-defined chromosome list and order: {args.plot_chr}")
        chrom_levels = args.plot_chr.split(',')

        # Filter the DataFrame to only these chromosomes
        tel_length_tidy = tel_length_tidy[tel_length_tidy['chr'].isin(chrom_levels)].copy()

    else:
        unique_chroms = tel_length_tidy['chr'].unique()

        chrom_levels = natsorted(unique_chroms)

    # Set categorical order for plotting
    tel_length_tidy['chr'] = pd.Categorical(
        tel_length_tidy['chr'],
        categories=chrom_levels,
        ordered=True
    )

    # Drop any rows that became NaN because their 'chr' was not in chrom_levels
    tel_length_tidy.dropna(subset=['chr'], inplace=True)

    tel_length_tidy['tel_arm'] = tel_length_tidy['tel_arm'].astype('category')

    dodge = position_dodge2(width=0.8,preserve='single')
    jitterdodge = position_jitterdodge(dodge_width=0.7, jitter_width=0.2)

    if args.arm_colors:
        logging.info(f"Using custom arm colors: {args.arm_colors}")
        color_values = args.arm_colors.split(',')
    else:
        color_values = ['#E41A1C', '#377EB8'] # Default colors

    common_layers = [
        geom_jitter(
            position=jitterdodge,
            shape='.',
            stroke=0,
            size=2,
            alpha=0.4
        ),
        geom_boxplot(
            position=dodge,
            outlier_shape='', # No outliers
            size=0.3, # Line width for boxplot borders and whiskers
            fatten=1, # Controls median line thickness
            fill='none', # No fill
            width=0.7, # Box width
        ),
        scale_color_manual(values=color_values),
        scale_y_continuous(
            breaks=breaks_extended(Q=[1, 2, 5, 10]),
            labels=lambda x: [f"{int(val/1000)}" for val in x]
        ),
        scale_x_discrete(),
        labs(x="Chromosome", y="Telomere length (kb)",color="Arm"),
        theme_matplotlib() +
        theme(
            panel_border=element_rect(color='black'),
            legend_position="right"
        )
    ]

    # --- 2b. Define vertical divider lines between chromosomes ---
    # Get the number of chromosome levels
    num_chroms = len(chrom_levels)

    # Calculate intercepts (at positions 1.5, 2.5, ..., N-1.5)
    vline_intercepts = [i + 1.5 for i in range(num_chroms-1)]

    # Create the vline layer
    chromosome_dividers = geom_vline(
        xintercept=vline_intercepts,
        linetype="dashed",
        color="black",
        alpha=0.2,
        size=0.3
    )

    # --- 3. Create Appropriate Plot Based on Data Type ---
    has_haplotype_data = 'yes' in tel_length_tidy['is_hap'].unique()

    # Get the first chromosome name from the *ordered* levels
    x_pos = chrom_levels[0] 
    y_min = tel_length_tidy['tel_length'].min()
    y_max = tel_length_tidy['tel_length'].max()
    y_pos = y_max - (y_max - y_min) * 0.05

    if has_haplotype_data:

        hap_medians_df = tel_length_tidy.groupby('haplotype')['tel_length'].median().reset_index()
        hap_medians_df['median_label'] = hap_medians_df['tel_length'].apply(lambda x: f"Genome-wide median: {x:.0f} bp")

        median_text_layer = geom_text(
            data=hap_medians_df,
            mapping=aes(label='median_label'),
            x=x_pos,
            y=y_pos,
            ha="left",
            va="top",
            color="black",
            size=8,
            inherit_aes=False
        )

        p = (
            ggplot(
                tel_length_tidy,
                aes(x='chr', y='tel_length', color='tel_arm')
            ) +
            chromosome_dividers +
            median_text_layer +
            common_layers +
            facet_grid('haplotype ~ .')
        )

        p.save(filename=output, width=args.width, height=args.height, dpi=300)

    else:

        global_median = tel_length_tidy['tel_length'].median()
        median_label = f"Genome-wide median: {global_median:.0f} bp"

        median_text_annotation = annotate(
            "text",
            x=x_pos,
            y=y_pos,
            label=median_label,
            ha="left",
            va="top",
            color="black",
            size=8
        )

        p = (
            ggplot(
                tel_length_tidy,
                aes(x='chr', y='tel_length', color='tel_arm')
                ) +
            chromosome_dividers +
            median_text_annotation +
            common_layers
        )

        p.save(filename=output, width=args.width, height=args.height, dpi=300, verbose=False)

    logging.info(f"Plot completed and saved to: {output}")

def setup_common_args(parser):

    group = parser.add_argument_group("Input/Output Options")
    group.add_argument("-i", "--length-input",
                       required=True,
                       help="Telomere length file from telox-length")
    group.add_argument(f"--out-prefix",
                       default="plot",
                       help="Output prefix [default: %(default)s].")
    group.add_argument(f"--outdir",
                       default=".",
                       help="Output directory [default: %(default)s].")

    return parser

def setup_length_plotting_args(parser):

    group_name = 'Length Plotting Options'
    group = parser.add_argument_group(group_name, 'Telomere legnth plotting parameters')

    group.add_argument("--arm-colors",
                       default=None,
                       help="Chromosome arm colors (comma-separated) [default: %(default)s]")
    group.add_argument("--plot-chr",
                       default=None,
                       help="A comma-separated list of chromosomes to plot in the given order (no 'chr' prefix, e.g. '1,2,3,X') [default: %(default)s]")
    group.add_argument("-W", "--width",
                       type=int,
                       default=8,
                       help="Plot width (inches) [default: %(default)s]")
    group.add_argument("-H", "--height",
                       type=int,
                       default=3,
                       help="Plot height (inches) [default: %(default)s]")

    return parser

def parse_arguments():

    parser = argparse.ArgumentParser(
        description="Telomere length plotting parameters."
    )

    parser = setup_common_args(parser)
    parser = setup_length_plotting_args(parser)

    return parser.parse_args()

def length_plot_cli():

    module = "length-plot"
    args = parse_arguments()

    os.makedirs(args.outdir, exist_ok=True)
    setup_basic_logging()

    #### Log startup information
    log_parameter_summary(args, module)

    #### Plot telomere legnth across chromosomes
    logging.info(f"\n> Plotting telomere legnth across chromosomes")
    run_length_plot(length_input=args.length_input,
                    outdir=args.outdir,
                    out_prefix=args.out_prefix,
                    args=args)

    logging.info(f"\n**** Plotting completed! ****")

if __name__ == '__main__':

    length_plot_cli()

