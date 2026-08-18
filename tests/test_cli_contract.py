from __future__ import annotations

from dataclasses import fields
from pathlib import Path
from tempfile import TemporaryDirectory
from unittest.mock import MagicMock, call, patch

import pandas as pd
from typer.main import get_command
from typer.testing import CliRunner

from telox import telox_align, telox_asm, telox_length, viz_length, viz_methyl, viz_reads, viz_tvr_hap
from telox.cli import app


RUNNER = CliRunner()


def _invoke_plot(command: str, arguments: list[str], target: str):
    logger = MagicMock()
    with TemporaryDirectory() as directory:
        with (
            patch("telox.cli.log_utils.init_logger", return_value=logger),
            patch(target) as plotter,
        ):
            result = RUNNER.invoke(
                app,
                [command, *arguments, "--outdir", directory],
            )
    assert result.exit_code == 0, result.output or repr(result.exception)
    assert plotter.call_count == 1
    return plotter.call_args.kwargs


def _option_names(command: str) -> set[str]:
    click_command = get_command(app).commands[command]
    return {
        option
        for parameter in click_command.params
        for option in (*parameter.opts, *parameter.secondary_opts)
    }


def _expect_exactly_one_error(command: str, arguments: list[str]) -> None:
    result = RUNNER.invoke(app, [command, *arguments])
    assert result.exit_code == 1
    assert "Provide exactly one of" in result.output


def _expect_value_error(callback, message: str) -> None:
    try:
        callback()
    except ValueError as exc:
        assert message in str(exc)
    else:
        raise AssertionError("Expected ValueError")


def test_plot_reads_defaults_come_from_config() -> None:
    kwargs = _invoke_plot(
        "plot-reads",
        ["--tvr", "reads.tsv"],
        "telox.viz_reads.plot_tvr_reads",
    )
    assert kwargs["config"] == viz_reads.TVRReadPlotter()


def test_plot_length_rasterize_cli_is_tristate() -> None:
    assert {"--rasterize", "--no-rasterize"} <= _option_names("plot-length")

    automatic = _invoke_plot(
        "plot-length",
        ["--input", "length.tsv"],
        "telox.viz_length.plot_length",
    )["config"]
    forced_raster = _invoke_plot(
        "plot-length",
        ["--input", "length.tsv", "--rasterize"],
        "telox.viz_length.plot_length",
    )["config"]
    forced_vector = _invoke_plot(
        "plot-length",
        ["--input", "length.tsv", "--no-rasterize"],
        "telox.viz_length.plot_length",
    )["config"]

    assert automatic == viz_length.LengthPlotter()
    assert automatic.rasterize is None
    assert forced_raster.rasterize is True
    assert forced_vector.rasterize is False


def test_plot_reads_forwards_every_cli_override() -> None:
    kwargs = _invoke_plot(
        "plot-reads",
        [
            "--tvr", "reads.tsv",
            "--consensus", "consensus.tsv",
            "--mods", "mods.tsv",
            "--motif", "TG{1,3}",
            "--top-tvr", "17",
            "--chroms", "chr2,chr3",
            "--pool-chroms",
            "--tvr-colors", "motifs.tsv",
            "--xlim=-123,456",
            "--read-thickness", "0.61",
            "--consensus-thickness", "2",
            "--legend-ncol", "7",
            "--show-outliers",
            "--show-unmethylated",
            "--sample-label", "READ_SENTINEL",
            "--width", "6.1",
            "--height", "4.2",
            "--format", "png",
            "--prefix", "reads-prefix",
        ],
        "telox.viz_reads.plot_tvr_reads",
    )
    config = kwargs["config"]
    assert config.canonical_motif == "TG{1,3}"
    assert config.top_tvr == 17
    assert config.chroms == "chr2,chr3"
    assert config.pool_chroms is True
    assert Path(config.tvr_colors) == Path("motifs.tsv")
    assert config.xlim == "-123,456"
    assert config.read_thickness == 0.61
    assert config.consensus_thickness == 2.0
    assert config.legend_ncol == 7
    assert config.show_outliers is True
    assert config.show_unmethylated_c is True
    assert config.sample_label == "READ_SENTINEL"
    assert config.width == 6.1
    assert config.height == 4.2
    assert config.output_format == "png"
    assert kwargs["tvr_file"] == Path("reads.tsv")
    assert kwargs["consensus_file"] == Path("consensus.tsv")
    assert kwargs["mods_file"] == Path("mods.tsv")
    assert kwargs["sample_sheet"] is None
    assert kwargs["prefix"] == "reads-prefix"


def test_plot_hap_defaults_come_from_config() -> None:
    kwargs = _invoke_plot(
        "plot-tvr-hap",
        ["--consensus", "consensus.tsv"],
        "telox.viz_tvr_hap.plot_tvr_haps",
    )
    assert kwargs["config"] == viz_tvr_hap.TVRHapPlotter()
    assert kwargs["config"].rasterize is None
    assert kwargs["write_similarity_score"] is False


def test_plot_hap_rasterize_cli_is_tristate() -> None:
    forced_raster = _invoke_plot(
        "plot-tvr-hap",
        ["--consensus", "consensus.tsv", "--rasterize"],
        "telox.viz_tvr_hap.plot_tvr_haps",
    )["config"]
    forced_vector = _invoke_plot(
        "plot-tvr-hap",
        ["--consensus", "consensus.tsv", "--no-rasterize"],
        "telox.viz_tvr_hap.plot_tvr_haps",
    )["config"]

    assert forced_raster.rasterize is True
    assert forced_vector.rasterize is False


def test_plot_xlim_none_selects_each_panel_full_data_span() -> None:
    reads = _invoke_plot(
        "plot-reads",
        ["--tvr", "reads.tsv", "--xlim", "none"],
        "telox.viz_reads.plot_tvr_reads",
    )
    haps = _invoke_plot(
        "plot-tvr-hap",
        ["--consensus", "consensus.tsv", "--xlim", "none"],
        "telox.viz_tvr_hap.plot_tvr_haps",
    )

    assert reads["config"].xlim == "none"
    assert haps["config"].xlim is None


def test_row_label_columns_parse_csv() -> None:
    assert viz_tvr_hap._row_label_columns(
        "chrom_phase, chrom, arm", multi_sample=False
    ) == ("chrom_phase", "chrom", "arm")


def test_annotation_columns_accept_none() -> None:
    assert viz_tvr_hap._annotation_columns("chrom, chrom_phase") == (
        "chrom", "chrom_phase"
    )
    assert viz_tvr_hap._annotation_columns("none") == ()

    kwargs = _invoke_plot(
        "plot-tvr-hap",
        ["--consensus", "consensus.tsv", "--annotation-columns", "none"],
        "telox.viz_tvr_hap.plot_tvr_haps",
    )
    assert kwargs["config"].annotation_columns == "none"


def test_plot_hap_forwards_every_cli_override() -> None:
    kwargs = _invoke_plot(
        "plot-tvr-hap",
        [
            "--consensus", "consensus.tsv",
            "--motif", "TTTAGGG",
            "--top-tvr", "13",
            "--chroms", "chr4,chr5",
            "--pool-arms",
            "--pool-chroms",
            "--tvr-colors", "motifs.tsv",
            "--hap-filter", "6,0.31,0.19",
            "--xlim", "1555",
            "--track-thickness", "0.67",
            "--no-cluster-rows",
            "--cluster-proximal-bp", "333",
            "--no-heatmap",
            "--heatmap-cmap", "viridis",
            "--row-label-columns", "chrom,arm",
            "--annotation-columns", "chrom,arm",
            "--annotation-file", "annotations.tsv",
            "--annotation-colors", "annotation-colors.tsv",
            "--legend-ncol", "8",
            "--write-similarity-score",
            "--rasterize",
            "--width", "6.2",
            "--height", "5.1",
            "--format", "svg",
            "--sample-label", "HAP_SENTINEL",
            "--prefix", "hap-prefix",
        ],
        "telox.viz_tvr_hap.plot_tvr_haps",
    )
    config = kwargs["config"]
    assert config.canonical_motif == "TTTAGGG"
    assert config.top_tvr == 13
    assert config.chroms == "chr4,chr5"
    assert config.pool_arms is True
    assert config.pool_chroms is True
    assert Path(config.tvr_colors) == Path("motifs.tsv")
    assert config.hap_filter == viz_tvr_hap.HaplotypeFilter(6, 0.31, 0.19)
    assert config.xlim == 1555.0
    assert config.track_thickness == 0.67
    assert config.cluster_rows is False
    assert config.cluster_proximal_bp == 333
    assert config.no_heatmap is True
    assert config.heatmap_cmap == "viridis"
    assert config.row_label_columns == "chrom,arm"
    assert config.annotation_columns == "chrom,arm"
    assert Path(config.annotation_file) == Path("annotations.tsv")
    assert Path(config.annotation_colors) == Path("annotation-colors.tsv")
    assert config.legend_ncol == 8
    assert config.rasterize is True
    assert config.width == 6.2
    assert config.height == 5.1
    assert config.output_format == "svg"
    assert config.sample_label == "HAP_SENTINEL"
    assert kwargs["consensus_file"] == Path("consensus.tsv")
    assert kwargs["sample_sheet"] is None
    assert kwargs["write_similarity_score"] is True
    assert kwargs["prefix"] == "hap-prefix"


def test_plot_methyl_defaults_come_from_config() -> None:
    kwargs = _invoke_plot(
        "plot-methyl",
        ["--mods", "mods.tsv"],
        "telox.viz_methyl.plot_arm_methyl",
    )
    assert kwargs["config"] == viz_methyl.MethylPlotter()


def test_plot_methyl_forwards_every_cli_override() -> None:
    kwargs = _invoke_plot(
        "plot-methyl",
        [
            "--mods", "mods.tsv",
            "--xlim=-200,400",
            "--bin-size", "20",
            "--smooth-window", "80",
            "--show-sd",
            "--no-heatmap",
            "--heatmap-cmap", "viridis",
            "--min-bin-reads", "7",
            "--min-bin-arms", "9",
            "--width", "6.3",
            "--height", "4.4",
            "--format", "png",
            "--prefix", "methyl-prefix",
        ],
        "telox.viz_methyl.plot_arm_methyl",
    )
    config = kwargs["config"]
    assert config.xlim == "-200,400"
    assert config.bin_size == 20
    assert config.smooth_window == 80
    assert config.show_sd is True
    assert config.show_heatmap is False
    assert config.heatmap_cmap == "viridis"
    assert config.min_bin_reads == 7
    assert config.min_bin_arms == 9
    assert config.width == 6.3
    assert config.height == 4.4
    assert config.output_format == "png"
    assert kwargs["mods_file"] == Path("mods.tsv")
    assert kwargs["prefix"] == "methyl-prefix"


def test_supported_sample_sheet_inputs_are_forwarded() -> None:
    reads = _invoke_plot(
        "plot-reads",
        ["--sample-sheet", "tvr-samples.tsv"],
        "telox.viz_reads.plot_tvr_reads",
    )
    assert reads["tvr_file"] is None
    assert reads["sample_sheet"] == Path("tvr-samples.tsv")

    hap = _invoke_plot(
        "plot-tvr-hap",
        ["--sample-sheet", "hap-samples.tsv"],
        "telox.viz_tvr_hap.plot_tvr_haps",
    )
    assert hap["consensus_file"] is None
    assert hap["sample_sheet"] == Path("hap-samples.tsv")

def test_paired_boolean_overrides_reach_configs() -> None:
    hap = _invoke_plot(
        "plot-tvr-hap",
        [
            "--consensus", "consensus.tsv",
            "--cluster-rows",
            "--no-rasterize",
        ],
        "telox.viz_tvr_hap.plot_tvr_haps",
    )["config"]
    assert hap.cluster_rows is True
    assert hap.rasterize is False

    methyl = _invoke_plot(
        "plot-methyl",
        ["--mods", "mods.tsv", "--no-sd", "--show-heatmap"],
        "telox.viz_methyl.plot_arm_methyl",
    )["config"]
    assert methyl.show_sd is False
    assert methyl.show_heatmap is True


def test_haplotype_filter_cli_preset_and_validation() -> None:
    config = _invoke_plot(
        "plot-tvr-hap",
        ["--consensus", "consensus.tsv", "--hap-filter", "hprc2"],
        "telox.viz_tvr_hap.plot_tvr_haps",
    )["config"]
    assert config.hap_filter == viz_tvr_hap.HaplotypeFilter(2, 0.4, 0.2)

    invalid = RUNNER.invoke(
        app,
        ["plot-tvr-hap", "--consensus", "consensus.tsv", "--hap-filter", "2,0.4"],
    )
    assert invalid.exit_code == 2
    assert "MAX,SINGLE,MULTI" in invalid.output


def test_haplotype_filter_logs_policy_and_row_count() -> None:
    logger = MagicMock()
    policy = viz_tvr_hap.HaplotypeFilter(2, 0.4, 0.2)
    haplotypes = pd.DataFrame({
        "sample_id": ["sample"] * 3,
        "chrom": ["chr1"] * 3,
        "chrom_phase": ["mat"] * 3,
        "arm": ["p"] * 3,
        "tvr_hap": ["C1", "C2", "C3"],
        "read_support": [10, 5, 1],
        "hap_frequency": [0.7, 0.25, 0.05],
        "hap_length": [100, 100, 100],
        "motif_blocks": ["TTAGGG:1"] * 3,
        "tokens": [(('TTAGGG', 1),)] * 3,
    })
    viz_tvr_hap.prepare_haplotypes(
        haplotypes,
        viz_tvr_hap.TVRHapPlotter(hap_filter=policy, no_heatmap=True),
        logger=logger,
    )
    logger.info.assert_has_calls([
        call("Haplotype filter (max=2, single>=0.4, multi>=0.2): 3 -> 2 haplotypes")
    ])


def test_removed_options_and_fields_are_absent() -> None:
    assert "--rasterize" not in _option_names("plot-reads")
    assert {"--chroms", "--smooth", "--rasterize", "--sample-sheet"}.isdisjoint(
        _option_names("plot-methyl")
    )
    assert {"--mm2-opts", "--min-mapq"}.isdisjoint(_option_names("reads"))
    assert {"--min-cov-per-arm", "--min-valid-arms"}.isdisjoint(
        _option_names("plot-methyl")
    )
    assert {"--min-bin-reads", "--min-bin-arms"}.issubset(
        _option_names("plot-methyl")
    )
    assert {
        "--filter-noise", "--max-haps-per-arm",
        "--min-freq-single", "--min-freq-multi",
    }.isdisjoint(_option_names("plot-tvr-hap"))
    assert "--hap-filter" in _option_names("plot-tvr-hap")
    assert {"--strip-chr-prefix", "--keep-chr-prefix"}.isdisjoint(
        _option_names("plot-tvr-hap")
    )
    assert {"--validate-lengths", "--no-validate-lengths"}.isdisjoint(
        _option_names("plot-tvr-hap")
    )
    assert "validate_lengths" not in {field.name for field in fields(viz_tvr_hap.TVRHapPlotter)}
    assert {field.name for field in fields(telox_length.LengthConfig)}.isdisjoint(
        {"width", "height"}
    )


def test_plot_reads_outlier_filter_is_configurable() -> None:
    with TemporaryDirectory() as directory:
        path = Path(directory) / "reads.tsv"
        pd.DataFrame({
            "read_id": ["clustered", "outlier"],
            "chrom": ["chr1", "chr1"],
            "chrom_phase": ["mat", "mat"],
            "arm": ["p", "p"],
            "tvr_hap": ["C1", "Outlier"],
            "tel_length": [100, 90],
            "tvrs": [".", "."],
        }).to_csv(path, sep="\t", index=False)

        hidden = viz_reads._load_tvr(path, "sample")
        shown = viz_reads._load_tvr(path, "sample", show_outliers=True)

    assert hidden["read_id"].tolist() == ["clustered"]
    assert shown["read_id"].tolist() == ["clustered", "outlier"]


def test_exactly_one_validation_for_input_sources() -> None:
    _expect_exactly_one_error("run", ["--ref", "ref.fa"])
    _expect_exactly_one_error(
        "run", ["--ref", "ref.fa", "--fastq", "reads.fq", "--bam", "reads.bam"]
    )
    _expect_exactly_one_error("align", ["--ref", "ref.fa"])
    _expect_exactly_one_error(
        "align", ["--ref", "ref.fa", "--fastq", "reads.fq", "--modbam", "mods.bam"]
    )
    _expect_exactly_one_error("reads", [])
    _expect_exactly_one_error("reads", ["--fastq", "reads.fq", "--bam", "reads.bam"])
    _expect_exactly_one_error("plot-reads", [])
    _expect_exactly_one_error(
        "plot-reads", ["--tvr", "reads.tsv", "--sample-sheet", "samples.tsv"]
    )
    _expect_exactly_one_error("plot-tvr-hap", [])
    _expect_exactly_one_error(
        "plot-tvr-hap", ["--consensus", "consensus.tsv", "--sample-sheet", "samples.tsv"]
    )


def test_plot_methyl_requires_mods() -> None:
    result = RUNNER.invoke(app, ["plot-methyl"])
    assert result.exit_code == 2
    assert "--mods is required" in result.output


def test_plot_reads_rejects_auxiliary_tracks_in_sample_sheet_mode() -> None:
    with TemporaryDirectory() as directory:
        result = RUNNER.invoke(
            app,
            [
                "plot-reads", "--sample-sheet", "samples.tsv",
                "--mods", "mods.tsv", "--outdir", directory,
            ],
        )
    assert result.exit_code == 1
    assert "sample-sheet mode supports TVR files only" in result.output


def test_align_accepts_modbam_as_the_only_input() -> None:
    logger = MagicMock()
    with TemporaryDirectory() as directory:
        with (
            patch("telox.cli.log_utils.init_logger", return_value=logger),
            patch("telox.telox_align.run_telogrep") as telogrep,
            patch("telox.telox_align.align_reads") as align_reads,
        ):
            result = RUNNER.invoke(
                app,
                [
                    "align", "--modbam", "mods.bam", "--ref", "ref.fa",
                    "--outdir", directory,
                ],
            )
    assert result.exit_code == 0, result.output or repr(result.exception)
    assert telogrep.call_args.kwargs["modbam_input"] == Path("mods.bam")
    assert telogrep.call_args.kwargs["fastq_input"] is None
    assert telogrep.call_args.kwargs["bam_input"] is None
    assert align_reads.call_count == 1


def test_presets_are_strict_in_cli_and_module_apis() -> None:
    invalid = RUNNER.invoke(
        app,
        ["align", "--preset", "unknown", "--fastq", "reads.fq", "--ref", "ref.fa"],
    )
    assert invalid.exit_code == 2
    assert "Invalid value" in invalid.output

    unsupported_asm = RUNNER.invoke(
        app,
        ["asm", "--preset", "human-r9", "--assembly", "assembly.fa"],
    )
    assert unsupported_asm.exit_code == 2
    assert "Invalid value" in unsupported_asm.output

    _expect_value_error(
        lambda: telox_align.AlignConfig.from_params(preset_name="unknown"),
        "Unknown preset",
    )
    _expect_value_error(
        lambda: telox_asm.AsmConfig.from_params(preset_name="human-r9"),
        "Unknown preset",
    )


if __name__ == "__main__":
    for name, test in sorted(globals().items()):
        if name.startswith("test_") and callable(test):
            test()
            print(f"PASS {name}")
