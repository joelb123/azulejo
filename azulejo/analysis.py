# -*- coding: utf-8 -*-
"""Data analysis and plotting."""
# standard library imports
import sys
from pathlib import Path

# third-party imports
import click
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import seaborn as sns

# first-party imports
from loguru import logger

# module imports
from . import cli
from . import click_loguru
from .common import cluster_set_name
from .common import homo_degree_dist_filename

# Global constants
NAME = "azulejo"
STATFILE_SUFFIX = f"-{NAME}_stats.tsv"
ANYFILE_SUFFIX = f"-{NAME}_ids-any.tsv"
ALLFILE_SUFFIX = f"-{NAME}_ids-all.tsv"
IDENT_LOG_MIN = -3
IDENT_LOG_MAX = 0
EPSILON = 0.000001
FILETYPE = "png"
MAX_BINS = 10


def make_histogram(dist, name, log10=False):
    """Do histogram plot with kernel density estimate."""
    dist = dist[dist < 10]
    mean = dist.mean()
    if log10:
        dist = np.log10(dist)
    # if len(dist) < MAX_BINS:
    #    bins = len(dist)
    # else:
    #    bins = MAX_BINS
    sns.distplot(
        dist,
        bins=None,
        rug=False,
        kde=False,
        norm_hist=True,
        rug_kws={"color": "b"},
        kde_kws={"color": "k", "linewidth": 1, "label": "KDE"},
        hist_kws={
            "histtype": "step",
            "linewidth": 2,
            "alpha": 1,
            "color": "b",
            #'range':(0,20)
        },
    )
    plt.title(f"{name} histogram of {len(dist):d} values, mean={mean:.1f}")
    if log10:
        plt.xlabel("log " + name)
    else:
        plt.xlabel(name)
    plt.ylabel("Frequency")
    # for ext in PLOT_TYPES:
    #    plt.savefig(LOG_PATH / ('%s-histogram.' % (name.rstrip('%')) + ext),
    #                bbox_inches='tight')
    plt.show()
    # plt.close('all')


def tick_function(tick_x):
    """Compute ticks."""
    tick_x = tick_x * 3.0 - 3
    vals = [
        (f"{v:f}").rstrip("0").rstrip(".")
        for v in (1.0 - 10 ** tick_x) * 100.0
    ]
    ticks = [f"{v}%" for v in vals]
    return ticks


def log_deriv(xvals, yvals):
    """Compute the logarithmic derivative."""
    log_x = -1.0 * np.log10(xvals + EPSILON)
    log_y = np.log10(yvals)
    return np.gradient(log_y) / np.gradient(log_x)


@cli.command()
@click_loguru.init_logger(logfile=False)
@click.argument("instemlist")
def analyze_clusters(
    dirname, instemlist, label, reference=None, on_id=None, match_type=None
):
    """Statistics of clustering as function of identity."""
    if match_type is None:
        matches = ["all", "any"]
    else:
        matches = [match_type]
    uniques = {}
    divergence = {}
    dirpath = Path(dirname)
    div_dist = {"all": {"ref": 0.0}, "any": {"ref": 0.0}}
    print("ref=", reference)
    for stem in instemlist:
        print("stem=", stem)
        paths = {
            "all": dirpath / (stem + ALLFILE_SUFFIX),
            "any": dirpath / (stem + ANYFILE_SUFFIX),
            "stat": dirpath / (stem + STATFILE_SUFFIX),
        }
        stats = pd.read_csv(paths["stat"], sep="\t", index_col=0)
        uniques[stem] = stats["unique_seqs"].iloc[0]
        divergence[stem] = stats["divergence"]
        if on_id is None:
            div_dist["all"][stem] = log_deriv(
                divergence[stem], stats["clusters"]
            )
            div_dist["any"][stem] = None
            if stem == reference:
                div_dist["all"]["ref"] = div_dist["all"][stem]
                div_dist["any"]["ref"] = None
        else:
            for match in ["any", "all"]:
                data = pd.read_csv(paths[match], sep="\t", index_col=0)
                try:
                    div_dist[match][stem] = log_deriv(
                        divergence[stem], data.loc[on_id]
                    )
                except KeyError:  # this label wasn't found
                    div_dist[match][stem] = None
                if stem == reference:
                    div_dist[match]["ref"] = div_dist[match][stem]
    #
    # Make the plots
    #
    plt.style.use("seaborn-whitegrid")
    axis_dict = {}
    fig, axes = plt.subplots(len(matches))
    try:
        for axis, i in enumerate(axes):
            axis_dict[matches[i]] = axis
            loweraxis = axes[1]
    except TypeError:
        axis_dict[matches[0]] = axes
        loweraxis = axes
    for stem in instemlist:
        for match in matches:
            if div_dist[match][stem] is None:
                continue
            axis_dict[match].plot(
                divergence[stem],
                div_dist[match][stem] - div_dist[match]["ref"],
                label=f"{stem.replace(label + '.', '')}",
            )
            # uniques[stem]/1000.))
    if reference is None:
        if on_id is None:
            title = f"{label} Divergence Distribution"
            outfilestem = f"{label}_divergence_dist."
        else:
            title = f'{label} Divergence Distribution on "{on_id}"'
            outfilestem = f"{label}_divergence_dist_{on_id}."
    else:
        if on_id is None:
            title = (
                f"{label}_Differential Divergence Distribution vs. {reference}"
            )
            outfilestem = f"{label}_divergence_dist_vs{reference}."
        else:
            title = (
                f'{label} Differential Divergence Distribution on "{on_id}"'
                f" vs. {reference}"
            )
            outfilestem = f"{label}_divergence_dist_on_{on_id}_vs_ref."
    if reference is None:
        fig.text(
            0.02,
            0.5,
            "Logarithmic Derivative on Clusters",
            ha="center",
            va="center",
            rotation="vertical",
        )
    else:
        fig.text(
            0.02,
            0.5,
            "Logarithmic Derivative Difference on Clusters",
            ha="center",
            va="center",
            rotation="vertical",
        )
    if len(matches) == 2:
        fig.text(0.5, 0.47, "All in Cluster", ha="center", va="center")
        fig.text(0.5, 0.89, "Any in Cluster", ha="center", va="center")
    else:
        fig.text(
            0.5,
            0.91,
            f"{matches[0].capitalize()} in Cluster",
            ha="center",
            va="center",
        )
    loweraxis.set(xlabel="Divergence on Sequence Identity")
    loweraxis.legend(loc="upper left")
    fig.suptitle(title)
    plt.xscale("log")
    limits = [0.001, 1.0]
    new_tick_locations = np.array([0.0, 1.0 / 3.0, 2.0 / 3.0, 1.0])
    # loweratick_xis.set_xlim(limits)
    axis_dict["second"] = loweraxis.twiny()
    axis_dict["second"].set_xlim(limits)
    axis_dict["second"].set_xticks(new_tick_locations)
    axis_dict["second"].set_xticklabels(tick_function(new_tick_locations))
    axis_dict["second"].set_xlabel("   ")
    # r'%Identity')
    # plt.ylim([-0.002,0.002])
    outfilename = outfilestem + f"{FILETYPE}"
    print(f"saving plot to {outfilename}")
    plt.savefig(dirpath / outfilename, dpi=200)
    plt.show()


def do_cuts(obs, high, low, label):
    """Cut at high and low levels."""
    if high > 0.0:
        hicuts = obs[obs > high]
        obs = obs[obs <= high]
        if len(hicuts) == 0:
            hifilename = label + "_hicuts.tsv"
            logger.info(
                "%d observations dropped by high-side cutoff of %.2f written"
                " to %s",
                len(hicuts),
                high,
                hifilename,
            )
            logger.info(hicuts)
    if low > 0.0:
        locuts = obs[obs > low]
        obs = obs[obs >= low]
        logger.info(
            "%d observations dropped by low-side cutoff of %.2f",
            len(locuts),
            low,
        )
        if len(locuts) == 0:
            lofilename = label + "_locuts.tsv"
            logger.info(
                "%d observations dropped by low-side cutoff of %.2f written"
                " to %s",
                len(locuts),
                low,
                lofilename,
            )
            logger.info(locuts)
    return obs


@cli.command()
@click_loguru.init_logger(logfile=False)
@click.option(
    "--hi_cutoff",
    default=2.0,
    show_default=True,
    help="Disregard above this value.",
)
@click.option(
    "--lo_cutoff",
    default=0.0,
    show_default=True,
    help="Disregard below this value.",
)
@click.argument("cluster_size")
@click.argument("combinedfile")
def outlier_length_dist(hi_cutoff, lo_cutoff, cluster_size, combinedfile):
    """Plot length distribution of outliers in clusters."""
    cluster_size = int(cluster_size)
    if cluster_size <= 0:
        logger.error("Positive cluster size must be specified")
        sys.exit(1)
    clusters = pd.read_csv(combinedfile, sep="\t", index_col=0)
    norm_lengths = []
    for unused_cluster_id, cluster in clusters.groupby(
        "hom.cluster"
    ):  # pylint: disable=unused-variable
        if cluster["siz"].iloc[0] != cluster_size:
            # not the right size
            continue
        if len(set(cluster["sub"])) != 2:
            # not just two subclusters
            continue
        if 1 not in set(cluster["sub_siz"]):
            # no singleton cluster
            continue
        singleton = cluster[cluster["sub_siz"] == 1]
        length = singleton["norm"].iloc[0]
        norm_lengths.append(length)
    norm_lengths = np.array(norm_lengths)
    norm_lengths = do_cuts(norm_lengths, hi_cutoff, lo_cutoff, "len")
    logger.info(
        "%d singleton outliers in clusters of size %d",
        len(norm_lengths),
        cluster_size,
    )
    logger.info("min:\t%.3f", min(norm_lengths))
    logger.info("maxes:\t%.3f", max(norm_lengths))
    logger.info("mean: %.3f", norm_lengths.mean())
    axis = sns.distplot(norm_lengths, bins=100, kde_kws={"label": "KDE"})
    axis.set_xlabel("Normalized Length of Singleton")
    plt.title(
        "Length distribution of %d singleton subclusters" % (len(norm_lengths))
    )
    outfilename = f"norm_len_dist.{FILETYPE}"
    logger.info("saving plot to %s", outfilename)
    # plt.yscale('log')
    plt.savefig(outfilename, dpi=200)
    plt.show()


@cli.command()
@click_loguru.init_logger(logfile=False)
@click.option(
    "--hi_cutoff",
    default=0.0,
    show_default=True,
    help="Disregard above this value.",
)
@click.option(
    "--lo_cutoff",
    default=0.0,
    show_default=True,
    help="Disregard below this value.",
)
@click.argument("cluster_size")
@click.argument("combinedfile")
def length_std_dist(cluster_size, hi_cutoff, lo_cutoff, combinedfile):
    """Plot length distribution of cluster singletons."""
    cluster_size = int(cluster_size)
    if cluster_size <= 0:
        logger.error("Positive cluster size must be specified")
        sys.exit(1)
    clusters = pd.read_csv(combinedfile, sep="\t", index_col=0)
    stds = []
    for unused_cluster_id, cluster in clusters.groupby(
        "hom.cluster"
    ):  # pylint: disable=unused-variable
        if cluster["siz"].iloc[0] != cluster_size:
            # not the right size
            continue
        if len(set(cluster["sub"])) != 1:
            # Only one subcluster
            continue
        val = cluster["std"].iloc[0]
        stds.append(val)
    stds = np.array(stds)
    pct_zeros = len(stds[stds == 0.0]) * 100 / len(stds)
    stds = do_cuts(stds, hi_cutoff, lo_cutoff, "stds")
    logger.info(
        "%d single-subgroup clusters of size %d", len(stds), cluster_size
    )
    logger.info("%.1f %% zeroes, max is %.2f", pct_zeros, max(stds))
    logger.info("mean is %.3f", stds.mean())
    logbins = np.logspace(0.7, 3, 100)
    axis = sns.distplot(stds, bins=logbins, kde=False)
    axis.set_xlabel("Standard Deviation of Single-Subgroup Clusters")
    title = "Length Standard Deviation distribution of %d clusters" % len(stds)
    plt.title(title)
    outfilename = f"std_dist.{FILETYPE}"
    logger.info("saving plot to %s", outfilename)
    plt.yscale("log")
    plt.xscale("log")
    plt.savefig(outfilename, dpi=200)
    plt.show()


@cli.command()
@click_loguru.init_logger(logfile=False)
@click.option(
    "--hi_cutoff",
    default=1000,
    show_default=True,
    help="Disregard above this value.",
)
@click.option(
    "--lo_cutoff",
    default=1,
    show_default=True,
    help="Disregard below this value.",
)
@click.option(
    "--identity",
    "-i",
    default=0.0,
    help="Minimum sequence ID (0-1). [default: lowest]",
)
@click.argument("setsize", type=click.IntRange(2, 100))
@click.argument("setname")
def plot_degree_dists(identity, hi_cutoff, lo_cutoff, setname, setsize):
    """Plot homology and synteny degree distributions."""
    set_path = Path(setname)
    homo_cluster_name = cluster_set_name(setname, identity)
    homo_degree = pd.read_csv(
        set_path / homo_degree_dist_filename(homo_cluster_name),
        index_col=0,
        sep="\t",
    )
    homo_degree.index.name = "size"
    homo_degree = homo_degree.rename(
        columns={"clusters": "clusts", "pct_total": "%clusts"}
    )
    homo_clusts = sum(homo_degree["clusts"])
    homo_genes = sum(homo_degree["clusts"] * homo_degree.index)

    homo_degree["%genes"] = (
        homo_degree["clusts"] * homo_degree.index * 100.0 / homo_genes
    )
    synteny_degree_path = set_path / "dagchainer" / "clusters-sizedist.tsv"
    synteny_degree = pd.read_csv(synteny_degree_path, index_col=0, sep="\t")
    synteny_clusts = sum(synteny_degree["clusts"])
    synteny_genes = sum(synteny_degree["clusts"] * synteny_degree.index)
    logger.info("  method       clusters  genes")
    logger.info(f"homology:\t{homo_clusts}\t{homo_genes}")
    logger.info(f" synteny:\t{synteny_clusts}\t{synteny_genes}")
    # Make plot
    plt.plot(homo_degree.index, homo_degree["%genes"], label="homology")
    plt.plot(synteny_degree.index, synteny_degree["%genes"], label="synteny")
    plt.style.use("seaborn-whitegrid")
    plt.xlabel("Cluster Size")
    plt.ylabel("% of Genes in Cluster")
    plt.title(
        f"Cluster size distribution of {homo_genes} genes in {setsize}"
        f" {setname} genomes"
    )
    outfilename = f"cluster_size_dist.{FILETYPE}"
    logger.info(f"saving plot to {outfilename}")
    plt.xlim([lo_cutoff, hi_cutoff])
    plt.legend()
    plt.yscale("log")
    plt.xscale("log", basex=setsize)
    plt.savefig(outfilename, dpi=200)
    plt.show()
