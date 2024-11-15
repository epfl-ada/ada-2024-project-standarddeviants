# todo: transfer defs from dataset_overview.ipynb
from itertools import combinations

import pandas as pd
import seaborn as sns
from matplotlib import pyplot as plt


def availability_per_group(df: pd.DataFrame, group: list) -> dict:
    """Caculates the percentage of available data (non-NA values) for each column in a specified group.

    Args:
        df (pd.DataFrame): dataframe containing the data to analyze.
        group (list): a list of column names of df for which we want to analyze the availability.

    Returns:
        dict: contains the column names from the group as keys and the corresponding percentages of available (non-NA) data for each.
    """
    percent_available = {key: 0 for key in group}
    for col in group:
        percent_available[col] = df[col].notna().mean()
    return percent_available


def plot_availability(df: pd.DataFrame, group: list, ax=None) -> None:
    """Plots a heatmap showing the availability (percentage of non-NA data) for a specified group of columns of df.

    Args:
        df (pd.DataFrame): dataframe containing the data of interest.
        group (list):  a list of column names of df for which we want to plot the availability.
        ax (matplotlib.axes.Axes, optional): Axes in which to draw the plot. If None, the plot is created on the current active axis. Defaults to None.
    """
    if ax:
        sns.heatmap(
            df[group].notna(),
            ax=ax,
            cbar=False,
            yticklabels=False,
        )

        percentages = availability_per_group(df, group)
        for i, col in enumerate(group):
            ax.text(
                i + 0.5,
                -0.5,
                "{:.2%}".format(percentages[col]),
                ha="center",
                va="bottom",
                fontsize=10,
                color="black",
            )


def plot_distributions(
    df: pd.DataFrame,
    variable: str,
    ax=None,
    bins: int = 25,
    log_scale: tuple[bool, bool] = (True, False),
    show_legend: bool = True,
):
    """Plot the data distribution of a variable of interest in df.

    Args:
        df (pd.DataFrame): dataframe containing the data of interest to be visualized.
        variable (str): name of the column in df to plot.
        ax (matplotlib.axes._subplots.AxesSubplot, optional): The matplotlib axis to draw the plot on. If None, the plot is created on the current active axis. Defaults to None.
        bins (int, optional): number of bins for the histogram. Defaults to 25.
        log_scale (tuple[bool, bool], optional): whether to use log axis for x-axis and y-axis. Defaults to (True, False).
        show_legend (bool, optional): whether to show the legend on the plot. Defaults to True.
    """
    stats = df[variable].describe()
    sns.histplot(
        df[variable],
        bins=bins,
        log_scale=log_scale,
        ax=ax,
        color="slategray",
        edgecolor=None,
    )
    if ax:
        ax.axvline(
            x=stats["25%"],
            color="lightcoral",
            linestyle="dashed",
            label="Lower and Upper Quartiles",
        )
        ax.axvline(x=stats["50%"], color="darkred", linestyle="dashed", label="Median")
        ax.axvline(x=stats["75%"], color="lightcoral", linestyle="dashed")
        if show_legend:
            ax.legend()
    else:
        plt.axvline(
            x=stats["25%"],
            color="lightcoral",
            linestyle="dashed",
            label="Lower and Upper Quartiles",
        )
        plt.axvline(x=stats["50%"], color="darkred", linestyle="dashed", label="Median")
        plt.axvline(x=stats["75%"], color="lightcoral", linestyle="dashed")
        if show_legend:
            plt.legend()


def plot_overlaps(
    df: pd.DataFrame,
    group: list[str] = None,
    ax=None,
    annot: bool = True,
):
    """Calculates and plots the overlap of available (non-NA) data between a specified group of columns of df.

    Args:
        df (pd.DataFrame): dataframe containing the data of interest.
        group (list[str], optional): a list of df column names for which we want to compare data availability overlap. If None, takes all df columns. Defaults to None.
        ax (matplotlib.axes.Axes, optional): Axes in which to draw the plot. If None, the plot is created on the current active axis. Defaults to None.
        annot (bool, optional): whether to show the percentage value on the plot. Defaults to True.
    """
    if not group:
        group = list(df.columns)
    pairs = combinations(group, 2)
    overlap = pd.DataFrame(
        index=group,
        columns=group.reverse(),
    )
    for item1, item2 in pairs:
        overlap.loc[item2, item1] = (df[item1].notna() & df[item2].notna()).mean()
        overlap.loc[item1, item2] = (
            df[item1].notna() & df[item2].notna()
        ).mean()  # just for symmetry

    max_overlaps = overlap.max()
    clim = max_overlaps.describe()["75%"]
    overlap = overlap.fillna(1)  # 100% overlap with self
    sns.heatmap(
        overlap,
        annot=annot,
        annot_kws={"size": 8},
        cmap="Reds",
        vmax=clim,
        fmt=".2%",
        ax=ax,
    )


def categorical_countplot(
    df: pd.DataFrame,
    category: str,
    N: int = 50,
    percentile: float | None = None,
    ax=None,
    x_scale: str = "log",
):
    """Plots the distribution of a categorical data of a dataframe column of interest.
    Args:
        df (pd.DataFrame): dataframe containing the data of interest.
        category (str): column name of df.
        N (int, optional): first N lines of the df column (category) are plotted. Defaults to 50.
        percentile (float | None, optional): percentile (between 0 and 1) to mark on the plot as a vertical dashed line. If None, no percentile line is plotted. Defaults to None.
        ax (matplotlib.axes.Axes, optional): Axes in which to draw the plot. If None, the plot is created on the current active axis. Defaults to None.
        x_scale (str, optional): scale for x-axis. Defaults to "log".
    """
    top = df[category].value_counts().head(N)
    top = pd.DataFrame(top).reset_index()
    sns.barplot(top, x="count", y=category, ax=ax, color="slategray")
    if x_scale in ["linear", "log"]:
        if ax is not None:
            ax.set_xscale(x_scale)
        else:
            plt.xscale(x_scale)
    else:
        if ax is not None:
            ax.set_xscale("log")
        else:
            plt.xscale(x_scale)
        raise UserWarning(
            "x_scale must be 'linear' or 'log'\nReverting to default: 'log'"
        )
    if ax is not None:
        ax.set_xlabel("Count")
    else:
        plt.xlabel("Count")

    if percentile is not None:
        if ax is not None:
            ax.axvline(
                x=df[category].value_counts().quantile(percentile),
                color="darkred",
                label="{:.0f}th percentile".format(100 * percentile),
                linestyle="dashed",
            )
            ax.legend()
        else:
            plt.axvline(
                x=df[category].value_counts().quantile(percentile),
                color="darkred",
                label="{:.0f}th percentile".format(100 * percentile),
                linestyle="dashed",
            )
            plt.legend()
