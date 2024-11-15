# todo: transfer defs from dataset_overview.ipynb
from itertools import combinations

import pandas as pd
import seaborn as sns
from matplotlib import pyplot as plt


def availability_per_group(df: pd.DataFrame, group: list) -> dict:
       """
    Calculate the availability of each column in a given group within a DataFrame.
    
    Parameters:
        df (pd.DataFrame): The DataFrame to analyze.
        group (list): List of column names to check availability for.
    
    Returns:
        dict: A dictionary where keys are column names and values are the fraction of non-null entries.
    """
    percent_available = {key: 0 for key in group}
    for col in group:
        percent_available[col] = df[col].notna().mean()
    return percent_available


def plot_availability(df: pd.DataFrame, group: list, ax=None) -> None:
    """
    Plot a heatmap of data availability for a specific group of columns and display availability percentages.
    
    Parameters:
        df (pd.DataFrame): The DataFrame to plot.
        group (list): List of column names to include in the heatmap.
        ax (matplotlib.axes._subplots.AxesSubplot, optional): Matplotlib axis to draw the plot on.
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
    """
    Plot the distribution of a variable with optional quartile lines and logarithmic scaling.
    
    Parameters:
        df (pd.DataFrame): The DataFrame containing the data.
        variable (str): Column name of the variable to plot.
        ax (matplotlib.axes._subplots.AxesSubplot, optional): Matplotlib axis to draw the plot on.
        bins (int): Number of bins for the histogram. Default is 25.
        log_scale (tuple[bool, bool]): Tuple indicating logarithmic scaling for x and y axes.
        show_legend (bool): Whether to show a legend for quartile lines. Default is True.
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
    group: list[str] = None,  # if none, takes all df columns
    ax=None,
    annot: bool = True,
):
    """
    Plot a heatmap of column overlaps (percentage of non-null values shared between columns) for a group.
    
    Parameters:
        df (pd.DataFrame): The DataFrame to analyze.
        group (list[str], optional): List of column names to consider for overlaps. Uses all columns if None.
        ax (matplotlib.axes._subplots.AxesSubplot, optional): Matplotlib axis to draw the plot on.
        annot (bool): Whether to annotate the heatmap cells with percentages. Default is True.
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
     """
    Plot a count plot for a categorical variable with optional percentile line and scaling.
    
    Parameters:
        df (pd.DataFrame): The DataFrame containing the data.
        category (str): Column name of the categorical variable to plot.
        N (int): Maximum number of unique categories to display. Default is 50.
        percentile (float | None): Quantile line to display (0 to 1) for the counts. If None, no line is shown.
        ax (matplotlib.axes._subplots.AxesSubplot, optional): Matplotlib axis to draw the plot on.
        x_scale (str): Scaling for the x-axis, either 'linear' or 'log'. Default is 'log'.
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
