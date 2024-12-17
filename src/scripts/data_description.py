# todo: transfer defs from dataset_overview.ipynb
from itertools import combinations

import numpy as np  # type: ignore
import pandas as pd  # type: ignore
import plotly.express as px  # type: ignore
import plotly.graph_objects as go  # type: ignore
import seaborn as sns  # type: ignore
from matplotlib import pyplot as plt  # type: ignore
from plotly.subplots import make_subplots  # type: ignore


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


def plot_overlaps_plotly(
    df: pd.DataFrame,
    group: list = None,
    annot: bool = True,
    title: str = "Overlap of Features",
):
    """
    Plots the overlap of available data between columns using Plotly.

    Parameters:
    - df: pandas DataFrame containing the data.
    - group: List of column names to include in the overlap calculation. If None, all columns are used.
    - annot: Boolean indicating whether to annotate the heatmap cells with overlap percentages.
    - title: Title of the heatmap.

    Returns:
    - fig: Plotly Figure object representing the overlap heatmap.
    """
    if group is None:
        group = list(df.columns)

    overlap = pd.DataFrame(index=group, columns=group, dtype=float)

    for item1, item2 in combinations(group, 2):
        overlap_value = (df[item1].notna() & df[item2].notna()).mean()
        overlap.loc[item1, item2] = overlap_value
        overlap.loc[item2, item1] = overlap_value
    np.fill_diagonal(overlap.values, 1.0)

    z = overlap.values.astype(float)

    if annot:
        text = np.vectorize(lambda x: f"{x:.2%}")(z)
    else:
        text = None

    # Threshold for accurate scale
    mask = ~np.eye(len(group), dtype=bool)
    non_diag_values = z[mask]
    clim = np.percentile(non_diag_values, 75) if len(non_diag_values) > 0 else 1.0

    heatmap = go.Heatmap(
        z=z,
        x=group,
        y=group,
        colorscale="Reds",
        zmin=0,
        zmax=clim,
        text=text,
        texttemplate="%{text}" if annot else None,
        textfont={"size": 10},
        hovertemplate="Feature X: %{x}<br>Feature Y: %{y}<br>Overlap: %{z:.2%}<extra></extra>",
        showscale=True,
        colorbar=dict(title="Overlap "),
    )

    fig = go.Figure(data=[heatmap])

    fig.update_layout(
        title="Overlap of Binding Kinetics Features",
        xaxis=dict(tickangle=45, autorange="reversed"),
        yaxis=dict(autorange="reversed"),
        template="plotly_dark",
        title_x=0.5,
        plot_bgcolor="rgb(34, 37, 41)",
        paper_bgcolor="rgb(34, 37, 41)",
        # width=800,
        # height=800
    )

    return fig


def availability_per_group(df: pd.DataFrame, group: list) -> dict:
    """Calculates the percentage of available data (non-NA) for each column."""
    percent_available = {col: df[col].notna().mean() for col in group}
    return percent_available


def plot_availability_plotly(df: pd.DataFrame, group: list, step=100):
    """
    Plots an aggregated heatmap of data availability using Plotly with a dark theme and a subtle colorscale.

    Parameters:
    - df: pandas DataFrame containing the data.
    - group: List of column names to include.
    - step: Number of rows to aggregate for density reduction.
    """

    binary_matrix = df[group].notna().astype(int)
    binary_matrix_aggregated = binary_matrix.groupby(
        binary_matrix.index // step
    ).max()  # only option to visualize with plotly

    percentages = availability_per_group(df, group)

    colorscale = [[0.0, "rgb(34, 37, 41)"], [1.0, "rgb(77, 166, 255)"]]

    heatmap = go.Heatmap(
        z=binary_matrix_aggregated.values,
        x=group,
        y=binary_matrix_aggregated.index.astype(str),
        colorscale=colorscale,
        hovertemplate="Feature: %{x}<br>Row: %{y}<br>Available: %{z}<extra></extra>",
        zmin=0,
        zmax=1,
        showscale=False,
    )

    fig = go.Figure(data=[heatmap])

    fig.update_yaxes(
        autorange="reversed",
        showgrid=False,
        tickfont=dict(color="white"),
        showticklabels=False,
    )
    fig.update_xaxes(showgrid=False, tickfont=dict(color="white"))

    annotations = []
    for i, col in enumerate(group):
        annotations.append(
            dict(
                x=col,
                y=1.03,  # Slightly above the plot
                xref="x",
                yref="paper",
                text="{:.2%}".format(percentages[col]),
                showarrow=False,
                xanchor="center",
                yanchor="top",
                font=dict(size=12, color="white"),
            )
        )

    fig.update_layout(
        title="Availability Matrix",
        xaxis_title="",
        yaxis_title="Observations",
        annotations=annotations,
        margin=dict(l=60, r=20, t=60, b=100),
        plot_bgcolor="rgb(34, 37, 41)",
        paper_bgcolor="rgb(34, 37, 41)",
        # width=800,
        # height=800,
        title_x=0.5,
        template="plotly_dark",
        font=dict(color="white"),
    )

    return fig


def categorical_countplot_plotly(
    df: pd.DataFrame,
    category: str,
    N: int = 50,
    percentile: float | None = None,
    x_scale: str = "log",
):
    """Plots the distribution of a categorical data of a dataframe column of interest using Plotly.
    Args:
        df (pd.DataFrame): dataframe containing the data of interest.
        category (str): column name of df.
        N (int, optional): first N lines of the df column (category) are plotted. Defaults to 50.
        percentile (float | None, optional): percentile (between 0 and 1) to mark on the plot. If None, no percentile line is plotted. Defaults to None.
        x_scale (str, optional): scale for x-axis. Defaults to "log".
    """
    # Get the top N categories and their counts
    top = df[category].value_counts().head(N).reset_index()
    top.columns = [category, "count"]

    # Create the bar plot
    fig = px.bar(
        top,
        x="count",
        y=category,
        orientation="h",
        title="Distribution of Researched Target Organisms",
        labels={"count": "Count", category: category},
    )

    if x_scale == "log":
        max_count = top["count"].max()
        min_count = top["count"].min()
        ticks = [
            10**i
            for i in range(int(np.log10(min_count)), int(np.log10(max_count)) + 1)
        ]
        tick_texts = [f"10^{int(np.log10(t))}" for t in ticks]

        fig.update_xaxes(
            type="log",
            tickvals=ticks,
            ticktext=tick_texts,
            title_text="Count",
        )
    elif x_scale == "linear":
        fig.update_xaxes(
            type="linear",
            title_text="Count",
        )
    else:
        raise UserWarning(
            "x_scale must be 'linear' or 'log'. Reverting to default: 'log'"
        )

    if percentile is not None:
        percentile_value = df[category].value_counts().quantile(percentile)
        fig.add_vline(
            x=percentile_value,
            line=dict(color="darkred", dash="dash"),
            annotation_text=f"{int(percentile * 100)}th Percentile",
            annotation_position="top right",
        )
        fig.add_trace(
            go.Scatter(
                x=[None],  # No points to display
                y=[None],
                mode="lines",
                line=dict(color="darkred", dash="dash"),
                name=f"{int(percentile * 100)}th Percentile",
            )
        )

        # Add percentile annotation if applicable
        annotations = [
            dict(
                x=percentile_value,
                y=1,
                xref="x",
                yref="paper",
                text=f"{int(percentile * 100)}th Percentile",
                showarrow=False,
                font=dict(color="darkred"),
            )
        ]
    else:
        annotations = []

    # Update layout for better readability
    fig.update_layout(
        xaxis_title="Count",
        yaxis_title=category,
        yaxis=dict(autorange="reversed"),
        plot_bgcolor="rgb(34, 37, 41)",
        paper_bgcolor="rgb(34, 37, 41)",
        # width=800,
        # height=600,
        title_x=0.5,
        template="plotly_dark",
        font=dict(color="white"),
        legend=dict(
            title="",
            itemsizing="constant",
            orientation="h",
            yanchor="bottom",
            y=-0.2,
            xanchor="right",
            x=1,
        ),
        annotations=annotations,
    )

    return fig


def plot_distributions_plotly(
    df: pd.DataFrame,
    variable: str,
    bins: int = 25,
    log_transform: bool = True,
    show_legend: bool = True,
    title: str = "Distributions",
):
    if variable not in df.columns:
        raise KeyError(f"Column '{variable}' not found in DataFrame")

    data = df[variable].dropna()
    if log_transform:
        data = data[data > 0]
        transformed_data = np.log10(data)
    else:
        transformed_data = data

    stats = data.describe()

    fig = go.Figure()
    fig.add_trace(
        go.Histogram(
            x=transformed_data,
            nbinsx=bins,
            marker_color="slategray",
            opacity=0.75,
            name="Distribution",
        )
    )

    counts, bin_edges = np.histogram(transformed_data, bins=bins)
    y_max = counts.max() if len(counts) > 0 else 0

    # Accurate upper bound y-axis
    if y_max > 0:
        y_upper = np.ceil(y_max / 10) * 10
        if y_upper == y_max:
            y_upper += 10
    else:
        y_upper = 10

    dtick = y_upper / 10.0
    line_top = y_upper

    q25, q50, q75 = stats["25%"], stats["50%"], stats["75%"]
    if log_transform:
        q25 = np.log10(q25) if q25 > 0 else None
        q50 = np.log10(q50) if q50 > 0 else None
        q75 = np.log10(q75) if q75 > 0 else None

    # Vertical lines for quartiles/median
    for q in [q25, q50, q75]:
        if q is not None:
            fig.add_trace(
                go.Scatter(
                    x=[q, q],
                    y=[0, line_top],
                    mode="lines",
                    line=dict(color="lightcoral", dash="dash")
                    if q in [q25, q75]
                    else dict(color="darkred", dash="dash"),
                    showlegend=False,
                )
            )

    # x-axis ticks for log scale
    if log_transform and len(transformed_data) > 0:
        min_log = np.floor(transformed_data.min())
        max_log = np.ceil(transformed_data.max())
        min_log = int(min_log)
        max_log = int(max_log)

        desired_ticks = 4
        if max_log > min_log:
            # Create evenly spaced ticks
            tickvals = np.linspace(min_log, max_log, desired_ticks)
            tickvals = np.unique([int(tv) for tv in tickvals])
            if len(tickvals) < 2:
                tickvals = [min_log, max_log]
        else:
            tickvals = [min_log, max_log]

        ticktext = [f"10^{tv}" for tv in tickvals]
        fig.update_xaxes(
            tickmode="array",
            tickvals=tickvals,
            ticktext=ticktext,
            title_text=f"log10({variable})",
            tickangle=-45,
        )
    else:
        fig.update_xaxes(title_text=variable, tickangle=-45, nticks=5)

    fig.update_yaxes(range=[0, y_upper], dtick=dtick)

    fig.update_layout(
        # height=1000, width=1000,
        title_text=title,
        title_x=0.5,
        plot_bgcolor="rgb(34, 37, 41)",
        paper_bgcolor="rgb(34, 37, 41)",
        font=dict(color="white"),
        showlegend=False,
    )

    return fig


# todo: transfer defs from dataset_overview.ipynb
from itertools import combinations

import numpy as np  # type: ignore
import pandas as pd  # type: ignore
import plotly.express as px  # type: ignore
import plotly.graph_objects as go  # type: ignore
import seaborn as sns  # type: ignore
from matplotlib import pyplot as plt  # type: ignore
from plotly.subplots import make_subplots  # type: ignore


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


def plot_overlaps_plotly(
    df: pd.DataFrame,
    group: list = None,
    annot: bool = True,
    title: str = "Overlap of Features",
):
    """
    Plots the overlap of available data between columns using Plotly.

    Parameters:
    - df: pandas DataFrame containing the data.
    - group: List of column names to include in the overlap calculation. If None, all columns are used.
    - annot: Boolean indicating whether to annotate the heatmap cells with overlap percentages.
    - title: Title of the heatmap.

    Returns:
    - fig: Plotly Figure object representing the overlap heatmap.
    """
    if group is None:
        group = list(df.columns)

    overlap = pd.DataFrame(index=group, columns=group, dtype=float)

    for item1, item2 in combinations(group, 2):
        overlap_value = (df[item1].notna() & df[item2].notna()).mean()
        overlap.loc[item1, item2] = overlap_value
        overlap.loc[item2, item1] = overlap_value
    np.fill_diagonal(overlap.values, 1.0)

    z = overlap.values.astype(float)

    if annot:
        text = np.vectorize(lambda x: f"{x:.2%}")(z)
    else:
        text = None

    # Threshold for accurate scale
    mask = ~np.eye(len(group), dtype=bool)
    non_diag_values = z[mask]
    clim = np.percentile(non_diag_values, 75) if len(non_diag_values) > 0 else 1.0

    heatmap = go.Heatmap(
        z=z,
        x=group,
        y=group,
        colorscale="Reds",
        zmin=0,
        zmax=clim,
        text=text,
        texttemplate="%{text}" if annot else None,
        textfont={"size": 10},
        hovertemplate="Feature X: %{x}<br>Feature Y: %{y}<br>Overlap: %{z:.2%}<extra></extra>",
        showscale=True,
        colorbar=dict(title="Overlap "),
    )

    fig = go.Figure(data=[heatmap])

    fig.update_layout(
        title="Overlap of Binding Kinetics Features",
        xaxis=dict(tickangle=45, autorange="reversed"),
        yaxis=dict(autorange="reversed"),
        template="plotly_dark",
        title_x=0.5,
        plot_bgcolor="rgb(34, 37, 41)",
        paper_bgcolor="rgb(34, 37, 41)",
        # width=800,
        # height=800
    )

    return fig


def availability_per_group(df: pd.DataFrame, group: list) -> dict:
    """Calculates the percentage of available data (non-NA) for each column."""
    percent_available = {col: df[col].notna().mean() for col in group}
    return percent_available


def plot_availability_plotly(df: pd.DataFrame, group: list, step=100):
    """
    Plots an aggregated heatmap of data availability using Plotly with a dark theme and a subtle colorscale.

    Parameters:
    - df: pandas DataFrame containing the data.
    - group: List of column names to include.
    - step: Number of rows to aggregate for density reduction.
    """

    binary_matrix = df[group].notna().astype(int)
    binary_matrix_aggregated = binary_matrix.groupby(
        binary_matrix.index // step
    ).max()  # only option to visualize with plotly

    percentages = availability_per_group(df, group)

    colorscale = [[0.0, "rgb(34, 37, 41)"], [1.0, "rgb(77, 166, 255)"]]

    heatmap = go.Heatmap(
        z=binary_matrix_aggregated.values,
        x=group,
        y=binary_matrix_aggregated.index.astype(str),
        colorscale=colorscale,
        hovertemplate="Feature: %{x}<br>Row: %{y}<br>Available: %{z}<extra></extra>",
        zmin=0,
        zmax=1,
        showscale=False,
    )

    fig = go.Figure(data=[heatmap])

    fig.update_yaxes(
        autorange="reversed",
        showgrid=False,
        tickfont=dict(color="white"),
        showticklabels=False,
    )
    fig.update_xaxes(showgrid=False, tickfont=dict(color="white"))

    annotations = []
    for i, col in enumerate(group):
        annotations.append(
            dict(
                x=col,
                y=1.03,  # Slightly above the plot
                xref="x",
                yref="paper",
                text="{:.2%}".format(percentages[col]),
                showarrow=False,
                xanchor="center",
                yanchor="top",
                font=dict(size=12, color="white"),
            )
        )

    fig.update_layout(
        title="Availability Matrix",
        xaxis_title="",
        yaxis_title="Observations",
        annotations=annotations,
        margin=dict(l=60, r=20, t=60, b=100),
        plot_bgcolor="rgb(34, 37, 41)",
        paper_bgcolor="rgb(34, 37, 41)",
        # width=800,
        # height=800,
        title_x=0.5,
        template="plotly_dark",
        font=dict(color="white"),
    )

    return fig


def categorical_countplot_plotly(
    df: pd.DataFrame,
    category: str,
    N: int = 50,
    percentile: float | None = None,
    x_scale: str = "log",
):
    """Plots the distribution of a categorical data of a dataframe column of interest using Plotly.
    Args:
        df (pd.DataFrame): dataframe containing the data of interest.
        category (str): column name of df.
        N (int, optional): first N lines of the df column (category) are plotted. Defaults to 50.
        percentile (float | None, optional): percentile (between 0 and 1) to mark on the plot. If None, no percentile line is plotted. Defaults to None.
        x_scale (str, optional): scale for x-axis. Defaults to "log".
    """
    # Get the top N categories and their counts
    top = df[category].value_counts().head(N).reset_index()
    top.columns = [category, "count"]

    # Create the bar plot
    fig = px.bar(
        top,
        x="count",
        y=category,
        orientation="h",
        title="Distribution of Researched Target Organisms",
        labels={"count": "Count", category: category},
    )

    if x_scale == "log":
        max_count = top["count"].max()
        min_count = top["count"].min()
        ticks = [
            10**i
            for i in range(int(np.log10(min_count)), int(np.log10(max_count)) + 1)
        ]
        tick_texts = [f"10^{int(np.log10(t))}" for t in ticks]

        fig.update_xaxes(
            type="log",
            tickvals=ticks,
            ticktext=tick_texts,
            title_text="Count",
        )
    elif x_scale == "linear":
        fig.update_xaxes(
            type="linear",
            title_text="Count",
        )
    else:
        raise UserWarning(
            "x_scale must be 'linear' or 'log'. Reverting to default: 'log'"
        )

    if percentile is not None:
        percentile_value = df[category].value_counts().quantile(percentile)
        fig.add_vline(
            x=percentile_value,
            line=dict(color="darkred", dash="dash"),
            annotation_text=f"{int(percentile * 100)}th Percentile",
            annotation_position="top right",
        )
        fig.add_trace(
            go.Scatter(
                x=[None],  # No points to display
                y=[None],
                mode="lines",
                line=dict(color="darkred", dash="dash"),
                name=f"{int(percentile * 100)}th Percentile",
            )
        )

        # Add percentile annotation if applicable
        annotations = [
            dict(
                x=percentile_value,
                y=1,
                xref="x",
                yref="paper",
                text=f"{int(percentile * 100)}th Percentile",
                showarrow=False,
                font=dict(color="darkred"),
            )
        ]
    else:
        annotations = []

    # Update layout for better readability
    fig.update_layout(
        xaxis_title="Count",
        yaxis_title=category,
        yaxis=dict(autorange="reversed"),
        plot_bgcolor="rgb(34, 37, 41)",
        paper_bgcolor="rgb(34, 37, 41)",
        # width=800,
        # height=600,
        title_x=0.5,
        template="plotly_dark",
        font=dict(color="white"),
        legend=dict(
            title="",
            itemsizing="constant",
            orientation="h",
            yanchor="bottom",
            y=-0.2,
            xanchor="right",
            x=1,
        ),
        annotations=annotations,
    )

    return fig


def availability_per_group(df: pd.DataFrame, group: list) -> dict:
    """Calculates the percentage of available data (non-NA) for each column."""
    percent_available = {col: df[col].notna().mean() for col in group}
    return percent_available


def plot_availability_plotly(df: pd.DataFrame, group: list, step=100):
    """
    Plots an aggregated heatmap of data availability using Plotly with a dark theme and a subtle colorscale.

    Parameters:
    - df: pandas DataFrame containing the data.
    - group: List of column names to include.
    - step: Number of rows to aggregate for density reduction.
    """

    binary_matrix = df[group].notna().astype(int)
    binary_matrix_aggregated = binary_matrix.groupby(
        binary_matrix.index // step
    ).max()  # only option to visualize with plotly

    percentages = availability_per_group(df, group)

    colorscale = [[0.0, "rgb(34, 37, 41)"], [1.0, "rgb(77, 166, 255)"]]

    heatmap = go.Heatmap(
        z=binary_matrix_aggregated.values,
        x=group,
        y=binary_matrix_aggregated.index.astype(str),
        colorscale=colorscale,
        hovertemplate="Feature: %{x}<br>Row: %{y}<br>Available: %{z}<extra></extra>",
        zmin=0,
        zmax=1,
        showscale=False,
    )

    fig = go.Figure(data=[heatmap])

    fig.update_yaxes(
        autorange="reversed",
        showgrid=False,
        tickfont=dict(color="white"),
        showticklabels=False,
    )
    fig.update_xaxes(showgrid=False, tickfont=dict(color="white"))

    annotations = []
    for i, col in enumerate(group):
        annotations.append(
            dict(
                x=col,
                y=1.03,  # Slightly above the plot
                xref="x",
                yref="paper",
                text="{:.2%}".format(percentages[col]),
                showarrow=False,
                xanchor="center",
                yanchor="top",
                font=dict(size=12, color="white"),
            )
        )

    fig.update_layout(
        title="Availability Matrix",
        xaxis_title="",
        yaxis_title="Observations",
        annotations=annotations,
        margin=dict(l=60, r=20, t=60, b=100),
        plot_bgcolor="rgb(34, 37, 41)",
        paper_bgcolor="rgb(34, 37, 41)",
        # width=800,
        # height=800,
        title_x=0.5,
        template="plotly_dark",
        font=dict(color="white"),
    )

    return fig


def plot_distributions_plotly(
    df: pd.DataFrame,
    variable: str,
    bins: int = 25,
    log_transform: bool = True,
    show_legend: bool = True,
    title: str = "Distributions",
):
    if variable not in df.columns:
        raise KeyError(f"Column '{variable}' not found in DataFrame")

    data = df[variable].dropna()
    if log_transform:
        data = data[data > 0]
        transformed_data = np.log10(data)
    else:
        transformed_data = data

    stats = data.describe()

    fig = go.Figure()
    fig.add_trace(
        go.Histogram(
            x=transformed_data,
            nbinsx=bins,
            marker_color="slategray",
            opacity=0.75,
            name="Distribution",
        )
    )

    counts, bin_edges = np.histogram(transformed_data, bins=bins)
    y_max = counts.max() if len(counts) > 0 else 0

    # Accurate upper bound y-axis
    if y_max > 0:
        y_upper = np.ceil(y_max / 10) * 10
        if y_upper == y_max:
            y_upper += 10
    else:
        y_upper = 10

    dtick = y_upper / 10.0
    line_top = y_upper

    q25, q50, q75 = stats["25%"], stats["50%"], stats["75%"]
    if log_transform:
        q25 = np.log10(q25) if q25 > 0 else None
        q50 = np.log10(q50) if q50 > 0 else None
        q75 = np.log10(q75) if q75 > 0 else None

    # Vertical lines for quartiles/median
    for q in [q25, q50, q75]:
        if q is not None:
            fig.add_trace(
                go.Scatter(
                    x=[q, q],
                    y=[0, line_top],
                    mode="lines",
                    line=dict(color="lightcoral", dash="dash")
                    if q in [q25, q75]
                    else dict(color="darkred", dash="dash"),
                    showlegend=False,
                )
            )

    # x-axis ticks for log scale
    if log_transform and len(transformed_data) > 0:
        min_log = np.floor(transformed_data.min())
        max_log = np.ceil(transformed_data.max())
        min_log = int(min_log)
        max_log = int(max_log)

        desired_ticks = 4
        if max_log > min_log:
            # Create evenly spaced ticks
            tickvals = np.linspace(min_log, max_log, desired_ticks)
            tickvals = np.unique([int(tv) for tv in tickvals])
            if len(tickvals) < 2:
                tickvals = [min_log, max_log]
        else:
            tickvals = [min_log, max_log]

        ticktext = [f"10^{tv}" for tv in tickvals]
        fig.update_xaxes(
            tickmode="array",
            tickvals=tickvals,
            ticktext=ticktext,
            title_text=f"log10({variable})",
            tickangle=-45,
        )
    else:
        fig.update_xaxes(title_text=variable, tickangle=-45, nticks=5)

    fig.update_yaxes(range=[0, y_upper], dtick=dtick)

    fig.update_layout(
        # height=1000, width=1000,
        title_text=title,
        title_x=0.5,
        plot_bgcolor="rgb(34, 37, 41)",
        paper_bgcolor="rgb(34, 37, 41)",
        font=dict(color="white"),
        showlegend=False,
    )

    return fig
