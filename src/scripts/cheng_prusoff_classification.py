from typing import Literal

import numpy as np
import pandas as pd
import seaborn as sns
import statsmodels.formula.api as smf
from matplotlib import pyplot as plt


def cheng_prusoff_data(df: pd.DataFrame):
    """From a df containing columns for Ki (nM) and IC50 (nM), prepares data for cheng-prusoff classification. Steps include finding overlapping data and applying log-transform.

    Args:
        df (pd.DataFrame): df containing at least columns with Ki (nM) and IC50 (nM) (other columns are dropped)

    Returns:
        pd.DataFrame: prepared data
    """
    both = df[df["IC50 (nM)"].notna() & df["Ki (nM)"].notna()][["IC50 (nM)", "Ki (nM)"]]
    both_log = both.apply(np.log)
    both_log = both_log.rename(columns={"Ki (nM)": "log_Ki", "IC50 (nM)": "log_IC50"})
    return both_log


def cheng_prusoff_model(df: pd.DataFrame, random_state=11):
    """From a dataframe containing prepared data (see cheng_prusoff_data) for log_Ki and log_IC50, return a fit ols linear model for log_Ki ~ log_IC50

    Args:
        df (pd.DataFrame): prepared data containing log_Ki and log_IC50
        random_state (int, optional): random state for model. Defaults to 11.

    Returns:
        RegressionResults: fit linear regression model for log_Ki ~ log_IC50
    """
    mod = smf.ols(formula="log_Ki ~ log_IC50", data=df)
    np.random.seed(random_state)
    res = mod.fit()
    return res


def cheng_prusoff_classifier(
    df: pd.DataFrame,
    min_seperator_interecpt: float = -5,
    max_seperator_interecpt: float = 10,
    metric: Literal["abs", "square"] = "square",
    show_evaluation: bool = True,
):
    """On the basis of binding kinetics offered by the Cheng-Prusoff equation, creates a linear classifier seperating data by binding kinetics.

    Args:
        df (pd.DataFrame): dataframe containing the prepared log_Ki and log_IC50 (see cheng_prusoff_data)
        min_seperator_interecpt (float, optional): minimum intercept for the linear classifier. Defaults to -5.
        max_seperator_interecpt (float, optional): maximum intecept for the linear classifier. Defaults to 10.
        metric (Literal['abs', 'square'], optional): metric to fit classifier. Defaults to "square".
        show_evaluation (bool, optional): show the fitting evaluation. Defaults to True.

    Raises:
        ValueError: in case of incorrect metric

    Returns:
        dict: key: alpha (score balancing parameter); value: intercept maximising the metric (optimal classifier per alpha)
    """
    # storage
    b_range = []
    slope_diff_range = []
    int_diff_range = []

    # linear seperator (y=mx+b)
    m = 1  # set to 1 by loglog cheng-prussof derivation
    for b in np.arange(min_seperator_interecpt, max_seperator_interecpt, 0.1):
        # assign labels; cluster 1 = below line
        df["label"] = (df["log_IC50"] < m * df["log_Ki"] + b).astype(int)
        fam0 = df.query("label == 0")
        fam1 = df.query("label == 1")

        # model per cluster
        res0 = cheng_prusoff_model(fam0)
        res1 = cheng_prusoff_model(fam1)
        m0, b0 = res0.params
        m1, b1 = res1.params

        b_range.append(b)
        if metric == "square":
            slope_diff_range.append(np.square(m0 - m1))
            int_diff_range.append(np.square(b0 - b1))
        elif metric == "abs":
            slope_diff_range.append(np.abs(m0 - m1))
            int_diff_range.append(np.abs(b0 - b1))
        else:
            raise ValueError("Metric must be either 'square' or 'abs'")

    dists = np.array([slope_diff_range, int_diff_range]).T
    dists = (dists - np.mean(dists, axis=0)) / np.std(dists, axis=0)
    dists[
        :, 0
    ] *= (
        -1
    )  # invert slope distance bcs need to maximise intercept diff and minimise slope dist
    distance_scores = {}
    b_max = {}
    # alpha: importance of having similar slopes (tradeoff)
    # 1-alpha: importance of having different interecepts
    for alpha in np.arange(0, 1, 0.01):
        distance_score = alpha * dists[:, 0] + (1 - alpha) * dists[:, 1]
        distance_scores[alpha] = distance_score
        b_max[alpha] = b_range[np.argmax(distance_score)]

    if show_evaluation:
        fig = plt.figure()
        fig.set_figheight(10)
        fig.set_figwidth(10)
        ax1 = plt.subplot2grid(shape=(2, 2), loc=(0, 0), colspan=1)
        ax2 = plt.subplot2grid(shape=(2, 2), loc=(0, 1), colspan=1)
        ax3 = plt.subplot2grid(shape=(2, 2), loc=(1, 0), colspan=2)
        axs = [ax1, ax2, ax3]

        axs[0].plot(b_range, dists[:, 0])
        axs[0].set_title(
            "Negative Standardised Distance\nin population Slopes", fontsize=10
        )
        axs[1].plot(b_range, dists[:, 1])
        axs[1].set_title("Standardised Distance\nin population Intercepts", fontsize=10)
        i = 0
        for alpha in np.arange(0, 1, 0.01):
            if i % 10 == 0:
                (line,) = axs[2].plot(
                    b_range, distance_scores[alpha], label=rf"$\alpha$ = {alpha:.2f}"
                )
                axs[2].vlines(
                    x=b_max[alpha],
                    ymin=-4,
                    ymax=2.3,
                    color=line.get_color(),
                    linestyle="--",
                )
            i += 1
        axs[2].set_title("Total Distance Score between Population")
        plt.legend()
        plt.show()
    return b_max


def cheng_prusoff_plot(
    families: list[pd.DataFrame],
    regs: list,
    ax=None,
):
    palettes = ["Reds_r", "Blues_r"]
    colors = ["maroon", "midnightblue"]
    n_families = len(families)

    for idx, subdata in enumerate(
        families
    ):  # subdata is one cluster from the original overlapping data
        # density lines
        sns.kdeplot(
            subdata[["log_Ki", "log_IC50"]].dropna(),
            y="log_Ki",
            x="log_IC50",
            levels=7,
            cmap=palettes[idx % n_families],
            fill=False,
            linewidths=0.7,
            ax=ax,
        )

        # data points
        sns.scatterplot(
            subdata[["log_Ki", "log_IC50"]].dropna(),
            y="log_Ki",
            x="log_IC50",
            s=1,
            color=colors[idx % n_families],
            ax=ax,
        )

        # model regression
        if ax is not None:
            x = np.linspace(-10, 25, num=50)
            b, a = regs[idx].params
            adjr2 = regs[idx].rsquared_adj
            ax.plot(
                x,
                a * x + b,
                color=colors[idx % n_families],
                linestyle="--",
                label=f"{a:.2f}x + {b:.2f} ($R^2_{{adj}}$ = {adjr2:.2f})",
            )


def classified_h_index(df_cit: pd.DataFrame):
    """calculate h-index per class if dataframe comes, highly classified: column1=class of element; column2=citation count of element.
    Inteded for use after Cheng-Prusoff classification, probably applicable elsewhere.
    For more complex dataframe types, see src.scripts.citations.calculate_H_index

    Args:
        df_cit (pd.DataFrame): column1=class of element, named "cluster"; column2=citation count of element, named "citation"

    Returns:
        dict: contains h-index of each class
    """
    classes = df_cit["cluster"].unique()
    h_per_class = {}
    for c in classes:
        temp = df_cit.query("cluster==@c")
        h = sum(
            x >= i + 1 for i, x in enumerate(sorted(temp["citation"], reverse=True))
        )
        h_per_class[c] = h
    return h_per_class
