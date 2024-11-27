import os
import pandas as pd 
import numpy as np 
import seaborn as sns 
import matplotlib.colors as mcolors
import matplotlib.lines as mlines
import matplotlib.pyplot as plt
import matplotlib.colors as mcolors
sns.set_theme(style="white", context="talk")

colors = set([
    "goldenrod",
    "olive",
    "indigo",
    "darkgreen",
    "lightblue",
    "blue",
    # "orange",
    "yellowgreen",
    "plum",
    "purple",
    "silver",
    "tomato",
    "turquoise",
    "pink",
    "brown",
])

def plot_heatmap(df):

    sns.set_theme(style="white")
    corr = df.corr()
    mask = np.triu(np.ones_like(corr, dtype=bool))
    f, ax = plt.subplots(figsize=(8, 6), layout="constrained")
    cmap = sns.diverging_palette(230, 20, as_cmap=True)

    sns.heatmap(corr, mask=mask, cmap=cmap, vmax=.3, center=0,
                square=True, linewidths=.5, cbar_kws={"shrink": .3},
                annot=True, fmt=".2f", ax=ax)
    plt.show()


def plot_boxplot(df, binned_col, col, ax):
    sns.set_theme(style="ticks")
    # f, ax = plt.subplots(figsize=(7, 6))
    # ax.set_xscale("log")

    boxplot = sns.boxplot(
        df,
        x=binned_col,
        y=col,
        hue=binned_col,
        legend=False,
        whis=[0, 100],
        width=0.6,
        palette="vlag",
        ax=ax,
    )
    sns.stripplot(df, x=binned_col, y=col, size=4, color=".3", jitter=True, ax=ax)

    means = df.groupby(binned_col, observed=True)[col].mean()
    x_positions = range(len(means))
    ax.scatter(x_positions, means, color="red", marker="s", s=50, label="Mean", zorder=5)
    ax.plot(x_positions, means, color="red", linestyle="-", linewidth=2, zorder=4)

    ax.set_xticks(x_positions)
    ax.set_xticklabels([str(x) for x in x_positions])
    ax.xaxis.grid(True)
    ax.set(ylabel="", xlabel="")
    custom_legend = mlines.Line2D([], [], linestyle="None", markersize=10, label=col)
    ax.legend(handles=[custom_legend], loc="best", frameon=True)
    sns.despine(trim=True, left=True)


def plot_barplot(df):

    methods = df["Methods"].unique()
    color_palette = {}
    for method in methods:
        color_palette[method] = colors.pop()
        print(method, color_palette[method])

    num_cols = len(df.columns[1:])
    ncols = 5
    nrows = int(np.ceil(num_cols / ncols))

    fig, axes = plt.subplots(
        nrows=nrows, 
        ncols=ncols, 
        # figsize=(4 * nrows, 3 * ncols), 
        sharex=False, 
        sharey=False,
        layout="constrained"
    )
    axes = axes.flatten() if num_cols > 1 else [axes]
    for ax, col in zip(axes, df.columns[1:]):
        # plot_boxplot(df, "bin_number", col, ax)
        if col in ["RefOG Fissions (%)", "Missing Genes (%)", "Missing Species (%)", "Entropy", "Effective Size", "Runtime", "Rank Score"]:
            df.sort_values(by=[col], inplace=True)
        elif col in ["Recall", "Precision",]:
            df.sort_values(by=[col], inplace=True, ascending=False)

        else:
            df.sort_values(by=["Rank Score"], inplace=True)
        sns.barplot(df, x="Methods", 
                    y=col, 
                    hue="Methods", 
                    # palette="deep", 
                    palette=color_palette,
                    # palette="Paired",
                    ax=ax,
                    width=1.0
                    )
        ax.axhline(0, color="k", clip_on=False)
        ax.set_xticks(range(len(df["Methods"])))
        ax.set_xticklabels(df["Methods"], rotation=90, ha="center", size=13)
        ax.set_xlabel("")
        ax.set_ylabel(col, size=13)
    for i in range(num_cols, len(axes)):
        fig.delaxes(axes[i])

    # plt.tight_layout()
    # output_path = os.path.join(d_output, f"{args.param}_boxplot_grid.png")
    # plt.savefig(output_path, dpi=300, bbox_inches="tight")
    plt.show()


input_path = r"./OrthoBench/scores_preprocessed_predicted_orthogroups/OrthoBench_global_scores.tsv"

df = pd.read_csv(input_path, sep="\t")
df.sort_values(
    by=[
        "Rank Score",
        "Recall",
        "Precision",
        "RefOG Fissions (%)",
        "RefOG Fusions (%)",
        "Entropy",
        "Missing Genes (%)",
        "Missing RefOGs (%)",
        "Runtime",
        "Effective Size",
    ],
    inplace=True,
    ascending=[
        True, 
        False, 
        False, 
        True, 
        True, 
        True, 
        True, 
        True, 
        True, 
        True
    ],
)

print(df)
plot_heatmap(df.set_index("Methods"))
plot_barplot(df)
