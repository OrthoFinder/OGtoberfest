import os
import pandas as pd 
import numpy as np 
import seaborn as sns 
import matplotlib.colors as mcolors
import matplotlib.lines as mlines
import matplotlib.pyplot as plt
sns.set_theme(style="white", context="talk")


def plot_heatmap(df):

    sns.set_theme(style="white")
    corr = df.corr()
    mask = np.triu(np.ones_like(corr, dtype=bool))
    f, ax = plt.subplots(figsize=(8, 6))
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

    num_cols = len(df.columns[1:])
    ncols = 4
    nrows = int(np.ceil(num_cols / ncols))

    fig, axes = plt.subplots(
        nrows=nrows, 
        ncols=ncols, 
        # figsize=(3 * nrows, 2 * ncols), 
        sharex=False, 
        sharey=False,
        layout="constrained"
    )
    axes = axes.flatten() if num_cols > 1 else [axes]
    for ax, col in zip(axes, df.columns[1:]):
        # plot_boxplot(df, "bin_number", col, ax)
        sns.barplot(df, x="Methods", y=col, hue="Methods", palette="deep", ax=ax)
        ax.axhline(0, color="k", clip_on=False)
        ax.set_xticklabels(ax.get_xticklabels(), rotation=90, ha="right")
        ax.set_xlabel("")
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
        "Avg Rank Score",
        "Weighted Avg Recall",
        "Weighted Avg Precision",
        "Fission (RefOG)",
        "Fussion (RefOG)",
        "Avg Entropy",
        "Missing Genes",
    ],
    inplace=True,
    ascending=[True, False, False, True, True, True, True],
)
print(df)
# plot_heatmap(df.set_index("Methods"))
plot_barplot(df)
