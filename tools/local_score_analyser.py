import pathlib
import os
import functools
import pandas as pd
import numpy as np
import seaborn as sns
import matplotlib.colors as mcolors
import matplotlib.lines as mlines
import matplotlib.pyplot as plt

# sns.set_theme(style="ticks")


def plot_barplot(df, colnames):

    num_cols = len(colnames)
    ncols = 4
    nrows = int(np.ceil(num_cols / ncols))

    fig, axes = plt.subplots(
        nrows=nrows,
        ncols=ncols,
        # figsize=(3 * nrows, 2 * ncols),
        sharex=False,
        sharey=False,
        layout="constrained",
    )
    axes = axes.flatten() if num_cols > 1 else [axes]
    for ax, col in zip(axes, colnames):
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



def plot_boxplot(df, binned_col, col, ax):
    sns.set_theme(style="ticks")
    boxplot = sns.boxplot(
        df, x=binned_col, y=col, hue=binned_col, legend=False,
        whis=[0, 100], width=.6, palette="vlag", ax=ax
    )
    sns.stripplot(df, x=binned_col, y=col, size=4, color=".3", jitter=True, ax=ax)

    means = df.groupby(binned_col, observed=True)[col].mean()
    x_positions = range(len(means))
    ax.scatter(x_positions, means, color='red', marker='s', s=50, label='Mean', zorder=5)
    ax.plot(x_positions, means, color='red', linestyle='-', linewidth=2, zorder=4)

    ax.set_xticks(x_positions)
    ax.set_xticklabels([str(x) for x in means.index])
    ax.xaxis.grid(True)
    ax.set(ylabel="", xlabel="")
    custom_legend = mlines.Line2D([], [], linestyle='None', markersize=10, label=col)
    ax.legend(handles=[custom_legend], loc='best', frameon=True)
    sns.despine(trim=True, left=True)


def corr_analysis(combined_df, title, d_output=None):
    num_cols = len(combined_df.columns[1:-2])
    ncols = 2  
    nrows = int(np.ceil(num_cols / ncols))

    fig, axes = plt.subplots(
        nrows=nrows, 
        ncols=ncols, 
        # figsize=(5 * nrows, 8 * ncols), 
        sharex=False, 
        sharey=True,
        layout="constrained",
        )
    axes = axes.flatten() if num_cols > 1 else [axes]
    for ax, col in zip(axes, combined_df.columns[1:-2]):
        plot_boxplot(combined_df, "avg_bin", col, ax)

    for i in range(num_cols, len(axes)):
        fig.delaxes(axes[i])
    fig.suptitle(title)
    # plt.tight_layout()
    # output_path = os.path.join(d_output, f"{title}.png")
    # plt.savefig(output_path, dpi=300, bbox_inches='tight') 
    plt.show()



local_score_dir = (
    r"./OrthoBench/scores_preprocessed_predicted_orthogroups/local_scores"
)
local_score_fig = (
    r"./OrthoBench/figs"
)
refogs_path = r"./OrthoBench/reference_orthogroups.tsv"

refogs_df = pd.read_csv(refogs_path, sep="\t", usecols=lambda x: x != "Orthogroups")
# print(refogs_df)

refog_colnames = [
    # "num_genes_per_species",
    # "avg_tree_length",
    # "mean_length",
    # "invariable_sites_proportion",
    "gamma_shape",
    # "total_tree_length",
]

predog_colnames = [
    # "avg_Recall (%)",
    # "avg_Precision (%)",
    "avg_F1-score (%)",
    # "EffectiveSize (JI_weighted)",
    # "MissingGenes (%)",
    # "Entropy",
]

# score_dfs_dict = {}
# for predog_col in predog_colnames:
#     score_dfs = []


for predog_col in predog_colnames:

    for eval_col in refog_colnames:
        if eval_col == "num_genes_per_species":
            refogs_df[eval_col] = 100. * refogs_df[eval_col]
        if eval_col == "gamma_shape":
            refogs_df.drop(refogs_df.tail(1).index, inplace=True)

        refogs_df["bin_col"] = pd.cut(refogs_df[eval_col], 20, duplicates="drop")
        refogs_df["avg_bin"] =  refogs_df["bin_col"].astype(str)
        bin_mean = refogs_df["avg_bin"].str.replace("\(|]", "", regex=True).str.split(", ", expand=True)
        refogs_df["avg_bin"] = bin_mean.astype(float).mean(axis=1).round(2)

        # refogs_df["bin_number"] = refogs_df["bin_col"].cat.codes + 1 

        score_dfs = [refogs_df[['RefOGs', eval_col, "avg_bin"]]]
        for file in pathlib.Path(local_score_dir).iterdir():
            # df = pd.read_csv(file, sep="\t", usecols=lambda x: x != "Missing_Genes")
            score_df = pd.read_csv(file, sep='\t', usecols=["RefOGs", predog_col])
            score_df.rename(columns={predog_col: file.name.rsplit(".", 1)[0]}, inplace=True)
            score_dfs.append(score_df)

        combined_df = functools.reduce(lambda left, right: pd.merge(left, right, on='RefOGs'), 
                                    score_dfs[::-1])
        
        combined_df.sort_values(by=[eval_col], inplace=True)
        corr_analysis(combined_df, predog_col + " vs. " + eval_col, local_score_fig)