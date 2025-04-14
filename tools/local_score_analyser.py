import pathlib
import os
import functools
import pandas as pd
import numpy as np
import seaborn as sns
import matplotlib.colors as mcolors
import matplotlib.lines as mlines
import matplotlib.pyplot as plt
import matplotlib.ticker as mticker
from scipy.spatial.distance import pdist, squareform
from scipy.stats import  wasserstein_distance
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



def plot_boxplot(df, x_col, col, ax, plot_average = True):
    sns.set_theme(style="ticks")
    boxplot = sns.boxplot(
        df, x=x_col, y=col, hue=x_col, legend=False,
        whis=[0, 100], width=.6, palette="vlag", ax=ax
    )
    sns.stripplot(df, x=x_col, y=col, size=4, color=".3", jitter=True, ax=ax)

    means = df.groupby(x_col, observed=True)[col].mean()
    x_positions = range(len(means))
    ax.scatter(x_positions, means, color='red', marker='s', s=50, label='Mean', zorder=5)
    if plot_average:
        ax.plot(x_positions, means, color='red', linestyle='-', linewidth=2, zorder=4)

    ax.set_xticks(x_positions)
    ax.set_xticklabels([str(x) for x in means.index])
    ax.xaxis.grid(True)

    # for spine_name in ['left', 'bottom', 'right', 'top']:
    #     ax.spines[spine_name].set_visible(True)  # Make all spines visible
    #     ax.spines[spine_name].set_linewidth(0.8)  # Optional: Adjust thickness
    #     ax.spines[spine_name].set_color('black')  # Optional: Set color to black


    # ax.spines['top'].set_visible(True)
    # ax.spines['right'].set_visible(True)
    # ax.spines['bottom'].set_visible(True)
    # ax.spines['left'].set_visible(True)

    ax.set(ylabel="", xlabel="")
    custom_legend = mlines.Line2D([], [], linestyle='None', markersize=10, label=col)
    ax.legend(handles=[custom_legend], loc='best', frameon=True, fontsize=14)
    sns.despine(trim=True, left=True)


def corr_analysis(combined_df, title, d_output=None):
    num_cols = len(combined_df.columns[1:-2])
    ncols = 3 
    nrows = int(np.ceil(num_cols / ncols))

    fig, axes = plt.subplots(
        nrows=nrows, 
        ncols=ncols, 
        # figsize=(5 * nrows, 8 * ncols), 
        sharex=False, 
        sharey=False,
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



def local_score_display(combined_df, x_col, title, d_output=None):
    num_cols = len(combined_df.columns[1:-1])
    ncols = 3 
    nrows = int(np.ceil(num_cols / ncols))

    fig, axes = plt.subplots(
        nrows=nrows, 
        ncols=ncols, 
        # figsize=(5 * nrows, 8 * ncols), 
        sharex=False, 
        sharey=False,
        layout="constrained",
        )
    x_positions = range(combined_df.shape[0])
    axes = axes.flatten() if num_cols > 1 else [axes]
    N = len(x_positions)
    # cmap = plt.get_cmap("viridis", N)
    # colors = [cmap(i) for i in range(N)]
    palette = sns.color_palette("deep", n_colors=N)

    # combined_df[x_col] = combined_df[x_col].astype(str)
    for ax, col in zip(axes, combined_df.columns[1:-1]):
        y_vals = combined_df[col].values
        ax.bar(x_positions, y_vals, color=palette[0])
        # sns.barplot(combined_df, x=x_col, y=col, ax=ax, errorbar=None)
        # ax.yaxis.set_major_locator(mticker.MaxNLocator(integer=True))
        # ax.axhline(0, color="k", clip_on=False)
        ax.set_xticks(x_positions)
        ax.set_xticklabels([str(round(combined_df[x_col].iloc[x], 2)) for x in x_positions])
        ax.xaxis.set_major_locator(mticker.MaxNLocator(10))
        ax.set(ylabel="", xlabel="")
        custom_legend = mlines.Line2D([], [], linestyle='None', markersize=10, label=col)
        ax.legend(handles=[custom_legend], loc='best', frameon=True, fontsize=14)

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
local_property_dir = r"./OrthoBench/scores_preprocessed_predicted_orthogroups/local_properties"
local_score_fig = (
    r"./OrthoBench/figs"
)
refogs_path = r"./OrthoBench/reference_orthogroups_stats.tsv"

refogs_df = pd.read_csv(refogs_path, sep="\t", usecols=lambda x: x != "Orthogroups")
# print(refogs_df)

refog_colnames = [
    "num_genes",
    # "num_genes_per_species",
    # "avg_tree_length",
    # "mean_length",
    # "invariable_sites_proportion",
    # "gamma_shape_alpha",
    # "total_tree_length",
]


predog_colnames = [
    # "Weighted Avg Recall",
    # "Weighted Avg Precision",
    # "Weighted Avg F1-score",
    # "Weighted Avg Fowlkes-Mallows Index",
    # "Weighted Avg Jaccard Index",
    # "Weighted Avg Dissimilarity",
    # "Missing Genes Percentage",
    # "Missing Genes Percentage",
    "Entropy",
]
property_colnames = [
    # "Missing Species Count",
    # "Missing Genes Count",
    "Effective Size",	
    # "Fusion Bool",	
    # "Fission Bool",
]

# score_dfs_dict = {}
# for predog_col in predog_colnames:
#     score_dfs = []

## -------------------------Corr Analysis -----------------------------------------

# for predog_col in predog_colnames:
    
#     for eval_col in refog_colnames:
#         if eval_col == "num_genes_per_species":
#             refogs_df[eval_col] = 100. * refogs_df[eval_col]
#         if eval_col == "gamma_shape":
#             refogs_df.drop(refogs_df.tail(1).index, inplace=True)
#         print(f"Processing {predog_col} vs. {eval_col}....")
#         refogs_df["bin_col"] = pd.cut(refogs_df[eval_col], 20, duplicates="drop")
#         refogs_df["avg_bin"] =  refogs_df["bin_col"].astype(str)
#         bin_mean = refogs_df["avg_bin"].str.replace("\(|]", "", regex=True).str.split(", ", expand=True)
#         refogs_df["avg_bin"] = bin_mean.astype(float).mean(axis=1).round(1)

#         # refogs_df["bin_number"] = refogs_df["bin_col"].cat.codes + 1 

#         score_dfs = [refogs_df[['RefOGs', eval_col, "avg_bin"]]]
#         for file in pathlib.Path(local_score_dir).iterdir():
            
#             if "global_score" in file.name or not file.is_file():
#                 continue
#             # df = pd.read_csv(file, sep="\t", usecols=lambda x: x != "Missing_Genes")
#             score_df = pd.read_csv(file, sep='\t', usecols=["RefOGs", predog_col])
#             score_df.rename(columns={predog_col: file.name.rsplit(".", 1)[0]}, inplace=True)
#             score_dfs.append(score_df)

#         combined_df = functools.reduce(lambda left, right: pd.merge(left, right, on='RefOGs'), 
#                                     score_dfs[::-1])
        
#         combined_df.sort_values(by=[eval_col], inplace=True)
#         corr_analysis(combined_df, eval_col + " vs. " + predog_col, local_score_fig)



## --------------------------------------


# for eval_col in refog_colnames:
#     if eval_col == "num_genes_per_species":
#         refogs_df[eval_col] = refogs_df[eval_col]
#     if eval_col == "gamma_shape":
#         refogs_df.drop(refogs_df.tail(1).index, inplace=True)
    
#     score_dfs = [refogs_df[['RefOGs', eval_col]]]
#     # for predog_col in predog_colnames:
#     for predog_col in property_colnames:
#         print(f"Processing {eval_col} vs. {predog_col}....")
#         # for file in pathlib.Path(local_score_dir).iterdir():
#         for file in pathlib.Path(local_property_dir).iterdir():
            
#             if "global_score" in file.name or not file.is_file():
#                 continue
#             # df = pd.read_csv(file, sep="\t", usecols=lambda x: x != "Missing_Genes")
#             score_df = pd.read_csv(file, sep='\t', usecols=["RefOGs", predog_col])
#             score_df.rename(columns={predog_col: file.name.rsplit(".", 1)[0]}, inplace=True)
#             score_dfs.append(score_df)

#         combined_df = functools.reduce(lambda left, right: pd.merge(left, right, on='RefOGs'), 
#                                     score_dfs[::-1])
        
#         combined_df.sort_values(by=[eval_col], inplace=True)
#         # print(combined_df.info())
#         local_score_display(combined_df, eval_col, eval_col + " vs. " + predog_col, local_score_fig)
        
## -----------------------------------------


    
for eval_col in refog_colnames:
    if eval_col == "num_genes_per_species":
        refogs_df[eval_col] = 100. * refogs_df[eval_col]
    if eval_col == "gamma_shape":
        refogs_df.drop(refogs_df.tail(1).index, inplace=True)
    score_dfs = [refogs_df[["RefOGs", eval_col]]]
    methods = ["RefOGs Size"]
    for predog_col in property_colnames:
        print(f"Processing {predog_col} vs. {eval_col}....")

        # refogs_df["bin_col"] = pd.cut(refogs_df[eval_col], 20, duplicates="drop")
        # refogs_df["avg_bin"] =  refogs_df["bin_col"].astype(str)
        # bin_mean = refogs_df["avg_bin"].str.replace("\(|]", "", regex=True).str.split(", ", expand=True)
        # refogs_df["avg_bin"] = bin_mean.astype(float).mean(axis=1).round(2)

        # refogs_df["bin_number"] = refogs_df["bin_col"].cat.codes + 1 

        for file in pathlib.Path(local_property_dir).iterdir():
            methods.append(file.name.split(".", 1)[0])
            # df = pd.read_csv(file, sep="\t", usecols=lambda x: x != "Missing_Genes")
            score_df = pd.read_csv(file, sep='\t', usecols=["RefOGs", predog_col])
            score_df.rename(columns={predog_col: file.name.rsplit(".", 1)[0]}, inplace=True)
            score_dfs.append(score_df)

        combined_df = functools.reduce(lambda left, right: pd.merge(left, right, on='RefOGs'), 
                                    score_dfs)
        
        combined_df.sort_values(by=[eval_col], inplace=True)

        combined_df.rename(columns={"num_genes": "RefOGs Size"}, inplace=True)
        # combined_df.drop(columns=[eval_col], inplace=True)
        combined_df = combined_df.iloc[:, 1:].astype(int)
        print(combined_df)
        combined_arr = combined_df.to_numpy().T

        dist_arr = pdist(
            combined_arr, 
            metric=lambda u, v: wasserstein_distance(u, v)
        )

        # Optionally, convert the condensed distance array to a square matrix
        # dist_matrix = squareform(dist_arr)

        # dist_arr = pdist(combined_arr, metric='jensenshannon')
        m = combined_arr.shape[0]
        dist_dict = {}
        i = 0
        for j in range(m):
            print(methods[j])
            dist_dict[methods[j]] = dist_arr[m * i + j - ((i + 2) * (i + 1)) // 2].round(2)
        # melted_df = combined_df.melt(id_vars="num_genes", var_name="Methods", value_name="Values")

        # # Create the violin plot
        # plt.figure(figsize=(12, 6))
        # sns.violinplot(data=melted_df, x="num_genes", y="Values", hue="Methods", palette="muted", split=True)
        # plt.title("Violin Plot of Values Against num_genes")
        # plt.show()
        # combined_df.set_index("RefOGs", inplace=True)
        # combined_df.columns = methods
        # corr_analysis(combined_df, predog_col + " vs. " + eval_col, local_score_fig)

        # sns.violinplot(data=combined_df, x="day", y="total_bill", hue="smoker",
        #        split=True, inner="quart", fill=False,
        #        palette={"Yes": "g", "No": ".35"})
print(dist_dict)
sorted_columns = sorted([*combined_df.columns], key=lambda col: dist_dict[col])

combined_df = combined_df[sorted_columns]
melted_df = combined_df.melt(var_name="Methods", value_name="Values")
print(melted_df)
# Create the violin plot for all distributions
plt.figure(figsize=(16, 8))
sns.violinplot(
    data=melted_df,
    x="Methods",    # All columns, including num_genes, will be treated as methods
    y="Values",     # Plot their values as distributions
    palette=sns.color_palette("muted", len(melted_df["Methods"].unique())),  # Custom palette
    dodge=False,
    hue="Methods",  # Use Methods for color differentiation
    legend=False,   # Disable legend if not needed
    split=True, 
    inner="quart",
    gap=0,
)

plt.title("Distributions of number of genes in RefOGs against different Methods", fontsize=18)
plt.xticks(rotation=45, fontsize=14)  # Rotate x-tick labels and increase font size
plt.yticks(fontsize=14)               # Increase y-tick font size
plt.ylabel("")     # Increase y-label font size
plt.xlabel("")    # Increase x-label font size
plt.tight_layout()

# Show the plot
plt.show()