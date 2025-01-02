import os
import pandas as pd 
import numpy as np 
import seaborn as sns 
import matplotlib.colors as mcolors
import matplotlib.lines as mlines
import matplotlib.pyplot as plt
import matplotlib.colors as mcolors
from met_brewer import met_brew
colors = met_brew(name="VanGogh2", n=13)#, brew_type="continuous") Signac
sns.set_theme(style="white", context="talk")
colors = colors[::-1]
# colors = set([
#     "goldenrod",
#     "olive",
#     "indigo",
#     "darkgreen",
#     "lightblue",
#     "blue",
#     # "orange",
#     "yellowgreen",
#     "plum",
#     "purple",
#     "silver",
#     "tomato",
#     "turquoise",
#     "pink",
#     "brown",
# ])
# '#e6194b', '#3cb44b', '#ffe119', '#4363d8', '#f58231', 
# '#911eb4', '#46f0f0', '#f032e6', '#bcf60c', '#fabebe', 
# '#008080', '#e6beff', '#9a6324', '#fffac8', '#800000', 
# '#aaffc3', '#808000', '#ffd8b1', '#000075', '#808080', '#ffffff', '#000000'



# colors = set([
#     # "#800000", # Maroon 1
#     "#9A6324", # Brown 2
#     "#808000", # Olive 3
#     "#469990", # Teal 4
#     "#000075", # Navy 5
#     # "#e6194B", # Red 6
#     # "#3cb44b", # Green 7
#     "#ffe119", # Yellow 8
#     "#4363d8", # Blue 9
#     "#f58231", # Orange 10
#     "#911eb4", # Purple 11
#     # "#42d4f4", # Cyan 12
#     "#f032e6", # Magenta 13
#     "#bfef45", # Line 14
#     # "#fabed4", # Pink 15
#     "#dcbeff", # Lavender 16
#     # "#fffac8", # Beige 17
#     "#aaffc3", # Mint 18
#     "#ffd8b1", # Apricot 19
#     # "#a9a9a9", # Grey 20
# ])

# color_palette = {
#     "OF3_Align_ST", "#C9A227",
#     "OF3_Linear", "#EDC531",
#     "OF3_DB", "#FFE169",
#     "OF2_Align_ST", "#A8A9AD",
#     "OF2_DB", "#D8D8D8",
#     "Broccoli", "#4C8C25",
#     "SP2_sens", "#7B3F00",
#     "SP2_def", "#954535",
#     "SP2_fast", "#B87333",
#     "ProteinOrtho", "#415D85",
#     "Hieranoid", "#FFAD98",
#     "FastOMA", "#FFD89C",
#     "OrthoMCL", "#9FCCAD",
# }


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
    # color_palette = {
    #     "OF3_Align_ST", "#C9A227",
    #     "OF3_Linear", "#EDC531",
    #     "OF3_DB", "#FFE169",
    #     "OF2_Align_ST", "#A8A9AD",
    #     "OF2_DB", "#D8D8D8",
    #     "Broccoli", "#4C8C25",
    #     "SP_sens", "#7B3F00",
    #     "SP_def", "#954535",
    #     "SP_fast", "#B87333",
    #     "ProteinOrtho", "#415D85",
    #     "Hieranoid", "#FFAD98",
    #     "FastOMA", "#FFD89C",
    #     "OrthoMCL", "#9FCCAD",
    # }

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

    # color_palette = {
    #     "OF3_Align_ST": "#C9A227",
    #     "OF3_Linear": "#EDC531",
    #     "OF3_DB": "#FFE169",
    #     "OF2_Align_ST": "#A8A9AD",
    #     "OF2_DB": "#D8D8D8",
    #     "Broccoli": "#4C8C25",
    #     "SP_sens": "#7B3F00",
    #     "SP_def": "#954535",
    #     "SP_fast": "#B87333",
    #     "ProteinOrtho": "#415D85",
    #     "Hieranoid": "#FFAD98",
    #     "FastOMA": "#FFD89C",
    #     "OrthoMCL": "#9FCCAD",
    # }

    num_cols = len(df.columns[1:])
    ncols = 4
    nrows = int(np.ceil(num_cols / ncols))

    fig, axes = plt.subplots(
        nrows=nrows, 
        ncols=ncols, 
        # figsize=(4 * nrows, 3 * ncols), 
        sharex=False, 
        sharey=False,
        layout="constrained"
    )

    axes = axes.flatten() if ncols > 1 else [axes]

    for ax, col in zip(axes, df.columns[1:]):
        # Sort the DataFrame based on the column
        if col in ["RefOG Fissions (%)", "Missing Genes (%)", "Missing Species (%)", "Entropy", "Effective Size", "Runtime", "Rank Score"]:
            df.sort_values(by=[col], inplace=True)
        elif col in ["Recall", "Precision"]:
            df.sort_values(by=[col], inplace=True, ascending=False)
        else:
            df.sort_values(by=["Rank Score"], inplace=True)

    #     sns.barplot(
    #         data=df, 
    #         y="Methods",  # Horizontal barplot requires swapping x and y
    #         x=col, 
    #         hue="Methods", 
    #         palette=color_palette,
    #         ax=ax,
    #         width=1.0,  # Adjust width for better visualization
    #         saturation=10
    #     )
        
    #     # Add a vertical line at x=0 for reference
    #     ax.axvline(0, color="k", clip_on=False)

    #     # Customize axis labels and ticks
    #     ax.set_yticks(range(len(df["Methods"])))
    #     ax.set_yticklabels(df["Methods"], ha="right", size=13)
    #     ax.set_ylabel("")
    #     ax.set_xlabel(col, size=13)

    # # Remove unused subplots
    # for i in range(len(df.columns[1:]), len(axes)):
    #     fig.delaxes(axes[i])
    # plt.show()

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
        if col in ["RefOG Fusions (%)", "RefOG Fissions (%)", "Missing Genes (%)", "Missing Species (%)", "Entropy", "Effective Size", "Runtime", "Rank Score"]:
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
# input_path = r"./Sim1k/scores_preprocessed_predicted_orthogroups/Sim1k_global_scores.tsv"

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
        # "Runtime",
        # "Effective Size",
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
        # True, 
        # True
    ],
)

print(df)
# plot_heatmap(df.set_index("Methods"))
plot_barplot(df)
