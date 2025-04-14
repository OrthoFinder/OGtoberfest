import matplotlib.pyplot as plt 
import numpy as np 
import pathlib
import functools
import pandas as pd 
import matplotlib.lines as mlines
import matplotlib.ticker as mticker
import seaborn as sns


def local_score_display(combined_df, x_col, title, d_output=None):
    num_cols = len(combined_df.columns[1])
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

    for ax, col in zip(axes, combined_df.columns[1]):
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




missing_genes_dir = r"./OrthoBench/scores_preprocessed_predicted_orthogroups/missing_genes"
fission_dir = r"./OrthoBench/scores_preprocessed_predicted_orthogroups/fission"

metric = "Missing Species" #"Fission Genes" # "Missing Genes" "Missing Species"
df_list = []
methods = []
for file in pathlib.Path(missing_genes_dir).iterdir():
    
    method = file.name.rsplit(".", 1)[0]
    if file.name.startswith("OF"):
    # if method in ["OF3_TEST", "OF3_Align_ST"]:
        df = pd.read_csv(file, sep="\t", usecols=[0, 1])
        df.rename(columns={metric: method}, inplace=True)
        df_list.append(df)
        methods.append(method)

combined_df = functools.reduce(
    lambda left, right: pd.merge(left, right, on='RefOGs'), df_list)


# df_all = combined_df.copy()
# df_all.set_index("RefOGs")
df_all = combined_df.dropna(how="all", subset=methods)
df_all.reset_index(drop=True, inplace=True)
df_all.to_csv(r"./missing_genes_of.txt", sep="\t", index=False)
# combined_df.to_csv(r"./missing_genes_of.txt", sep="\t", index=False)
# rows_as_dicts = combined_df.to_dict(orient='records')

missing_genes_dict = {}
missing_genes_dict_count = {}
missing_genes_dict_symdiff = {}
missing_genes_dict_count_symdiff = {}
for row in combined_df.itertuples():
    row_dict = row._asdict()
    mg_set = set()
    # for method in ["OF3_TEST", "OF3_Align_ST"]:
    for method in methods:
        if row_dict[method] is np.nan:
            row_dict[method] = set()
        else:
            row_dict[method] = set(row_dict[method].split(", "))
        # if method in ["OF3_TEST", "OF3_Align_ST"] and len(row_dict[method]) != 0:
        #     print(row_dict['RefOGs'], method, row_dict[method])
        mg_set.update(row_dict[method])
    # if len(mg_set) != 0:
    #     print(row_dict['RefOGs'], mg_set.pop())
    missing_genes_dict[row_dict['RefOGs']] = mg_set
    missing_genes_dict_count[row_dict['RefOGs']] = len(mg_set)

    # missing_genes_dict_symdiff[row_dict['RefOGs']] = \
    #     row_dict["OF3_TEST"].symmetric_difference(row_dict["OF3_Align_ST"])
    # missing_genes_dict_count_symdiff[row_dict['RefOGs']] = \
    #     len(missing_genes_dict_symdiff[row_dict['RefOGs']])


# for refog_key, mg_set in missing_genes_dict.items():
#     mg_count = missing_genes_dict_count[refog_key]
#     if mg_count != 0:
#         print(refog_key, mg_count)
#         print(mg_set)

# mg_count_df = pd.DataFrame.from_dict(missing_genes_dict_count_symdiff, orient='index')
mg_df = pd.DataFrame.from_dict(missing_genes_dict, orient='index')
mg_df.dropna(inplace=True)
mg_df.reset_index(inplace=True)
mg_df.columns = ["RefOGs", f"{metric} Count"]
print(mg_df)

mg_count_df = pd.DataFrame.from_dict(missing_genes_dict_count, orient='index')
mg_count_df.reset_index(inplace=True)
mg_count_df.reset_index(drop=True, inplace=True)
mg_count_df.columns = ["RefOGs", f"{metric} Count"]

mg_count_df.sort_values(by=f"{metric} Count", ascending=False, inplace=True)
mg_count_df = mg_count_df[mg_count_df[f"{metric} Count"] != 0]
print(mg_count_df)

fig, ax = plt.subplots(layout="constrained")

x_positions = range(mg_count_df.shape[0])
N = len(x_positions)
palette = sns.color_palette("deep", n_colors=N)
sns.barplot(mg_count_df, x="RefOGs", 
            y=f"{metric} Count", 
            hue="RefOGs", 
            palette="deep", 
            # palette=color_palette,
            # palette="Paired",
            ax=ax,
            width=1.0
            )
ax.axhline(0, color="k", clip_on=False)
ax.set_xticks(range(len(mg_count_df["RefOGs"])))
ax.set_xticklabels(mg_count_df["RefOGs"], rotation=90, ha="center", size=15)
ax.set_xlabel("")
ax.set_ylabel(f"{metric} Count", size=15)
plt.show()
