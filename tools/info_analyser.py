import os
import pathlib
import pandas as pd 
import numpy as np 
import seaborn as sns 
import matplotlib.colors as mcolors
import matplotlib.lines as mlines
import matplotlib.pyplot as plt
import matplotlib.colors as mcolors
from scipy.spatial.distance import squareform
from scipy.cluster.hierarchy import linkage
sns.set_theme(style="white", context="talk")


input_path = r"./OrthoBench/scores_preprocessed_predicted_orthogroups/local_score_distance"

for file in pathlib.Path(input_path).iterdir():
    # print(file)
    if file.name == "avg_VI.tsv":
        df = pd.read_csv(file, sep="\t")
        df.set_index(" ", inplace=True)

        df = (df + df.T) / 2
        np.fill_diagonal(df.values, 0)

        condensed_matrix = squareform(df)
        linkage_matrix = linkage(condensed_matrix, method="ward")


        ax = sns.clustermap(
            df,
            annot=True, 
            fmt=".2f", 
            annot_kws={"size": 15}, 
            cmap="cubehelix", 
            row_cluster=True,  # Enable clustering
            col_cluster=True,  # Enable clustering
            row_linkage=linkage_matrix,  # Use precomputed linkage for rows
            col_linkage=linkage_matrix,   # Use precomputed linkage for columns
            cbar_pos=None
        )

        ax.ax_heatmap.tick_params(axis='both', which='both', labelsize=15)
        ax.ax_heatmap.set_title(
            "Variation of Information Distance", 
            # file.name.split(".", 1)[0],  # Use your desired title string here
            fontsize=16,
            pad=20  # Adjust the padding between the title and the plot
        )
        plt.tight_layout()
        plt.show()
