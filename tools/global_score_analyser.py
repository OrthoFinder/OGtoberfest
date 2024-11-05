import pandas as pd 
import numpy as np 
import seaborn as sns 
import matplotlib.pyplot as plt

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


input_path = r"./OrthoBench/scores_preprocessed_predicted_orthogroups/OrthoBench_global_scores.tsv"

df = pd.read_csv(input_path, sep="\t")

plot_heatmap(df.set_index("Methods"))
