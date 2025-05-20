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

overlap_predog_dir = r"./OrthoBench/scores_preprocessed_predicted_orthogroups/overlap_predogs"

metric = "Genes" #"Fission Genes" # "Missing Genes" "Missing Species"
predog_df_dict = {}
methods = []
left_method = "OF3_LV" # "OF3_Align_ST" # "OF3_TEST" "OF3_LV"
right_method = "OF3_N0_DM" # "OF3_DB" # "OF3_TEST" # "OF3_Align_ST" "OF3_N0_DM"
for file in pathlib.Path(overlap_predog_dir).iterdir():
    method = file.name.rsplit(".", 1)[0]
    if file.name.startswith("OF"):
    # if method in ["OF3_TEST", "OF3_Align_ST"]:
        predog_df_dict[method] = {}
        with open(file) as reader:
            for line in reader:
                predog_key, predog = line.strip().split(": ")
                predog_set = set(predog.split(", "))
                predog_df_dict[method][predog_key] = predog_set

        # df.rename(columns={metric: method}, inplace=True)
        # df_list.append(df)
        
        if method != left_method:
            methods.append(method)


# for i in range(len(methods)):
#     for j in range(i+1, len(methods)):
#         print(methods[i], methods[j])

of3_test = predog_df_dict[left_method]
min_left_diff_dict = {}
min_right_diff_dict = {}
predog_left_key_dict = {}
predog_right_key_dict = {}
predog_left_diff_dict = {}
predog_right_diff_dict = {}
predog_left_diff_count_dict = {}
predog_right_diff_count_dict = {}
of3_test_count_dict = {}

for test_predog_key, test_predog in of3_test.items():
    of3_test_count_dict[test_predog_key] = len(test_predog)
    min_left_diff_dict[test_predog_key] = {}
    predog_left_diff_dict[test_predog_key] = {}
    predog_left_diff_count_dict[test_predog_key] = {}
    predog_left_key_dict[test_predog_key] = {}
    

    # min_right_diff_dict[test_predog_key] = {}
    # predog_right_diff_dict[test_predog_key] = {}
    # predog_right_diff_count_dict[test_predog_key] = {}
    # predog_right_key_dict[test_predog_key] = {}

    for method in methods:
        min_left_diff_dict[test_predog_key][method] = {}
        predog_left_diff_dict[test_predog_key][method] = {}
        predog_left_diff_count_dict[test_predog_key][method] = {}
        predog_left_key_dict[test_predog_key][method] = {}

        # min_right_diff_dict[test_predog_key][method] = {}
        # predog_right_diff_dict[test_predog_key][method] = {}
        # predog_right_diff_count_dict[test_predog_key][method] = {}
        # predog_right_key_dict[test_predog_key][method] = {}
        min_diff = len(test_predog)
        
        for predog_key, predog in predog_df_dict[method].items():
            left_n_diff = len(test_predog - predog)
            right_n_diff = len(predog - test_predog)
            if left_n_diff < min_diff:
                min_left_diff_dict[test_predog_key][method] = left_n_diff
                predog_left_diff_dict[test_predog_key][method] = test_predog - predog
                predog_left_diff_count_dict[test_predog_key][method] = len(predog)
                predog_left_key_dict[test_predog_key][method] = predog_key

            # if right_n_diff < min_diff:
            #     min_right_diff_dict[test_predog_key][method] =  right_n_diff
            #     predog_right_diff_dict[test_predog_key][method] = predog - test_predog
            #     predog_right_diff_count_dict[test_predog_key][method] = len(predog)
            #     predog_right_key_dict[test_predog_key][method] = predog_key
    
    # for method in methods:
    
    if min_left_diff_dict[test_predog_key][right_method] != 0:
        print(left_method, test_predog_key, of3_test_count_dict[test_predog_key])
        print(right_method, predog_left_key_dict[test_predog_key][method], predog_left_diff_count_dict[test_predog_key][right_method])
        print(min_left_diff_dict[test_predog_key][right_method], predog_left_diff_dict[test_predog_key][right_method])
        print()
        

    # if min_right_diff_dict[test_predog_key]["OF3_Align_ST"] != 0 :
    #     print(min_right_diff_dict[test_predog_key]["OF3_Align_ST"])
    #     print(predog_right_key_dict[test_predog_key]["OF3_Align_ST"])
    #     # print(predog_right_diff_dict[test_predog_key]["OF3_Align_ST"])
    #     print(of3_test_count_dict[test_predog_key], predog_right_diff_count_dict[test_predog_key]["OF3_Align_ST"])


    