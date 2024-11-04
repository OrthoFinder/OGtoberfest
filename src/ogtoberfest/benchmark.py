#!/usr/bin/env python3
"""
This script calculates the benchmarks for an input set of orthogroups

Instructions:
1. Predict the complete set of orthogroups for the genes in the "Input/" directory
2. Write Orthogroups to a file, one orthogroup per line (header lines and commented
   out lines starting with '#' are allowed). 
3. Download this script and the accompanying RefOGs directory
3. Call the script with the the orthogroup filename as the only argument

By default the script with use regular expressions to extract the genes from the 
additional text on each line. 

You can also specify the option '-b' to use the more basic file reader which requires 
the following format:
- First line is a header and is ignored. 
- One orthogroup per line
- Genes can be separated by commas, spaces or tabs. 
- Lines starting with '#' are comments and are ignored
- Optionally, each line can start with the name of the orthogroup followed by a colon.
   E.g. "OG0000001: Gene1, Gene2"
"""
import os
import re
import sys
import glob
import argparse
import pathlib
import numpy as np
import pandas as pd
import math
import traceback
from typing import Optional, List, Tuple, Dict, Set, Iterator
from itertools import product
import matplotlib.pyplot as plt
import seaborn as sns
import scorefuncs as sf 


COMPILE = re.compile(r"[_\s]+")

delims = " |,|\t"
expected_gene_names_base = {"ENSP000": "Homo sapiens", 
    "ENSRNO":"Rattus norvegicus", 
    "ENSCAF":"Canis familiaris", 
    "ENSMUS":"Mus musculus", 
    "ENSMOD":"Monodelphis domestica", 
    "ENSGAL":"Gallus gallus", 
    "ENSCIN":"Ciona intestinalis",
    "ENSTNI":"Tetraodon nigroviridis",
    "ENSPTR":"Pan troglodytes",
    "ENSDAR":"Danio rerio",
    "FBpp":"Drosophila melanogaster",
    "WBGene":"C elegans"}
n_genes_total = 251378


def plot_clustermap(symmetric_df,  data_dirname=None, output_filename=None,
                 gc="orthogroup", analysis_type="intersection",
                 save_to="orthogroups_compare", title="Variation of Information"):

    # plt.figure(figsize=(12, 12))
    # ax = sns.heatmap(symmetric_df, annot=True, fmt=".0f", annot_kws={"size": 15}, cmap="viridis")
    ax = sns.clustermap(symmetric_df, annot=True, fmt=".2f", annot_kws={"size": 15}, cmap="viridis")
    ax.tick_params(axis='both', which='both', labelsize=15)

    # ax.set_yticklabels(ax.get_yticklabels(), rotation=0, ha='right')
    # cbar = ax.collections[0].colorbar
    
    # cbar.formatter.set_powerlimits((0, 0))
    # cbar.ax.tick_params(labelsize=15)
    # cbar.update_ticks()

    # if not title:
    #     plt.title(f"{gc}: {analysis_type} %", fontsize = 12)
    # else:
    #     plt.title(title, fontsize = 12)
    plt.tight_layout()
    plt.show()
    # output_dir = pathlib.Path.cwd() / "_".join((data_dirname, save_to))

    # if not os.path.exists(output_dir):
    #     os.makedirs(output_dir, exist_ok=True)
    
    # output_path = output_dir / f"{output_filename}.png"
    # plt.savefig(output_path)
    # plt.close()


def exit():
    sys.exit()


def read_hierarchical_orthogroup(infile):
    ogs = []
    iSkip = 3 
    for line in infile:
        genes = [g for s in line.rstrip().split("\t")[iSkip:] for g in s.split(", ") if g != ""]
        ogs.append(set(genes))
    return ogs


def read_orthogroups_smart(file_path, 
                           skip_nrows: Optional[int] = None,
                           skip_ncols: Optional[int] = None, 
                           use_reg: bool = True, 
                           outgroups: Optional[List[str]] = None):
    ogs = {}
        
    gene_pat = re.compile("WBGene00\d+\.1|ENSCAFP\d+|ENSCINP\d+|ENSDARP\d+|FBpp0\d+|ENSGALP\d+|ENSP000\d+|ENSMODP\d+|ENSMUSP\d+|ENSPTRP\d+|ENSRNOP\d+|ENSTNIP\d+")
    base_predog_format = "PredOG%07d"
    with open(file_path, 'r') as infile:
        for i, l in enumerate(infile):
            # if l.startswith("#"):
            #     continue
            if use_reg:
                genes = re.findall(gene_pat, l)
                # print(genes)
                # sys.exit()
                if len(genes) > 0:
                    if ":" in l:
                        og_label = l.split(":")[0]
                    else:
                        og_label = l.split(None, 1)[0]
                    ogs[og_label] = set(genes)
                    # ogs.append(set(genes))
            else:
                if skip_nrows is not None:
                    if skip_nrows == (i + 1):
                        print("Skipping the following row!")
                        print(l)
                        continue

                else:
                    print("Please provide the number of colnums to skip")
                    Fail()
                if ":" in l:
                    og_label = l.split(":")[0]
                    genes = re.split(delims, l.split(":")[1].strip())
                    genes = [g.strip().split("|")[-1] for g in genes if g.strip() != "" ]
                else:
                    # item = l.strip().split(None, skip_ncols)
                    # og_label, t = item[0], item[-1]
                    if skip_ncols is not None:
                        if skip_ncols == 0:
                            t = l.strip()
                        else:
                            t = l.strip().split(None, skip_ncols)[-1]
                    else:
                        print("Please provide the number of colnums to skip")
                        Fail()

                    og_label = base_predog_format % i
                    
                    genes = re.split(delims, t)
                    genes = [g.strip().split("|")[-1] for g in genes if g.strip() != "" ]
                    
                    # genes = [g for g in genes if g != ""]
                    # og_label = l.split(None, 1)[0]
                    # genes = l.split(None, 3)[-1].strip().split()
                # print(og_label, set(genes))
                ogs[og_label] = set(genes)
    return ogs

def read_orthogroups(fn, exp_genes, n_col_skip=1):
    """
    Read the orthogroups from a file formatted as specified above
    """
    ogs = []
    q_past_header = False
    with open(fn, 'r') as infile:
        for l in infile:
            if l.startswith("#"):
                continue
            t = l.rstrip().split(None, n_col_skip)[-1]
            genes = re.split(delims, t)
            # remove white spaces
            genes = [g.strip() for g in genes]
            genes = set([g for g in genes if g != ""])
            if q_past_header or (exp_genes & genes):
                q_past_header = True
                if len(genes) > 1: 
                    ogs.append(genes)
    return ogs


def get_n_col_skip(fn):
    with open(fn, 'r') as infile:
        header = next(infile)
        if header.startswith("HOG\tOG\tGene Tree Parent Clade"):
            return 3
    return 1


def check_orthogroups(ogs, exp_genes):
    all_pred_genes = set([g for og in ogs for g in ogs[og]])
    # all_pred_genes = set([g for og in ogs for g in og])
    x = all_pred_genes.difference(exp_genes)
    if len(x) != 0:
        print("ERROR: found extra genes in input file, check its formatting is correct and there are no incorrect genes")
        print("Examples:")
        for g in list(x)[:10]:
            print(g)
    x = exp_genes.difference(all_pred_genes)
    if len(x) != 0:
        print("Examples of genes not in file:")
        for g in list(x)[:3]:
            print(g)

    n_genes = sum([len(ogs[og]) for og in ogs])
    # n_genes = sum([len(og) for og in ogs])
    if n_genes < 0.5 * n_genes_total:
        print("ERROR: Too many missing genes in predicted orthogroups.")
        print("Orthogroups should contain at least 50% of all genes but")
        print("orthogroups file only contained %d genes" % n_genes)
        exit()

    all_genes = [g for og in ogs for g in ogs[og]]
    # all_genes = [g for og in ogs for g in og]
    n_genes_no_dups = len(set(all_genes))
    if n_genes_no_dups != n_genes:
        print("ERROR: Some genes appear in multiple orthogroups, benchmark are meaningless with such data.")
        print("with such data")
        print((n_genes_no_dups, n_genes))
        from collections import Counter
        c = Counter(all_genes)
        for g, n in c.most_common(10):
            print("%d: %s" % (n,g))
        raise Exception()
    # checked genes from each of the expected species are present
    for g_pat, sp in expected_gene_names_base.items():
        if not any(g.startswith(g_pat) for g in all_genes):
            print("ERROR: No genes found from %s" % sp)
    p = 100.*n_genes/float(n_genes_total)
    print("%d genes found in predicted orthogroups, this is %0.1f%% of all genes" % (n_genes, p))


# def calculate_benchmarks_pairwise(ref_ogs: Dict[str, str], 
#                                   pred_ogs: Dict[str, str], 
#                                   uncert_genes: Optional[Dict] = None, 
#                                   q_even: bool = True, 
#                                   q_remove_uncertain: bool = True,
#                                   ):
#     referenceOGs = ref_ogs
#     predictedOGs = pred_ogs
#     totalFP = 0.
#     totalFN = 0.
#     totalTP = 0.
#     totalGroundTruth = 0.
#     n_exact = 0
#     n_splits = []
#     # so as not to count uncertain genes either way remove them from the 
#     # expected and remove them from any predicted OG (as though they never existed!)
#     missing_genes_dict = {}
#     intersection_dict = {}
#     union_dict = {}
#     FN_genes_num_dict = {}
#     FP_genes_num_dict = {}
#     FP_genes_dict = {}
#     FN_genes_dict = {}
#     TP_genes_dict = {}
#     refogs_num_dict = {}
#     purity_og_total = []
#     purity_og_max = []
#     predOg_size = {}
#     totalGenes = 0.
#     overlap_predOGs = {}
#     refOGs_set = set()
#     predOGs_overlap_set = set()
#     predOGs_all_set = set()
    
#     total_nPredOGs = len(predictedOGs)
#     unique_predOGs = {}

#     entropy_pred_all = 0.0
#     total_genes_nPredOGs = np.sum([len(predictedOGs[predog_key]) for predog_key in predictedOGs])
#     for predog_key in predictedOGs:
#         predOGs_all_set = predOGs_all_set.union(predictedOGs[predog_key])
#         p_pred_all = len(predictedOGs[predog_key]) / total_genes_nPredOGs
#         entropy_pred_all += np.sum(-p_pred_all * np.log2(p_pred_all))
    
#     for refog_key in referenceOGs:
#         refOg = referenceOGs[refog_key]
#         if uncert_genes is not None:
#             uncert = uncert_genes[refog_key]

#         thisFP = 0.
#         thisFN = 0.
#         thisTP = 0.
#         this_split = 0
#         if q_remove_uncertain and uncert_genes is not None:
#             refOg = refOg.difference(uncert)

#         refOGs_set = refOGs_set.union(refOg)
#         nRefOG = len(refOg)
#         totalGenes += nRefOG
#         not_present = set(refOg)
#         intersection = 0

#         refogs_num_dict[refog_key] = nRefOG
#         intersection_dict[refog_key] = {}
#         FN_genes_num_dict[refog_key] = {}
#         FP_genes_num_dict[refog_key] = {}
#         FP_genes_dict[refog_key] = {}
#         FN_genes_dict[refog_key] = {}
#         TP_genes_dict[refog_key] = {}
#         union_dict[refog_key] = {}
#         predOg_size[refog_key] = {}
#         overlap_predOGs[refog_key] = []
#         predOg_genes_set = set()
#         overlap_list = []
#         predognum = 0

#         for predog_key in predictedOGs:
#             predOg = predictedOGs[predog_key]
#             overlap = len(refOg.intersection(predOg))

#             intersection += overlap

#             if overlap > 0:
#                 if q_remove_uncertain and uncert_genes is not None:
#                     predOg = predOg.difference(uncert)   # I.e. only discount genes that are uncertain w.r.t. this RefOG
#                 overlap = len(refOg.intersection(predOg))

#             if overlap > 0:
#                 this_split += 1
#                 not_present = not_present.difference(predOg)
                
#                 thisTP += overlap * (overlap - 1)/2    # n-Ch-2
#                 thisFP += overlap * (len(predOg) - overlap)
#                 thisFN += (nRefOG - overlap) * overlap
                
#                 if predog_key not in intersection_dict[refog_key]:
#                     if predog_key not in unique_predOGs:
#                         unique_predOGs[predog_key] = len(predOg)
#                     predOGs_overlap_set = predOGs_overlap_set.union(predOg)

#                     overlap_predOGs[refog_key].append(predog_key)
#                     predOg_genes_set = predOg_genes_set.union(predOg)
#                     union_dict[refog_key][predog_key] = len(refOg | predOg)
#                     overlap_list.append(overlap)
#                     intersection_dict[refog_key][predog_key] = overlap
#                     TP_genes_dict[refog_key][predog_key]  = refOg & predOg

#                     FN_genes_num_dict[refog_key][predog_key] = len(refOg - predOg)
#                     FN_genes_dict[refog_key][predog_key] = refOg - predOg
                    
#                     FP_genes_num_dict[refog_key][predog_key] = len(predOg - refOg)
#                     FP_genes_dict[refog_key][predog_key] = predOg - refOg
#                     predOg_size[refog_key][predog_key] = len(predOg)

#             else:
#                 overlap_list.append(0)
    
#         missing_genes_dict[refog_key] = refOg - predOg_genes_set
#         purity_og_max.append(overlap_list)    
#         # print(intersection_dict[refognum])

#         # Are FNs more from splintered OGs or missing genes?
#         # print("%f\t%f" % (thisFN/2./(nRefOG-1), len(not_present)*(nRefOG-1)/2./(nRefOG-1))) 
#         # finally, count all the FN pairs from those not in any predicted OG
#         thisFN += len(not_present)*(nRefOG-1)
#         # don't count 'orphan genes' as splits, it's more informative only to count 
#         # clusters that this orthogroup has been split into. Recall already counts
#         #  the missing genes, this should give a distinct measure.
#         # this_split += len(not_present)     
#         # All FN have been counted twice
#         assert(thisFN % 2 == 0)
#         n_splits.append(this_split)
#         # print(this_split)      # Orthogroup fragments
#         # print(len(not_present))  # Unclustered genes 
#         thisFN /= 2 
#         # sanity check
#         nPairs1 = thisTP + thisFN
#         nPairs2 = nRefOG * (nRefOG - 1) / 2
#         if nPairs1 != nPairs2:
#             print("ERROR: %d != %d" % (nPairs1,nPairs2))    

#         totalGroundTruth += nPairs1
#         if thisFN == 0 and thisFP == 0:
#             n_exact += 1
#         if q_even:
#             N = float(len(refOg)-1)
#             totalFN += thisFN/N
#             totalFP += thisFP/N
#             totalTP += thisTP/N
#         else:
#             totalFN += thisFN
#             totalFP += thisFP
#             totalTP += thisTP
#         # print("%d\t%d\t%d" % (thisTP, thisFN, thisFP))  

#     purity_og_max_arr = np.amax(np.stack(purity_og_max).T, axis=1)

#     # print(np.sum(purity_og_max_arr[purity_og_max_arr>0]), totalGenes)
#     matched_col = [
#         "RefOG", 
#         "RefOG_Size", 
#         "max_overlap", 
#         "max_overlap_FN", #"max_overlap_FN (%)", 
#         "max_overlap_FP", #"max_overlap_FP (%)", 
#         "max_overlap_Precision (%)", "max_overlap_Recall (%)", "max_overlap_F-score (%)", 
#         "avg_Precision (%)", "avg_Recall (%)", "avg_F-score (%)", 
#         "Entropy", 
#         "Distance",
#         "nPredictedOGs", 
#         "TotalMissingGenes", 
#         "EffectiveSize", 
#         "PredOG:Size:Intersection:JarccardIndex:Precision:Recall", 
#         "Missing_Genes",
#         "Overlapping_Genes"
#     ]                   
                      
#     matched_results = []
#     missing_genes = ""
#     totalEntropy_list = []
#     total_avg_presision = 0.0
#     total_avg_recall = 0.0
#     total_avg_fscore = 0.0
#     round_precision = 1
#     all_missing_genes = 0
#     totalEffectiveSize = 0
#     nRefOG_list = []
#     nPredOG_effective_size_list = []
#     total_overlap_predogs = 0
#     total_overlaps = 0
#     unique_total_overlap_predogs = set()
#     total_distance_list = []
#     total_RefOGs = len(refogs_num_dict)
#     unidentified_refOGs = []
#     p_ref_pred_dict = {}
#     p_ref_dict = {}
#     total_overlap_genes_nPredOGs = np.sum([*unique_predOGs.values()])
#     p_pred_dict = {predog_key: val / total_overlap_genes_nPredOGs for predog_key, val in unique_predOGs.items()}
#     entropy_ref = 0.0
#     for refog_key in refogs_num_dict:
#         nRefOG = refogs_num_dict[refog_key]
#         p_ref = nRefOG / totalGenes
#         entropy_ref += np.sum(-p_ref * np.log2(p_ref))
#         p_ref_dict[refog_key] = p_ref 
#         p_ref_pred_dict[refog_key] = {}

#         total_overlap_predogs += len(overlap_predOGs[refog_key])
#         unique_total_overlap_predogs = unique_total_overlap_predogs | set(overlap_predOGs[refog_key])
        
#         # overlaps_list = []
#         if len(intersection_dict[refog_key]) > 0:
#             nRefOG_list.append(nRefOG)

#             jaccard_index_list = [(key, val / union_dict[refog_key][key]) for key, val in intersection_dict[refog_key].items()]
#             sorted_jaccard_index_list = sorted(jaccard_index_list, reverse=True, key=lambda x: x[1])
#             item_list = [(key, intersection_dict[refog_key][key], jaccard_index) for key, jaccard_index in sorted_jaccard_index_list]
            
#             ogs, ogs_val, _ = zip(*item_list)
#             overlap_ogs = ",".join([og + ":" + str(predOg_size[refog_key][og]) + ":" \
#                                     + str(val) + ":" + str(np.round(100 * ja, round_precision)) + ":" \
#                                     + str(np.round(100 * val / predOg_size[refog_key][og], round_precision)) + ":" \
#                                     + str(np.round(100 * val / nRefOG, round_precision)) \
#                                           for og, val, ja in item_list])
#             num_predogs = len(ogs)
#             ogs_val = np.array(ogs_val)
#             # ogs_val_per = np.round(100 * ogs_val / nRefOG, round_precision)

#             max_overlap = ogs_val[0]
#             # max_overlap_per = ogs_val_per[0]

#             max_overlap_og = ogs[0]
#             overlaps = np.sum([*intersection_dict[refog_key].values()])
#             missing_genes = ""
#             if overlaps < nRefOG:
#                 missing_genes = ",".join(missing_genes_dict[refog_key])

#             total_missing_genes = len(missing_genes_dict[refog_key])
#             all_missing_genes += total_missing_genes
            
#             # refEntropy = 0
#             # for val in ogs_val:
#             #     refEntropy += -(val/nRefOG) * np.log2(val/nRefOG)
            
#             # totalEntropy += refEntropy * nRefOG / totalGenes
            
#             avg_jaccard_index = 0.0
#             avg_dismilarity = 0.0
#             avg_precision = 0.0
#             avg_recall = 0.0
#             avg_fscore = 0.0
#             Overlapping_Genes = ""
#             effective_size = 0
#             refEntropy = 0.0
#             if total_missing_genes != 0:
#                 refEntropy += -(total_missing_genes/nRefOG) * np.log2(total_missing_genes/nRefOG)
            
#             total_intersection = np.sum([*intersection_dict[refog_key].values()])

#             disimilarity_list = []
            
#             for predog_key in intersection_dict[refog_key]:
#                 overlap = intersection_dict[refog_key][predog_key]
#                 if overlap != 0:
#                     refEntropy += -(overlap / nRefOG) * np.log2( overlap / nRefOG)
                
#                 total_overlaps += overlap
                
#                 p_ref_pred_dict[refog_key][predog_key] = overlap / totalGenes

#                 jaccard_index = overlap / union_dict[refog_key][predog_key]
#                 avg_jaccard_index +=  jaccard_index * overlap / total_intersection
#                 avg_dismilarity += (1 - jaccard_index) * overlap / total_intersection
#                 disimilarity_list.append(1 - jaccard_index)
#                 # overlaps_list.append(overlap)

#                 nPredOG = predOg_size[refog_key][predog_key]
#                 predOg_precision = overlap / nPredOG
#                 predOg_recall = overlap / nRefOG
#                 predOg_fscore = 2 * predOg_precision * predOg_recall / (predOg_precision + predOg_recall)

#                 avg_precision += predOg_precision * overlap / total_intersection
#                 avg_recall += predOg_recall * overlap / total_intersection
#                 avg_fscore += predOg_fscore * overlap / total_intersection

#                 effective_size += int(overlap * overlap / nPredOG)
#                 # effective_size += int(nPredOG * jaccard_index)

#                 Overlapping_Genes += predog_key + ":" + ",".join(TP_genes_dict[refog_key][predog_key])
#                 Overlapping_Genes += "||"
            
#             Overlapping_Genes.rstrip("||")
#             nPredOG_effective_size_list.append(max(1, effective_size))
#             totalEffectiveSize += effective_size
#             totalEntropy_list.append(refEntropy * nRefOG / totalGenes)
#             # #==========================
#             total_avg_presision += avg_precision * nRefOG / totalGenes
#             total_avg_recall += avg_recall * nRefOG / totalGenes
#             total_avg_fscore += avg_fscore * nRefOG / totalGenes
#             # #================================
#             # total_avg_presision += avg_precision**2 #* nRefOG / totalGenes
#             # total_avg_recall += avg_recall**2 #* nRefOG / totalGenes
#             # total_avg_fscore += avg_fscore**2 #* nRefOG / totalGenes

#             disimilarity_arr = np.array(disimilarity_list)
            
#             distances = expon.ppf(disimilarity_arr, scale=nPredOG)
#             rms_total_distance = np.sqrt(np.mean(distances**2))
#             total_distance_list.append(rms_total_distance)

#             max_overlap_fn = FN_genes_num_dict[refog_key][max_overlap_og]
#             # max_overlap_fn_per = max_overlap_fn / nRefOG

#             max_overlap_fp = FP_genes_num_dict[refog_key][max_overlap_og]
#             # max_overlap_fp_per =  max_overlap_fp / nRefOG
            
#             max_overlap_pres = max_overlap / (max_overlap + max_overlap_fp)
#             max_overlap_recall = max_overlap / (max_overlap + max_overlap_fn) 
#             max_overlap_f =  2 * max_overlap_pres * max_overlap_recall / (max_overlap_pres + max_overlap_recall)
#         else:
#             p_ref_pred_dict[refog_key][predog_key] = 0.0
#             max_overlap = 0
#             max_overlap_fn = 0
#             max_overlap_fp = 0
#             max_overlap_pres = 0.0
#             max_overlap_recall = 0.0
#             avg_precision = 0.0
#             avg_fscore = 0.0
#             refEntropy = np.nan
#             nRefOG_list.append(0)
#             totalEntropy_list.append(refEntropy * nRefOG / totalGenes)
#             num_predogs = 0
#             total_missing_genes = nRefOG
#             effective_size = 0
#             overlap_ogs = ""
#             missing_genes = ""
#             Overlapping_Genes = ""
#             nPredOG_effective_size_list.append(effective_size)
#             avg_recall = 0.0
#             avg_precision = 0.0
#             avg_fscore = 0.0
#             max_overlap_f = 0.0
#             max_overlap_fn = 0.0
#             max_overlap_fp = 0.0
#             rms_total_distance = np.nan
#             total_distance_list.append(rms_total_distance)
#             totalEffectiveSize += 0
#             unidentified_refOGs.append(refog_key)
#             # overlaps_list.append(nRefOG)

#         matched_results.append((refog_key, nRefOG, 
#                                 max_overlap,
#                                 max_overlap_fn, 
#                                 # np.round(100 * max_overlap_fn_per, round_precision), 
#                                 max_overlap_fp,
#                                 # np.round(100 * max_overlap_fp_per, round_precision), 
#                                 np.round(100 * max_overlap_pres, round_precision), 
#                                 np.round(100 * max_overlap_recall, round_precision),
#                                 np.round(100 * max_overlap_f, round_precision),
#                                 np.round(100 * avg_precision, round_precision), 
#                                 np.round(100 * avg_recall, round_precision), 
#                                 np.round(100 * avg_fscore, round_precision), 
#                                 np.round(refEntropy, round_precision + 1), 
#                                 np.round(rms_total_distance, round_precision + 1), 
#                                 # refEntropy,
#                                 num_predogs, 
#                                 total_missing_genes, 
#                                 effective_size,
#                                 overlap_ogs, 
#                                 missing_genes,
#                                 Overlapping_Genes,
#                                 ))
#     nPredOG_effective_size_porportion = np.array(nPredOG_effective_size_list) / totalEffectiveSize
#     nRefOG_proportion_arr = np.array(nRefOG_list) / np.sum(nRefOG_list)

#     totalEntropy = np.nansum(totalEntropy_list)
#     total_avg_distance = np.nanmean(total_distance_list)

#     ## MI and VI  --------------------------------------------
#     MI = 0.0
#     epsilon = 1e-10
#     for refog_key in p_ref_pred_dict:
#         pi = p_ref_dict[refog_key]
#         for predog_key, pij in p_ref_pred_dict[refog_key].items():
#             pj = p_pred_dict[predog_key]
#             if pij > epsilon :
#                 MI += pij * np.log2(pij / (pi * pj))

#     overlap_entropy_ref = np.sum([-pi * np.log2(pi) for pi in p_ref_dict.values()])
#     overlap_entropy_pred = np.sum([-pj * np.log2(pj) for pj in p_pred_dict.values()])
#     VI = overlap_entropy_pred + overlap_entropy_ref - 2 * MI
#     ## -------------------------------------------------------
    
#     overlap_genes_nPredOGs = np.sum([*unique_predOGs.values()])
#     entropy_pred = np.sum([
#         -(val/overlap_genes_nPredOGs) * np.log2(val / overlap_genes_nPredOGs) for val in unique_predOGs.values()
#     ])

#     DKL = np.sum([
#         nRefOG_proportion_arr[i] * np.log2(nRefOG_proportion_arr[i] / nPredOG_effective_size_porportion[i])
#         if nPredOG_effective_size_porportion[i] != 0.0 else 0.0
#         for i in range(len(nRefOG_proportion_arr))
#     ])
#     print(DKL)
    
#     corss_entropy = np.sum([
#         -nRefOG_proportion_arr[i] * np.log2(nPredOG_effective_size_porportion[i])
#         if nPredOG_effective_size_porportion[i] != 0 else 0
#         for i in range(len(nRefOG_proportion_arr))
#     ])
    
#     matched_df = pd.DataFrame.from_records(matched_results, columns=matched_col)
#     # print(matched_df[["RefOG_intersection_PredOG", "min_FN", "FP"]])
#     matched_df.sort_values(by=["max_overlap_Recall (%)", 
#                                "max_overlap", 
#                                "max_overlap_FN", 
#                                "max_overlap_FP"], \
#                            inplace=True, ascending=[False, False, True, True])

#     exact_matched = matched_df.loc[(matched_df["max_overlap_Recall (%)"] == 100) & (matched_df["max_overlap_FP"] == 0), : ]
#     exact_matched_ogs = exact_matched.shape[0]
#     assert(n_exact == exact_matched_ogs)
    
#     # n_rows, n_cols = 3, 3
#     # fig, axs = plt.subplots(n_rows, n_cols, figsize=(15, n_rows * 5))  
#     # data_cols = [
#     #     "avg_Precision (%)", 
#     #     "avg_Recall (%)", 
#     #     "avg_F-score (%)", 
#     #     "Entropy", 
#     #     "Distance",
#     #     "RefOG_Size", 
#     #     "EffectiveSize"
#     # ]

#     # for ii, col in enumerate(data_cols):
#     #     ax = axs[ii // n_cols, ii % n_cols] 
#     #     sns.histplot(data=matched_df[col], alpha=0.5, bins=30, ax=ax, label=col, color='blue')
#     #     ax.set_xlabel(col)
#     #     ax.set_ylabel('Frequency')
#     # # ax.set_title(f'{fn}')
#     # ax.legend()
#     # ax.grid(True)
#     # plt.tight_layout()
#     # plt.show()

#     # # ----------------------------------------------
#     overlap_all = len(refOGs_set & predOGs_all_set)
#     assert overlap_all == total_overlaps, f"overlap_all {overlap_all} and otal_overlaps {total_overlaps} should be the same"
#     # FN_all = len(refOGs_set - predOGs_all_set)
#     # FP_all = len(predOGs_all_set - refOGs_set)
#     # FP_overlap = len(predOGs_overlap_set - refOGs_set)
#     Union_all = len(refOGs_set | predOGs_all_set)
#     Union_overlap = len(predOGs_overlap_set | refOGs_set)
#     recall_all = total_overlaps / totalGenes
#     precision_all = total_overlaps / total_genes_nPredOGs
#     recall_overlap = recall_all 
#     precision_overlap = total_overlaps / overlap_genes_nPredOGs
#     fscore_all = 2 * (recall_all * precision_all) / (recall_all + precision_all)
#     fscore_overlap = 2 * (recall_overlap * precision_overlap) / (recall_overlap + precision_overlap)
#     jaccard_index_all = overlap_all / Union_all
#     jaccard_index_overlap = overlap_all / Union_overlap
#     disimilarity_all = 1 - jaccard_index_all
#     disimilarity_overlap = 1 - jaccard_index_overlap
#     distance_all = expon.ppf(disimilarity_all, scale=total_RefOGs)
#     distance_overlap = expon.ppf(disimilarity_overlap, scale=total_RefOGs)
#     ## -----------------------------------------------   

    
#     TP, FP, FN = (totalTP, totalFP, totalFN)
#     # purity = np.sum(purity_og_max) / np.sum(purity_og_total)
#     purity = np.sum(purity_og_max_arr) / totalGenes
#     print("%d Correct gene pairs" % TP)
#     print("%d False Positives gene pairs" % FP)
#     print("%d False Negatives gene pairs\n" % FN)
#     pres = TP/(TP+FP)
#     recall = TP/(TP+FN)
#     f = 2*pres*recall/(pres+recall)

#     print("%0.1f%%  F-score" % (100.*f))
#     print("%0.1f%%  Precision" % (100.*pres))
#     print("%0.1f%%  Recall\n" % (100.*recall))
    
#     print("NEW SCORES")
#     print("%0.1f%%  F-score (weighted avg)" % (100. * total_avg_fscore))
#     print("%0.1f%%  Precision (weighted avg)" % (100. * total_avg_presision))
#     print("%0.1f%%  Recall (weighted avg)" % (100. * total_avg_recall))

#     print("%0.1f%%  F-score (overlap)" % (100. * fscore_overlap))
#     print("%0.1f%%  Precision (overlap)" % (100. * precision_overlap))
#     print("%0.1f%%  Recall (overlap)" % (100. * recall_overlap))

#     print("%0.1f%%  F-score (all)" % (100. * fscore_all))
#     print("%0.1f%%  Precision (all)" % (100. * precision_all))
#     # print("%0.1f%%  Recall (all)" % (100. * recall_all))

#     # print("%0.1f%%  F-score (avg)" % (100. * total_avg_fscore / len(intersection_dict)))
#     # print("%0.1f%%  Precision (avg)" % (100. * total_avg_presision / len(intersection_dict)))
#     # print("%0.1f%%  Recall (avg)" % (100. * total_avg_recall / len(intersection_dict)))

#     # print("%0.1f%%  F-score (rms)" % (100. * np.sqrt(total_avg_fscore /len(intersection_dict))))
#     # print("%0.1f%%  Precision (rms)" % (100. * np.sqrt(total_avg_presision / len(intersection_dict))))
#     # print("%0.1f%%  Recall (rms)" % (100. * np.sqrt(total_avg_recall / len(intersection_dict))))

#     print("%0.1f%%  Purity" % (100.*purity))
#     print("%0.2f   Entropy" % (totalEntropy))
#     print("%0.2f   Cross Entropy" % (corss_entropy))
#     print("%0.3f   KL-Divergence" % (DKL))
#     print("%0.3f   Overlap MI" % (MI))
#     print("%0.3f   Overlap VI" % (VI))

#     print("%0.3f   Overlap predOGs Entropy" % (entropy_pred))
#     print("%0.3f   Total predOGs Entropy" % (entropy_pred_all))

#     print("%0.2f   Total Average Distance" % (total_avg_distance))
#     print("%0.2f   Overlap Distance" % (distance_overlap))
#     print("%0.2f   Total Distance" % (distance_all))
#     print("%d      Total overlapped predOGs" % total_overlap_predogs)
#     print("%d / %d  Unidentified RefOgs / Total RefOGs" % (len(unidentified_refOGs), total_RefOGs))
#     print("%0.1f%%  Unidentified RefOgs" % (100. * len(unidentified_refOGs) / total_RefOGs))
#     print("%d       Total unique overlapped predOGs" % len(unique_total_overlap_predogs))
#     print("%d / %d  Total Missing Genes / Total Genes" % (all_missing_genes, totalGenes))
#     print("%0.1f%%  Missing Genes\n" % (100 * all_missing_genes / totalGenes))

#     print("%d Orthogroups exactly correct\n" % n_exact)
    
#     scores = [
#         100.*recall, 100.*pres, 100.*f, \
#         100.*recall_all, 100.*precision_all, 100.*fscore_all, \
#         100.*recall_all, 100.*precision_overlap, 100.*fscore_overlap, \
        
#         100.*total_avg_recall, \
#         100.*total_avg_presision , \
#         100.*total_avg_fscore, \

#         # 100.*total_avg_recall  / len(intersection_dict), \
#         # 100.*total_avg_presision  / len(intersection_dict) , \
#         # 100.*total_avg_fscore / len(intersection_dict), \
        
#         totalEntropy, \
#         corss_entropy, \
#         DKL, \
#         MI, \
#         VI, \
#         entropy_pred, \
#         entropy_pred_all, \
#         100 * purity, \
#         100 * all_missing_genes / totalGenes, \
#         100. * len(unidentified_refOGs) / total_RefOGs, \
#         total_overlap_predogs - len(unique_total_overlap_predogs), \
#         total_avg_distance,
#         distance_overlap,
#         distance_all
#     ]
    
#     return scores, matched_df



def calculate_benchmarks_pairwise(ref_ogs: Dict[str, str], 
                                  pred_ogs: Dict[str, str], 
                                  uncert_genes: Optional[Dict] = None, 
                                  q_even: bool = True, 
                                  q_remove_uncertain: bool = True,
                                  ):
    referenceOGs = ref_ogs
    predictedOGs = pred_ogs
    totalFP = 0.
    totalFN = 0.
    totalTP = 0.
    totalGroundTruth = 0.
    n_exact = 0
    n_splits = []
    # so as not to count uncertain genes either way remove them from the 
    # expected and remove them from any predicted OG (as though they never existed!)
    totalGenes = 0.

    for refog_key in referenceOGs:
        refOg = referenceOGs[refog_key]
        if uncert_genes is not None:
            uncert = uncert_genes[refog_key]

        thisFP = 0.
        thisFN = 0.
        thisTP = 0.
        this_split = 0
        if q_remove_uncertain and uncert_genes is not None:
            refOg = refOg.difference(uncert)

        nRefOG = len(refOg)
        totalGenes += nRefOG
        not_present = set(refOg)
        intersection = 0

        for predog_key in predictedOGs:
            predOg = predictedOGs[predog_key]
            overlap = len(refOg.intersection(predOg))

            intersection += overlap
            
            if overlap > 0:
                if q_remove_uncertain and uncert_genes is not None:
                    predOg = predOg.difference(uncert)   # I.e. only discount genes that are uncertain w.r.t. this RefOG
                overlap = len(refOg.intersection(predOg))

            if overlap > 0:
                this_split += 1
                not_present = not_present.difference(predOg)
                
                thisTP += overlap * (overlap - 1)/2    # n-Ch-2
                thisFP += overlap * (len(predOg) - overlap)
                thisFN += (nRefOG - overlap) * overlap
                
        # print(intersection_dict[refognum])

        # Are FNs more from splintered OGs or missing genes?
        # print("%f\t%f" % (thisFN/2./(nRefOG-1), len(not_present)*(nRefOG-1)/2./(nRefOG-1))) 
        # finally, count all the FN pairs from those not in any predicted OG
        thisFN += len(not_present)*(nRefOG-1)
        # don't count 'orphan genes' as splits, it's more informative only to count 
        # clusters that this orthogroup has been split into. Recall already counts
        #  the missing genes, this should give a distinct measure.
        # this_split += len(not_present)     
        # All FN have been counted twice
        assert(thisFN % 2 == 0)
        n_splits.append(this_split)
        # print(this_split)      # Orthogroup fragments
        # print(len(not_present))  # Unclustered genes 
        thisFN /= 2 
        # sanity check
        nPairs1 = thisTP + thisFN
        nPairs2 = nRefOG * (nRefOG - 1) / 2
        if nPairs1 != nPairs2:
            print("ERROR: %d != %d" % (nPairs1,nPairs2))    

        totalGroundTruth += nPairs1
        if thisFN == 0 and thisFP == 0:
            n_exact += 1
        if q_even:
            N = float(len(refOg)-1)
            totalFN += thisFN/N
            totalFP += thisFP/N
            totalTP += thisTP/N
        else:
            totalFN += thisFN
            totalFP += thisFP
            totalTP += thisTP
    
    TP, FP, FN = (totalTP, totalFP, totalFN)
    if TP != 0 or FP !=0:
        pres = TP/(TP+FP)
    else:
        pres = 0.0
    if TP != 0 or FN !=0:
        recall = TP/(TP+FN)
    else:
        recall = 0.0
    if pres !=0 or recall !=0:
        f = 2*pres*recall/(pres+recall)
    else:
        f = 0.0
    
    return TP, FP, FN, recall, pres, f




def read_refogs(d_refogs: str, check_genes: bool = True) -> Dict[str, str]:
    
    n_expected = 1945
    n_genes = 0
    refogs = {}
    try:
        if os.path.isdir(d_refogs):
            
            for file in pathlib.Path(d_refogs).iterdir():
                if not os.path.isfile(file):
                    continue
                if file.name.startswith("."):
                    continue
                key = file.name.rsplit(".", 1)[0]
                with open(file, 'r') as infile:
                    refogs[key] = set([g.rstrip() for g in infile.readlines()])
                    n_genes += len(refogs[key])
        elif os.path.isfile(d_refogs):
            with open(d_refogs, "r") as infile:
                for l in infile:
                    if l.startswith("#") or "Orthogroup" in l:
                        continue
                    key, t = l.strip().split(None, 1)
                    genes = re.split(delims, t)
                    genes = [g.strip() for g in genes if g.strip() != "" ]
                    # genes = set([g for g in genes if g != ""])
                    refogs[key] = set(genes)

    except:
        print("ERROR: RefOG file not found in: %s" % d_refogs)
        sys.exit(1)

    if check_genes:
        assert n_expected == n_genes, "ERROR: There are genes missing from the RefOG files. Found %d, expected %d" % (n_genes, n_expected)
    
    return refogs


def filter_orthogroups(ref_ogs: Dict[str, Set], pred_ogs: Dict[str, Set]):

    ref_ogs_set = set()
    for og in ref_ogs:
        ref_ogs_set = ref_ogs_set | set(ref_ogs[og])
    
    for pog in pred_ogs:
        diff = pred_ogs[pog] - ref_ogs_set
        if len(diff) > 0:
            pred_ogs[pog].difference_update(diff)

    return pred_ogs


def local_scores(ref_ogs, 
                 V_prime,
                 weighted_recall_dict,
                 weighted_precision_dict, 
                 weighted_f1score_dict, 
                 weighted_fowlkes_mallows_dict,  
                 weighted_JI_dict, 
                 weighted_dissimilarity_dict, 
                 weighted_distance_dict,  
                 entropy_dict,
                 missing_genes_dict, 
                 missing_genes_count_dict, 
                 missing_genes_proportion_dict,
                 effective_size_precision_weighted_dict,
                 effective_size_JI_weighted_dict,
                 effective_size_JI_refog_weighted_dict,
                 all_score_dict
                ):
    round_precision = 1
    data_dict = [
        (
            refog_key,
            len(refog),
            len(V_prime[refog_key]),
            np.round(100. * weighted_recall_dict[refog_key], round_precision),
            np.round(100. * weighted_precision_dict[refog_key], round_precision),
            np.round(100. * weighted_f1score_dict[refog_key], round_precision),
            np.round(100. * weighted_fowlkes_mallows_dict[refog_key], round_precision),
            np.round(100. * weighted_JI_dict[refog_key], round_precision),
            np.round(100. * weighted_dissimilarity_dict[refog_key], round_precision), 
            np.round(weighted_distance_dict[refog_key], 2), 
            np.round(entropy_dict[refog_key], 2), 
            missing_genes_count_dict[refog_key],
            np.round(100. * missing_genes_proportion_dict[refog_key], round_precision), 
            int(effective_size_precision_weighted_dict[refog_key]),
            int(effective_size_JI_weighted_dict[refog_key]),
            int(effective_size_JI_refog_weighted_dict[refog_key]),
            missing_genes_dict[refog_key],
        )
        for refog_key, refog in ref_ogs.items()
    ]
    

    colnames = [
        "RefOGs",
        "RefOG_Size", 
        "nPredictedOGs", 
        "avg_Recall (%)",
        "avg_Precision (%)", 
        "avg_F1-score (%)", 
        "avg_Fowlkes_Mallows (%)", 
        "avg_Jaccard_Index (%)", 
        "avg_Dissimilarity (%)",
        "avg_Distance",
        "Entropy", 
        "TotalMissingGenes", 
        "MissingGenes (%)",
        "EffectiveSize (precision_weighted)", 
        "EffectiveSize (JI_weighted)", 
        "EffectiveSize (JI_refog_weighted)",         
        "Missing_Genes",
    ]                   
            
    df = pd.DataFrame.from_records(data_dict, columns=colnames)
    # print(matched_df[["RefOG_intersection_PredOG", "min_FN", "FP"]])
    df.sort_values(by=["avg_Recall (%)", 
                      "avg_Precision (%)",
                      "Entropy"], \
                    inplace=True, ascending=[False, False, True])

    return df


def benchmark(ogs_filename, d_refogs, d_input, 
              skip_nrows,
              skip_ncols,
              d_uncertain=None, 
              check_ogs=True, 
              use_reg=True, 
              filter_genes=False,
              outgroups=None,
              additional_species=None,
              input_species = None):
    print("\nReading RefOGs from: %s" % d_refogs)

    ref_ogs = read_refogs(d_refogs, check_genes=False)
    ref_ogs_uncertain = None
    if d_uncertain is not None:
        ref_ogs_uncertain = read_refogs(d_uncertain, check_genes=False)
        ref_ogs_uncertain = {key: ref_ogs_uncertain[key] if key in ref_ogs_uncertain else set() for key in ref_ogs}
    

    # exp_genes = get_expected_genes()
    print("\nReading predicted orthogroups from: %s" % ogs_filename)
    # if q_basic_read: # not sure when to use this
    #     n_col_skip = get_n_col_skip(ogs_filename)
    #     pred_ogs = read_orthogroups(ogs_filename, exp_genes, n_col_skip)
    # else:

    pred_ogs = read_orthogroups_smart(ogs_filename, 
                                      skip_nrows, 
                                      skip_ncols, 
                                      use_reg=use_reg, 
                                      outgroups=outgroups)

    if filter_genes:
        pred_ogs = filter_orthogroups(ref_ogs, pred_ogs)
    
    if check_ogs:
        exp_genes, species_of_interest = get_expected_genes(d_input, 
                                                            outgroups, 
                                                            additional_species,
                                                            input_species)
        check_orthogroups(pred_ogs, exp_genes)

    
    V_prime = sf.V_raw_to_V_prime(ref_ogs, pred_ogs)
    V = sf.V_prime_to_V(V_prime)

    N = np.sum([len(refog) for refog in ref_ogs.values()])
    M = np.sum([len(predog) for predog in V.values()])
    M_raw = np.sum([len(predog) for predog in pred_ogs.values()])

    print("\nCalculating benchmarks:")
    print("%d  Number of RefOGs" % len(ref_ogs))
    print("%d  Number of genes in RefOGs" % N)
    print("%d  Number of PredOGs" % len(pred_ogs))
    print("%d  Number of genes in PredOGs (overlap)" % M)
    print("%d  Number of genes in PredOGs" % M_raw)
    print()
    
    # print(os.path.basename(ogs_filename))
    gene_pair_TP, gene_pair_FP, gene_pair_FN, \
    gene_pair_recall, gene_pair_precision, gene_pair_f1score = calculate_benchmarks_pairwise(ref_ogs, 
                                                                                            pred_ogs, 
                                                                                            ref_ogs_uncertain)
    
    print("%d  Correct gene pairs" % gene_pair_TP)
    print("%d  False Positives gene pairs" % gene_pair_FP)
    print("%d  False Negatives gene pairs\n" % gene_pair_FN)
    
    print("%0.1f%%  Gene Pair Recall" % (100. * gene_pair_recall))
    print("%0.1f%%  Gene Pair Precision" % (100. * gene_pair_precision))
    print("%0.1f%%  Gene Pair F1-score\n" % (100. * gene_pair_f1score))    
    
    missing_predogs_list = sf.missing_predogs(V_prime)
    print("%0.1f%%  Missing PrefOGs" % (100. * len(missing_predogs_list) / len(ref_ogs)))

    total_missing_genes, total_missing_genes_set, \
    missing_genes_dict, missing_genes_count_dict, \
    missing_genes_proportion_dict = sf.check_missing_genes(ref_ogs, V_prime)

    print("%0.1f%%  Missing Genes" % (100. * total_missing_genes / N))

    fussion_refog_set, fussion_predog_dict = sf.fussion(V_prime, V)
    fussion_refog_score = 100. * len(fussion_refog_set) / len(ref_ogs)
    if len(V) != 0:
        fussion_predog_score = 100. * len(fussion_predog_dict) / len(V)
    else:
        fussion_predog_score = 0.0
    print("%0.1f%%  Fussion (RefOG)" % fussion_refog_score)
    print("%0.1f%%  Fussion (PredOG)" % fussion_predog_score)


    fission_refog_set, fission_predog_set = sf.fission(ref_ogs, V_prime)
    fission_refog_score = 100. * len(fission_refog_set) / len(ref_ogs)
    if len(V) != 0:
        fission_predog_score = 100. * len(fission_predog_set) / len(V)
    else:
        fission_predog_score = 0.0
    print("%0.1f%%  Fission (RefOG)" % fission_refog_score)
    print("%0.1f%%  Fission (PredOG)" % fission_predog_score)
    print()

    micro_recall, micro_precision, \
    micro_f1score, micro_fowlkes_mallows_index, \
    micro_jaccard_index, micro_dissimilarity, micro_distance = sf.micro_scores(ref_ogs, V_prime, N)
    
    print("%0.1f%%  micro Recall" % (100. * micro_recall))
    print("%0.1f%%  micro Precision" % (100. * micro_precision))
    print("%0.1f%%  micro F1-score" % (100. * micro_f1score))
    print("%0.1f%%  micro Fowlkes Mallows Index" % (100. * micro_fowlkes_mallows_index))
    print("%0.1f%%  micro Jaccard Index" % (100. * micro_jaccard_index))
    print("%0.1f%%  micro Disimilarity" % (100. * micro_dissimilarity))
    print("%0.2f   micro Distance" % (micro_distance))
    print()

    total_weighted_recall, macro_recall, weighted_recall_dict, \
    total_weighted_precision, macro_precision, weighted_precision_dict, \
    total_weighted_f1score, macro_f1score, weighted_f1score_dict, \
    total_weighted_fowlkes_mallows, macro_fowlkes_mallows, weighted_fowlkes_mallows_dict, \
    total_weighted_JI, macro_JI, weighted_JI_dict,  \
    total_weighted_dissimilarity, macro_dissimilarity, weighted_dissimilarity_dict, \
    total_weighted_distance, macro_distance, weighted_distance_dict, \
    macro_dissimilarity_distance, \
    all_score_dict, \
    all_TP_dict, \
    effective_size_precision_weighted_dict, \
    effective_size_JI_weighted_dict, \
    effective_size_JI_refog_weighted_dict = sf.macro_and_weighted_avg_scores(ref_ogs, V_prime, N)
    
    print("%0.1f%%  macro Recall" % (100. * macro_recall))
    print("%0.1f%%  macro Precision" % (100. * macro_precision))
    print("%0.1f%%  macro F1-score" % (100. * macro_f1score))
    print("%0.1f%%  macro Fowlkes Mallows Index" % (100. * macro_fowlkes_mallows))
    print("%0.1f%%  macro Jaccard Index" % (100. * macro_JI))
    print("%0.1f%%  macro Disimilarity" % (100. * macro_dissimilarity))
    print("%0.2f   macro Distance" % (macro_distance))
    print("%0.2f   macro Disimilarity Distance" % (macro_dissimilarity_distance))
    print()

    print("%0.1f%%  Weighted Avg Recall" % (100. * total_weighted_recall))
    print("%0.1f%%  Weighted Avg Precision" % (100. * total_weighted_precision))
    print("%0.1f%%  Weighted Avg F1-score" % (100. * total_weighted_f1score))
    print("%0.1f%%  Weighted Avg Fowlkes Mallows Index" % (100. * total_weighted_fowlkes_mallows))
    print("%0.1f%%  Weighted Avg Jaccard Index" % (100. * total_weighted_JI))
    print("%0.1f%%  Weighted Avg Disimilarity" % (100. * total_weighted_dissimilarity))
    print("%0.2f   Weighted Avg Distance" % (total_weighted_distance))
    print()

    total_entropy, entropy_dict = sf.entropy_score(ref_ogs, V_prime, N)
    print("%0.2f  Entropy" % (total_entropy))

    precision_weighted_KLD = sf.kl_divergence(ref_ogs, effective_size_precision_weighted_dict, N, M)
    JI_weighted_KLD = sf.kl_divergence(ref_ogs, effective_size_JI_weighted_dict, N, M)
    JI_refog_weighted_KLD = sf.kl_divergence(ref_ogs, effective_size_JI_refog_weighted_dict, N, M)

    print("%0.2f  Precision Weighted KLD" % (precision_weighted_KLD))
    print("%0.2f  Jaccard Index Weighted KLD" % (JI_weighted_KLD))
    print("%0.2f  Jaccard Index refOG Weighted KLD" % (JI_refog_weighted_KLD))
    print()

    # V_complete = sf.V_to_V_complete(V, total_missing_genes_set)
    # AMI, AVI, ARI = sf.prime_information(ref_ogs, V_complete, N)
    
    # # print("%0.2f  Reduced Mutual Information (global)" % (NRMI))
    # print("%0.2f  Adjusted Mutual Information (global)" % (AMI))
    # print("%0.2f  Adjusted Variation of Information (global)" % (AVI))
    # print("%0.2f  Ajusted Rand Index (global)" % (ARI))
    # print()

    local_score_df = local_scores(ref_ogs, 
                                V_prime,
                                weighted_recall_dict,
                                weighted_precision_dict, 
                                weighted_f1score_dict, 
                                weighted_fowlkes_mallows_dict,  
                                weighted_JI_dict, 
                                weighted_dissimilarity_dict, 
                                weighted_distance_dict,
                                entropy_dict,  
                                missing_genes_dict, 
                                missing_genes_count_dict, 
                                missing_genes_proportion_dict,
                                effective_size_precision_weighted_dict,
                                effective_size_JI_weighted_dict,
                                effective_size_JI_refog_weighted_dict,
                                all_score_dict
                                )

    measures = [
        # len(V),
        # len(pred_ogs),
        # M,
        # M_raw,
        # gene_pair_TP,
        # gene_pair_FP,
        # gene_pair_FN,
        # 100. * gene_pair_recall,
        # 100. * gene_pair_precision,
        # 100. * gene_pair_f1score,
        100. * len(missing_predogs_list) / len(ref_ogs),
        100. * total_missing_genes / N,
        fussion_refog_score,
        # fussion_predog_score,
        fission_refog_score,
        # fission_predog_score,
        # 100. * micro_recall,
        # 100. * micro_precision,
        # 100. * micro_f1score,
        # 100. * micro_fowlkes_mallows_index,
        # 100. * micro_jaccard_index,
        # 100. * micro_dissimilarity,
        # micro_distance,
        # 100. * macro_recall,
        # 100. * macro_precision,
        # 100. * macro_f1score,
        # 100. * macro_fowlkes_mallows,
        # 100. * macro_JI,
        # 100. * macro_dissimilarity,
        # macro_distance,
        # macro_dissimilarity_distance,
        100. * total_weighted_recall,
        100. * total_weighted_precision,
        # 100. * total_weighted_f1score,
        # 100. * total_weighted_fowlkes_mallows,
        # 100. * total_weighted_JI,
        # 100. * total_weighted_dissimilarity,
        # total_weighted_distance,
        total_entropy,
        # precision_weighted_KLD,
        # JI_weighted_KLD,
        # JI_refog_weighted_KLD,
        # AMI, 
        # AVI, 
        # ARI
    ]

    return measures, local_score_df, all_TP_dict


def Fail():
    sys.stderr.flush()
    print(traceback.format_exc())
    sys.exit(1)


def iter_dir(d: Optional[str] = None) -> Iterator[str]:
    if d is None: 
        Fail()
    with os.scandir(d) as entries:
        for entry in entries:
            if entry.is_file():
                yield entry.name

def curate_labels(label: str):

    label = COMPILE.sub(" ", label).split()
    label = "_".join(label)

    return label

def get_expected_genes(input_dir: str, 
                       outgroups: Optional[List[str]] = None, 
                       additional_specieces: Optional[List[str]] = None,
                       input_species: Optional[List[str]] = None, 
                       ) -> Tuple[Set[str], Set[str]]:

    all_genes = set()
    species_of_interest = set()
    for fn in iter_dir(input_dir):
        species = fn.split(".")[0]
        species = curate_labels(species)
        
        if outgroups is not None:
            if species.lower() in outgroups:
                continue

        if additional_specieces is not None:
            if species.lower() in additional_specieces:
                continue

        if input_species is not None and outgroups is None and additional_specieces is None:
            species_of_interest.add(species)
        species_of_interest.add(species)
        
        fpath = os.path.join(input_dir, fn)
        with open(fpath, 'r') as infile:
            for l in infile:
                if l.startswith(">"):
                    all_genes.add(l[1:].rstrip())

    # assert len(all_genes) == n_genes_total, f"species_of_interest {len(species_of_interest)} all_genes {len(all_genes)} vs. n_genes_tatal {n_genes_total}"

    return all_genes, species_of_interest


def GetDirectoryArgument(arg: str) -> str:
    directory = os.path.abspath(arg)
    if not os.path.exists(directory):
        print("Reference orthogroups directory doesn't exist: %s" % directory)
        Fail()
    if not os.path.isfile(directory) and directory[-1] != os.sep: 
        directory += os.sep
    return directory

def GetFileArgument(arg: str) -> str:
    file_path = os.path.abspath(arg)
    if not os.path.isfile(file_path):
        print("Orthogroups file doesn't exist: %s" % file_path)
        Fail()
    return file_path

def save_all_TPs(all_TP_dict, filename):
    
    file_path = filename + "_overlap" + ".tsv"
    with open(file_path, "w") as writer:
        for refog_key, overlap_dict in all_TP_dict.items():
            
            writer.write(refog_key + "    ")
            predog_overlaps = []
            for predog_key, overlaps in overlap_dict.items():
                predog_overlaps.append(predog_key + ":" + overlaps)
            predog_overlaps_str = "    ".join(predog_overlaps)
            writer.write(predog_overlaps_str + "\n")




if __name__ == "__main__":
    # python3 benchmark.py -f test_files  -rd RefOGs
    # python3 benchmark.py  -rd simulated_data_1k/Simulated_Orthogroups_1k.txt -f simulation_1k  --skip-ncols "OF2:3,OF3:3,Broccoli:1,SP2:1,proteinortho:3,swiftortho:0,OrthoMCL:1" --skip-nrows "OF2:1,OF3:1,Broccoli:1,SP2:1,proteinortho:1,swiftortho:0,OrthoMCL:0" -c -reg
    parser = argparse.ArgumentParser()
    parser.add_argument("-f", help="File containing orthogroups.")
    parser.add_argument("--skip-ncols", 
                        default="OF2:3,OF3:3,Broccoli:1,SP2:1,proteinortho:3,swiftortho:0,OrthoMCL:1:Hieranoid:1", 
                        help="Skip the number of columns")
    parser.add_argument("--skip-nrows", 
                        default="OF2:1,OF3:1,Broccoli:1,SP2:1,proteinortho:1,swiftortho:0,OrthoMCL:0:Hieranoid:0", 
                        help="Skip the number of rows")
    parser.add_argument("-rd", default="RefOGs", help="Reference orthogroups directory")
    parser.add_argument("-d", default="orthobench_proteomes", help="Data directory")
    parser.add_argument("--low-uncertainty-assignment", help="low uncertainty assignment in refOGs")
    parser.add_argument("--additional-species", help="Additional species to help identification")
    parser.add_argument("--outgroups",  help="Outgroups")
    parser.add_argument("--input-species", help="Species of insterest")
    parser.add_argument("--output", help="Output file name")
    parser.add_argument("-ud", help="Low certainty reference orthogroups directory")
    parser.add_argument("-c", action="store_false", help="Check the query orthogroups")
    parser.add_argument("-filter", action="store_true", help="Filter the genes in the orthogroups")
    parser.add_argument("-reg", action="store_false", help="Using regex in refOGs")
    parser.add_argument("-b", "--basic", action="store_true", help="Basic delimited orthogroup file reader")
    args = parser.parse_args()
    print(args)

    if args.f is None:
        print("Missing query orthogroups file!\n")
        Fail()

    if os.path.isdir(args.f):
        og_file_dir = GetDirectoryArgument(args.f)
        og_files_dict = {fn.rsplit(".", 1)[0]: os.path.join(og_file_dir, fn) for fn in iter_dir(og_file_dir)}
    elif os.path.isfile(args.f):
        og_files_dict = {args.f.rsplit(".", 1)[0]: GetFileArgument(args.f)}

    print(args.rd)
    if os.path.isdir(args.rd):
        d_refogs = GetDirectoryArgument(args.rd)
        ref_dir = os.path.basename(d_refogs.rstrip(os.sep))
    elif os.path.isfile(args.rd):
        d_refogs = GetFileArgument(args.rd)
        ref_dir = os.path.basename(os.path.dirname(d_refogs))
    
    d_input = GetDirectoryArgument(args.d)

    if args.output is None:
        d_output_dir = os.path.join(os.path.dirname(d_refogs.rstrip(os.sep)), f"{args.f}_outputs")
    else:
        d_output_dir = GetDirectoryArgument(args.output)

    if not os.path.exists(d_output_dir):
        os.makedirs(d_output_dir, exist_ok=True)

    d_output_dict = {}
    for fn in og_files_dict:
        d_output_dict[fn] =  os.path.join(d_output_dir, fn + ".tsv")
    
    if args.ud is not None:
        u_refogs = GetDirectoryArgument(args.ud)
        args.outgroups = args.outgroups.strip().split(",")

    elif args.rd == "RefOGs":
        args.low_uncertainty_assignment = GetDirectoryArgument(ref_dir) + "low_certainty_assignments" + os.sep
        args.outgroups = ["Mnemiopsis leidyi", "Trichoplax adhaerens", "Nematostella vectensis"]
        args.additional_species = ["Branchiostoma lanceolatum", "Schistosoma mansoni"]

    if args.outgroups is not None:
        args.outgroups = [curate_labels(outgroup).lower() for outgroup in args.outgroups]

    if args.additional_species is not None:
        args.additional_species = [curate_labels(additional_species).lower() for additional_species in args.additional_species]
    
    if args.input_species is not None:
        args.input_species = [curate_labels(input_species).lower() for input_species in args.input_species]

    if len(args.skip_ncols.strip().split(",")) > 1:
        skip_ncols_dict = {item.split(":")[0].lower(): int(item.split(":")[1]) for item in args.skip_ncols.strip().split(",")}
    else:
        skip_ncols_dict = {args.skip_ncols.split(":")[0].lower(): int(args.skip_ncols.split(":")[1])}
    

    if len(args.skip_nrows.strip().split(",")) > 1:
        skip_nrows_dict = {item.split(":")[0].lower(): int(item.split(":")[1]) for item in args.skip_nrows.strip().split(",")}
    else:
        skip_nrows_dict = {args.skip_nrows.split(":")[0].lower(): int(args.skip_nrows.split(":")[1])}
    
    # old_scores_dict = {}
    new_scores_dict = {}
    overlapped_genes_df_dict = {}
    other_info_dict = {}
    missing_genes_dict = {}
   
    # fig, ax = plt.subplots(figsize=(10, 6)) 
    # colors = sns.color_palette("husl", len(og_files_dict))
    n_rows, n_cols = 3, 3
    fig, axs = plt.subplots(n_rows, n_cols, figsize=(15, n_rows * 5))  
    for ii, (fn, og_file) in enumerate(og_files_dict.items()):
        print()
        print(f"------------------------- {fn} ------------------------------")
        method = fn.split("_", 1)[0].lower()
        skip_ncols = skip_ncols_dict.get(method)
        skip_nrows = skip_nrows_dict.get(method)

        new_scores, local_score_df, all_TP_dict = benchmark(og_file, d_refogs, d_input, 
                                                            skip_nrows,
                                                            skip_ncols,
                                                            d_uncertain=args.low_uncertainty_assignment, 
                                                            check_ogs=args.c,
                                                            use_reg=args.reg, 
                                                            filter_genes=args.filter,
                                                            outgroups=args.outgroups,
                                                            additional_species=args.additional_species,
                                                            input_species=args.input_species)
        print(fn)
        local_score_df.to_csv(d_output_dict[fn], sep='\t', index=False, header=True)
        save_all_TPs(all_TP_dict, fn)
        new_scores_dict[fn] = new_scores

    colnames = [
        # "nPredOGs (overlap)",
        # "nPredOGs",
        # "totalGenes PredOGs (overlap)",
        # "totalGenes PredOGs",
        # "TP Gene Pairs",
        # "FP Gene Pairs",
        # "FN Gene Pairs",
        # "Gene Pair Recall",
        # "Gene Pair Precision",
        # "Gene Pair F1-score",
        "Missing PredOGs",
        "Missing Genes",
        "Fussion (RefOG)",
        # "Fussion (PredOG)",
        "Fission (RefOG)",
        # "Fission (PredOG)",
        # "micro Recall",
        # "micro Precision",
        # "micro F1-score",
        # "micro Fowlkes Mallows Index",
        # "micro Jaccard Index",
        # "micro Disimilarity",
        # "micro Distance",
        # "macro Recall",
        # "macro Precision",
        # "macro F1-score",
        # "macro Fowlkes Mallows Index",
        # "macro Jaccard Index",
        # "macro Disimilarity",
        # "macro Distance",
        # "macro Disimilarity Distance",
        "Weighted Avg Recall",
        "Weighted Avg Precision",
        # "Weighted Avg F1-score",
        # "Weighted Avg Fowlkes Mallows Index",
        # "Weighted Avg Jaccard Index",
        # "Weighted Avg Disimilarity",
        # "Weighted Avg Distance",
        "Entropy",
        # "Precision Weighted KLD",
        # "Jaccard Index Weighted KLD",
        # "Jaccard Index refOG Weighted KLD",
        # "Adjusted Mutual Information (global)",
        # "Adjusted Variation of Information (global)",   
        # "Ajusted Rand Index (global)"                    
    ]
    new_scores_df = pd.DataFrame.from_dict(new_scores_dict, 
                                           orient="index", 
                                           columns=colnames).round(3)

    
    # old_scores_df.reset_index(inplace=True)
    # old_scores_df.rename(columns={"index": "Filename"}, inplace=True)
    new_scores_df.reset_index(inplace=True)
    new_scores_df.rename(columns={"index": "Filename"}, inplace=True)

    # print(old_scores_df)
    print(new_scores_df)

    # old_scores_df.to_csv(os.path.join(d_output_dir, "_".join((args.f, 'old_scores.csv'))), index=False, header=True)  
    new_scores_df.to_csv(os.path.join(d_output_dir, "_".join((args.f, 'new_scores.csv'))), index=False, header=True)
    


        

        
