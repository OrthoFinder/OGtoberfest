import numpy as np
import scipy.stats as ss
from scipy.spatial import distance
from typing import Dict, List, Set, Optional
from ogtoberfest.utils import util
from ogtoberfest.orthogroups import scorefuncs as sf
from ogtoberfest.orthogroups import orthogroups_analyser as opa
import itertools

def combine_scores(scores_list, 
                   score_names, 
                   precision=3, 
                   avg_method="RMS Score",
                   weights=None):
    if weights is None:
        weights = {
            name: 1 / len(score_names)
            for name in score_names
        }
    else:
        weights = {
            
            
        }
    score_dict = dict(zip(score_names, scores_list))
    if avg_method == "Avg Score":
        combined_scores = 0.0
        for name, score in score_dict.items():
            if "recall" in name.lower() \
                or "precision" in name.lower() \
                    or "f1-score" in name.lower() \
                        or "index" in name.lower():
                combined_scores += weights[name] * score / 100
            elif "entropy" in name.lower() or "disimilarity" in name.lower():
                combined_scores += weights[name] * (1 - score)
            else:
                combined_scores += weights[name] * (1 - score / 100)

        # combined_scores /= len(scores_list)
        combined_scores /= np.sum([*weights.values()])

    elif avg_method == "RMS Score":
        combined_scores = 0.0
        for name, score in score_dict.items():
            if "recall" in name.lower() \
                or "precision" in name.lower() \
                    or "f1-score" in name.lower() \
                        or "index" in name.lower():
                combined_scores += weights[name] * (score / 100) ** 2
            elif "entropy" in name.lower() or "disimilarity" in name.lower():
                combined_scores += weights[name] * (1 - score) ** 2
            else:
                combined_scores += weights[name] * (1 - score / 100) ** 2

        # combined_scores /= len(scores_list)
        combined_scores /= np.sum([*weights.values()])
        combined_scores = np.sqrt(combined_scores)
    return np.round(combined_scores, precision)


def rank_score(global_score_dict, score_names, precision=3, rank_method="min"):
    method, scores = zip(*global_score_dict.items())
    filtered_score = []
    for score in scores:
        score_dict = dict(zip(score_names, score))
        score_list = []
        for name, score in score_dict.items():
            if "recall" in name.lower() \
                or "precision" in name.lower() \
                    or "f1-score" in name.lower() \
                        or "index" in name.lower():
                score_list.append(1 - score / 100)
            elif "entropy" in name.lower() \
                or "disimilarity" in name.lower() \
                or "size" in name.lower() \
                or "time" in name.lower():
                score_list.append(score)
            else:
                score_list.append(score / 100)

        filtered_score.append(score_list)

    score_arr = np.array(filtered_score)
    rank_arr = ss.rankdata(score_arr, method=rank_method, axis=0)

    mean_rank = np.mean(rank_arr, axis=1).round(precision)
    for i, score in enumerate(scores):
        score.append(mean_rank[i])

    ranked_global_score_dict = dict(zip(method, scores))
    global_scores_rank_dict = dict(zip(method, rank_arr))
    return ranked_global_score_dict, global_scores_rank_dict

def z_score(global_score_dict, score_names, precision=3):
    method, scores = zip(*global_score_dict.items())
    filtered_score = []
    for score in scores:
        score_dict = dict(zip(score_names, score))
        score_list = []
        for name, score in score_dict.items():
            if "recall" in name.lower() \
                or "precision" in name.lower() \
                    or "f1-score" in name.lower() \
                        or "index" in name.lower():
                score_list.append(score / 100)
            elif "entropy" in name.lower() or "disimilarity" in name.lower():
                score_list.append(1 - score)
            else:
                score_list.append(1 - score / 100)

        filtered_score.append(score_list)

    score_arr = np.array(filtered_score)
    z_score_arr = ss.zscore(score_arr, axis=0, ddof=1, nan_policy='omit')
    z_score_arr = np.nan_to_num(z_score_arr)
    z_score_sum = np.sum(z_score_arr, axis=1).round(3)
    for i, score in enumerate(scores):
        score.append(z_score_sum[i])

    ranked_global_score_dict = dict(zip(method, scores))
    return ranked_global_score_dict

def get_scores(
        ref_ogs: Dict[str, Set[str]],
        pred_ogs: Dict[str, Set[str]],
        precision: int,
        nthreads: int
    ):

    V_prime = sf.V_raw_to_V_prime(ref_ogs, pred_ogs)
    V = sf.V_prime_to_V(V_prime)

    N = np.sum([len(refog) for refog in ref_ogs.values()])
    M = np.sum([len(predog) for predog in V.values()])
    M_raw = np.sum([len(predog) for predog in pred_ogs.values()])

    refog_species_dict = opa.get_species_dict(ref_ogs)
    predog_species_dict_raw = opa.get_species_dict(pred_ogs)
    predog_species_dict_overlap = opa.get_species_dict(V)
    predog_species_dict_prime = opa.get_predog_species_dict(V_prime)
    #
    num_species_refog_dict = opa.get_num_species_dict(refog_species_dict)
    num_species_predog_dict_raw = opa.get_num_species_dict(predog_species_dict_raw)
    num_species_predog_dict_overlap = opa.get_num_species_dict(predog_species_dict_overlap)
    num_species_predog_dict_prime = opa.get_predog_num_species_dict(predog_species_dict_prime)
    
    missing_species_dict, missing_species_count_dict, missing_species_percent_dict = \
        sf.check_missing_species(refog_species_dict, predog_species_dict_prime, precision)
    missing_species_refogs_num = np.sum([
        1 if num != 0 else 0
        for num in missing_species_count_dict.values()
    ])

    missing_species_refogs = missing_species_refogs_num / len(ref_ogs)
    total_num_species_refog = len(set(util.flatten_list_of_list([*refog_species_dict.values()])))

    total_num_species_predog_raw = len(
        set(util.flatten_list_of_list([*predog_species_dict_raw.values()]))
    )

    total_num_species_predog_overlap = len(
        set(util.flatten_list_of_list([*predog_species_dict_overlap.values()]))
    )
    if nthreads == 1:
        print("\nCalculating benchmarks:")
        print("%d  Number of RefOGs" % len(ref_ogs))
        print("%d  Number of species in RefOGs" % total_num_species_refog)
        print("%d  Number of genes in RefOGs" % N)
        print("%d  Number of PredOGs" % len(pred_ogs))
        print("%d  Number of species in PredOGs" % total_num_species_predog_raw)
        print("%d  Number of genes in PredOGs" % M_raw)
        print("%d  Number of PredOGs (overlap)" % len(V))
        print("%d  Number of species in PredOGs (overlap)" % total_num_species_predog_overlap)
        print("%d  Number of genes in PredOGs (overlap)" % M)

        print()

    (
        gene_pair_TP,
        gene_pair_FP,
        gene_pair_FN,
        gene_pair_recall,
        gene_pair_precision,
        gene_pair_f1score,
    ) = sf.calculate_benchmarks_pairwise(ref_ogs, pred_ogs)
    
    # print("%d  Correct gene pairs" % gene_pair_TP)
    # print("%d  False Positives gene pairs" % gene_pair_FP)
    # print("%d  False Negatives gene pairs\n" % gene_pair_FN)

    # print("%0.1f%%  Gene Pair Recall" % (100. * gene_pair_recall))
    # print("%0.1f%%  Gene Pair Precision" % (100. * gene_pair_precision))
    # print("%0.1f%%  Gene Pair F1-score\n" % (100. * gene_pair_f1score))


    missing_predogs_list = sf.missing_predogs(V_prime)
    missing_predogs = len(missing_predogs_list) / len(ref_ogs)

    if nthreads == 1:
        print("%0.1f%%  Missing RefOGs" % (100.0 * missing_predogs))
        print("%0.1f%%  Missing Speciess" % (100.0 * missing_species_refogs))

    (
        total_missing_genes,
        total_missing_genes_set,
        missing_genes_dict,
        missing_genes_count_dict,
        missing_genes_per_dict,
    ) = sf.check_missing_genes(ref_ogs, V_prime, precision)

    missing_genes = total_missing_genes / N
    if nthreads == 1:
        print("%0.1f%%  Missing Genes" % (100.0 * missing_genes))

    fusion_refog_set, fusion_predog_dict, fusion_refog_bool_dict = sf.fusion(ref_ogs, V_prime, V)
    fusion_refog_score = len(fusion_refog_set) / len(ref_ogs)
    if len(V) != 0:
        fusion_predog_score = len(fusion_predog_dict) / len(V)
    else:
        fusion_predog_score = 0.0

    if nthreads == 1:
        print("%0.1f%%  RefOG Fusions" % (100.0 * fusion_refog_score))
    # print("%0.1f%%  PredOG Fusion" % (100.0 * fusion_predog_score))

    fission_refog_set, fission_predog_set, fission_refog_bool_dict = sf.fission(ref_ogs, V_prime)
    fission_refog_score = len(fission_refog_set) / len(ref_ogs)
    if len(V) != 0:
        fission_predog_score = len(fission_predog_set) / len(V)
    else:
        fission_predog_score = 0.0
    if nthreads == 1:
        print("%0.1f%%  RefOG Fissions" % (100.0 * fission_refog_score))
    # print("%0.1f%%  PredOG Fissions" % (100.0 * fission_predog_score))
    # print()

    micro_recall, micro_precision, \
    micro_f1score, micro_fowlkes_mallows_index, \
    micro_jaccard_index, micro_dissimilarity, micro_distance = sf.micro_scores(ref_ogs, V_prime, N)

    # print("%0.1f%%  micro Recall" % (100. * micro_recall))
    # print("%0.1f%%  micro Precision" % (100. * micro_precision))
    # print("%0.1f%%  micro F1-score" % (100. * micro_f1score))
    # print("%0.1f%%  micro Fowlkes Mallows Index" % (100. * micro_fowlkes_mallows_index))
    # print("%0.1f%%  micro Jaccard Index" % (100. * micro_jaccard_index))
    # print("%0.1f%%  micro Disimilarity" % (100. * micro_dissimilarity))
    # print("%0.2f   micro Distance" % (micro_distance))
    # print()

    (
        avg_weighted_recall,
        macro_recall,
        weighted_recall_dict,
        avg_weighted_precision,
        macro_precision,
        weighted_precision_dict,
        avg_weighted_f1score,
        macro_f1score,
        weighted_f1score_dict,
        avg_weighted_fowlkes_mallows,
        macro_fowlkes_mallows,
        weighted_fowlkes_mallows_dict,
        avg_weighted_JI,
        macro_JI,
        weighted_JI_dict,
        avg_weighted_dissimilarity,
        macro_dissimilarity,
        weighted_dissimilarity_dict,
        avg_weighted_distance,
        macro_distance,
        weighted_distance_dict,
        macro_dissimilarity_distance,
        all_score_dict,
        all_TP_dict,
        effective_size_precision_weighted_dict,
        effective_size_JI_weighted_dict,
        effective_size_JI_refog_weighted_dict,
    ) = sf.macro_and_weighted_avg_scores(ref_ogs, V_prime, N, precision)

    fusion_refog_genes_dict = {
        refog_key: {
            predog_key: tp
            for predog_key, tp in tp_dict.items()
            if predog_key in fusion_predog_dict
        }
        if refog_key in fusion_refog_set else {}
        for refog_key, tp_dict in all_TP_dict.items()
    }

    fission_refog_genes_dict = {
        refog_key: {
            predog_key: tp
            for predog_key, tp in tp_dict.items()
        }
        if refog_key in fission_refog_set else {}
        for refog_key, tp_dict in all_TP_dict.items()
    }

    total_entropy, entropy_dict = sf.entropy_score(ref_ogs, V_prime, N)
    if nthreads == 1:
        # print("%0.1f%%  macro Recall" % (100. * macro_recall))
        # print("%0.1f%%  macro Precision" % (100. * macro_precision))
        # print("%0.1f%%  macro F1-score" % (100. * macro_f1score))
        # print("%0.1f%%  macro Fowlkes Mallows Index" % (100. * macro_fowlkes_mallows))
        # print("%0.1f%%  macro Jaccard Index" % (100. * macro_JI))
        # print("%0.1f%%  macro Disimilarity" % (100. * macro_dissimilarity))
        # print("%0.2f   macro Distance" % (macro_distance))
        # print("%0.2f   macro Disimilarity Distance" % (macro_dissimilarity_distance))
        # print()

        print("%0.1f%%  Weighted Avg Recall" % (100.0 * avg_weighted_recall))
        print("%0.1f%%  Weighted Avg Precision" % (100.0 * avg_weighted_precision))
        print("%0.1f%%  Weighted Avg F1-score" % (100.0 * avg_weighted_f1score))
        # print("%0.1f%%  Weighted Avg Fowlkes Mallows Index" % (100. * avg_weighted_fowlkes_mallows))
        # print("%0.1f%%  Weighted Avg Jaccard Index" % (100. * avg_weighted_JI))
        # print("%0.1f%%  Weighted Avg Disimilarity" % (100. * avg_weighted_dissimilarity))
        # print("%0.2f   Weighted Avg Distance" % (avg_weighted_distance))
        # print()
        print("%0.2f  Entropy" % (total_entropy))
        print()

    # precision_weighted_KLD = sf.kl_divergence(ref_ogs, effective_size_precision_weighted_dict, N, M)
    # JI_weighted_KLD = sf.kl_divergence(ref_ogs, effective_size_JI_weighted_dict, N, M)
    # JI_refog_weighted_KLD = sf.kl_divergence(ref_ogs, effective_size_JI_refog_weighted_dict, N, M)

    # print("%0.2f  Precision Weighted KLD" % (precision_weighted_KLD))
    # print("%0.2f  Jaccard Index Weighted KLD" % (JI_weighted_KLD))
    # print("%0.2f  Jaccard Index refOG Weighted KLD" % (JI_refog_weighted_KLD))
    

    global_stats_dict = {
        "nrefogs": len(ref_ogs),
        "n_species_refogs": total_num_species_refog,
        "n_genes_refogs": N,
        "npredogs": len(pred_ogs),
        "n_species_predogs": total_num_species_predog_raw,
        "n_genes_predogs": M_raw,
        "n_overlaped_species_predogs": total_num_species_predog_overlap,
        "n_overlaped_genes_predogs": M,
        "correct_gene_pairs": gene_pair_TP,
        "false_positive_gene_pairs": gene_pair_FP,
        "false_negative_gene_pairs": gene_pair_FN
    }

    local_stats_dict = {
        "n_species_refogs": num_species_refog_dict,
        "n_species_predogs": num_species_predog_dict_raw, 
        "n_overlaped_species_predogs": num_species_predog_dict_overlap,
        "n_species_predogs_refogs": num_species_predog_dict_prime,
    }

    local_score_dict = {
        "Entropy": entropy_dict,
        "Weighted Avg Recall": weighted_recall_dict,
        "Weighted Avg Precision": weighted_precision_dict,
        "Weighted Avg F1-score": weighted_f1score_dict,
        "Weighted Avg Fowlkes-Mallows Index": weighted_fowlkes_mallows_dict,
        "Weighted Avg Jaccard Index": weighted_JI_dict,
        "Weighted Avg Dissimilarity": weighted_dissimilarity_dict,
        "Missing Genes Percentage": missing_genes_per_dict,
        "Missing Species Percentage": missing_species_percent_dict,
    }

    local_properties_dict = {
        "Missing Species Count": missing_species_count_dict,
        "Missing Genes Count": missing_genes_count_dict,
        "Effective Size": effective_size_JI_weighted_dict,
        "Fusion Bool": fusion_refog_bool_dict,
        "Fission Bool": fission_refog_bool_dict,
    }

    global_score_dict = {
        "Missing RefOGs (%)": np.round(100.0 * missing_predogs, precision),
        "Missing Species (%)": np.round(100.0 * missing_species_refogs, precision),
        "Missing Genes (%)": np.round(100.0 * missing_genes, precision),
        "RefOG Fusions (%)": np.round(100.0 * fusion_refog_score, precision),
        "RefOG Fissions (%)": np.round(100.0 * fission_refog_score, precision),
        "PredOG Fusions (%)": np.round(100.0 * fusion_predog_score, precision),
        "PredOG Fissions (%)": np.round(100.0 * fission_predog_score, precision),
        "Recall": np.round(100.0 * avg_weighted_recall, precision),
        "Precision": np.round(100.0 * avg_weighted_precision, precision),
        "F1-score": np.round(100.0 * avg_weighted_f1score, precision),
        "Fowlkes-Mallows Index": np.round(100. * avg_weighted_fowlkes_mallows, precision),
        "Jaccard Index": np.round(100. * avg_weighted_JI, precision),
        "Disimilarity": np.round (100. * avg_weighted_dissimilarity, precision),
        "Entropy": np.round(total_entropy, precision),
        "Macro Recall": np.round(100. * macro_recall, precision),
        "Macro Precision": np.round(100. * macro_precision, precision),
        "Macro F1-score": np.round(100. * macro_f1score, precision),
        "Macro Fowlkes-Mallows Index": np.round(100. * macro_fowlkes_mallows, precision),
        "Macro Jaccard Index": np.round(macro_JI, precision),
        "Macro Disimilarity": np.round(100. * macro_dissimilarity, precision),
        "Micro Recall": np.round(100. * micro_recall, precision),
        "Micro Precision": np.round(100. * micro_precision, precision),
        "Micro F1-score": np.round(100. * micro_f1score, precision),
        "Micro Fowlkes-Mallows Index": np.round(100. * micro_fowlkes_mallows_index, precision),
        "Micro Jaccard Index": np.round(100. * micro_jaccard_index, precision),
        "Micro Disimilarity": np.round(100. * micro_dissimilarity, precision),
        "Gene Pair Recall": np.round(100. * gene_pair_recall, precision),
        "Gene Pair Precision": np.round(100. * gene_pair_precision, precision),
        "Gene Pair F1-score": np.round(100. * gene_pair_f1score, precision),
    }

    other_info = {
        "Missing Genes": missing_genes_dict,
        "Missing Genes Set": total_missing_genes_set,
        "PredOGs Info": all_score_dict,
        "Overlapped Genes": all_TP_dict,
        "Missing Species": missing_species_dict,
        "Fusion Genes": fusion_refog_genes_dict,
        "Fission Genes": fission_refog_genes_dict,
        "Overlap PredOGs": V,
    }

    return global_stats_dict, local_stats_dict, global_score_dict, \
        local_properties_dict, local_score_dict, other_info



def corr_vi_analysis(
        method_predogs_dict, 
        refogs_ngenes_dict,
        filewriter,
        precision = 3,
):

    methods = [*method_predogs_dict.keys()]
    methods_product = [*itertools.product(methods, methods)]
    VI_dict = {m1: {m2: {} for m2 in methods} for m1 in methods}
    avg_VI_dict = {m1: {m2: {} for m2 in methods} for m1 in methods}
    for m1, m2 in methods_product:
        VI_dict[m1][m2] = {}
        VI_list = [] 
        nrefog_list = []
        for refog_key, nrefog in refogs_ngenes_dict.items():
            predog_1 = method_predogs_dict[m1][refog_key]
            predog_2 = method_predogs_dict[m2][refog_key]
            list_products = [*itertools.product(predog_1.items(), predog_2.items())]
            VI = sf.variation_of_information(
                predog_1, 
                predog_2, 
                nrefog, 
                list_products
            )

            VI_dict[m1][m2][refog_key] = np.round(VI, precision)
            VI_list.append(VI)
            nrefog_list.append(nrefog)
        VI_arr = np.array(VI_list)
        nrefog_arr = np.array(nrefog_list)
        avg_VI_dict[m1][m2] = np.round(np.average(VI_arr, weights=nrefog_arr), precision)
        

    filewriter.save_global_VI_scores(
        avg_VI_dict, 
    )
    
    filewriter.save_local_VI_scores(
        VI_dict
    )

def corr_dist_analysis(
    local_scores_dict,
    local_scores,
    filewriter,
    metric,
    precision = 3,
):
    
    methods = [*local_scores_dict.keys()]

    metric_dict = {
        score_name: np.array([
            [*local_scores_dict[m1][score_name].values()] for m1 in methods
            ])
        for score_name in local_scores
    }

    for score, score_arr in metric_dict.items():
        dist_arr = distance.squareform(distance.pdist(score_arr, metric)).round(precision)
        filewriter.save_dist_metrics(metric, score, np.nan_to_num(dist_arr), methods)


def distribution_analyser(ref_dict, compare_dict, metric):

    dist_funcs = {name: getattr(distance, name) for name in dir(distance) if "__" not in name}

    ref_list = []
    compare_list = []
    for ref_key, ref_val in ref_dict.items():
        ref_list.append(ref_val)
        compare_val = compare_dict.get(ref_key)
        compare_val = compare_val if compare_val is not None else 0
        compare_list.append(compare_val)

    dfunc = dist_funcs.get(metric)
    dist_val = np.nan
    if dfunc is not None:
        dist_val = dfunc(ref_list, compare_list)
    return dist_val
        

    
    
        


            








