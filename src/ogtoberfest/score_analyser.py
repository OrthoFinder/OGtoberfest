import numpy as np
import scipy.stats as ss
from typing import Dict, List, Set, Optional
from ogtoberfest import utils
from ogtoberfest import scorefuncs as sf
from ogtoberfest import orthogroups_analyser as opa


def combine_scores(scores_list, score_names, precision=3, avg_method="mean"):

    score_dict = dict(zip(score_names, scores_list))
    if avg_method == "mean":
        combined_scores = 0.0
        for name, score in score_dict.items():
            if "recall" in name.lower() \
                or "precision" in name.lower() \
                    or "f1-score" in name.lower() \
                        or "index" in name.lower():
                combined_scores += score / 100
            elif "entropy" in name.lower() or "disimilarity" in name.lower():
                combined_scores += 1 - score
            else:
                combined_scores += 1 - score / 100

        combined_scores /= len(scores_list)

    elif avg_method == "rms":
        combined_scores = 0.0
        for name, score in score_dict.items():
            if "recall" in name.lower() \
                or "precision" in name.lower() \
                    or "f1-score" in name.lower() \
                        or "index" in name.lower():
                combined_scores += (score / 100) ** 2
            elif "entropy" in name.lower() or "disimilarity" in name.lower():
                combined_scores += (1 - score) ** 2
            else:
                combined_scores += (1 - score / 100) ** 2

        combined_scores /= len(scores_list)
        combined_scores = np.sqrt(combined_scores)
    return np.round(combined_scores, precision)


def rand_score(global_score_dict, score_names, precision=3, rank_method="min"):
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
            elif "entropy" in name.lower() or "disimilarity" in name.lower():
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

    total_num_species_refog = len(set(utils.flatten_list_of_list([*refog_species_dict.values()])))

    total_num_species_predog_raw = len(
        set(utils.flatten_list_of_list([*predog_species_dict_raw.values()]))
    )

    total_num_species_predog_overlap = len(
        set(utils.flatten_list_of_list([*predog_species_dict_overlap.values()]))
    )

    print("\nCalculating benchmarks:")
    # print("%d  Number of RefOGs" % len(ref_ogs))
    # print("%d  Number of species in RefOGs" % total_num_species_refog)
    # print("%d  Number of genes in RefOGs" % N)
    # print("%d  Number of PredOGs" % len(pred_ogs))
    # print("%d  Number of species in PredOGs" % total_num_species_predog_raw)
    # print("%d  Number of genes in PredOGs" % M_raw)
    # print("%d  Number of species in PredOGs (overlap)" % total_num_species_predog_overlap)
    # print("%d  Number of genes in PredOGs (overlap)" % M)

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
    print("%0.1f%%  Missing PrefOGs" % (100.0 * missing_predogs))

    (
        total_missing_genes,
        total_missing_genes_set,
        missing_genes_dict,
        missing_genes_count_dict,
        missing_genes_per_dict,
    ) = sf.check_missing_genes(ref_ogs, V_prime, precision)

    missing_genes = total_missing_genes / N
    print("%0.1f%%  Missing Genes" % (100.0 * missing_genes))

    fusion_refog_set, fusion_predog_dict = sf.fusion(V_prime, V)
    fusion_refog_score = len(fusion_refog_set) / len(ref_ogs)
    if len(V) != 0:
        fusion_predog_score = len(fusion_predog_dict) / len(V)
    else:
        fusion_predog_score = 0.0
    print("%0.1f%%  Fussion (RefOG)" % (100.0 * fusion_refog_score))
    # print("%0.1f%%  Fussion (PredOG)" % (100.0 * fusion_predog_score))

    fision_refog_set, fision_predog_set = sf.fision(ref_ogs, V_prime)
    fision_refog_score = len(fision_refog_set) / len(ref_ogs)
    if len(V) != 0:
        fision_predog_score = len(fision_predog_set) / len(V)
    else:
        fision_predog_score = 0.0
    print("%0.1f%%  Fission (RefOG)" % (100.0 * fision_refog_score))
    # print("%0.1f%%  Fission (PredOG)" % (100.0 * fision_predog_score))
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

    total_entropy, entropy_dict = sf.entropy_score(ref_ogs, V_prime, N)
    print("%0.2f  Entropy" % (total_entropy))

    # precision_weighted_KLD = sf.kl_divergence(ref_ogs, effective_size_precision_weighted_dict, N, M)
    # JI_weighted_KLD = sf.kl_divergence(ref_ogs, effective_size_JI_weighted_dict, N, M)
    # JI_refog_weighted_KLD = sf.kl_divergence(ref_ogs, effective_size_JI_refog_weighted_dict, N, M)

    # print("%0.2f  Precision Weighted KLD" % (precision_weighted_KLD))
    # print("%0.2f  Jaccard Index Weighted KLD" % (JI_weighted_KLD))
    # print("%0.2f  Jaccard Index refOG Weighted KLD" % (JI_refog_weighted_KLD))
    print()

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
        "Weighted Avg Fowlkes-Mallows Index":weighted_fowlkes_mallows_dict,
        "Weighted Avg Jaccard Index": weighted_JI_dict,
        "Weighted Avg Dissimilarity": weighted_dissimilarity_dict,
        "Effective Size": effective_size_JI_weighted_dict,
        "Missing Genes Count": missing_genes_count_dict,
        "Missing Genes Percentage": missing_genes_per_dict,
    }

    global_score_dict = {
        "Missing PredOGs": np.round(100.0 * missing_predogs, precision),
        "Missing Genes": np.round(100.0 * missing_genes, precision),
        "Fusion (RefOGs)": np.round(100.0 * fusion_refog_score, precision),
        "Fision (RefOGs)": np.round(100.0 * fision_refog_score, precision),
        "Fusion (PredOGs)": np.round(100.0 * fusion_predog_score, precision),
        "Fision (PredOGs)": np.round(100.0 * fision_predog_score, precision),
        "Weighted Avg Recall": np.round(100.0 * avg_weighted_recall, precision),
        "Weighted Avg Precision": np.round(100.0 * avg_weighted_precision, precision),
        "Weighted Avg F1-score": np.round(100.0 * avg_weighted_f1score, precision),
        "Weighted Avg Fowlkes-Mallows Index": np.round(100. * avg_weighted_fowlkes_mallows, precision),
        "Weighted Avg Jaccard Index": np.round(100. * avg_weighted_JI, precision),
        "Weighted Avg Disimilarity": np.round (100. * avg_weighted_dissimilarity, precision),
        "Avg Entropy": np.round(total_entropy, precision),
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
    predogs_info = {
        "Missing Genes": missing_genes_dict,
        "Missing Genes Set": total_missing_genes_set,
        "PredOGs Info": all_score_dict,
        "True Positives": all_TP_dict,
    }
    # scores = [
    #     np.round(100.0 * missing_predogs, precision),
    #     np.round(100.0 * missing_genes, precision),
    #     np.round(100.0 * fusion_refog_score, precision),
    #     np.round(100.0 * fision_refog_score, precision),
    #     np.round(100.0 * avg_weighted_recall, precision),
    #     np.round(100.0 * avg_weighted_precision, precision),
    #     np.round(total_entropy, 2),
    # ]

    # if isinstance(global_additional_scores, list):
    return global_stats_dict, local_stats_dict, global_score_dict, local_score_dict, predogs_info


def local_scores(
        ref_ogs,
        V_prime,
        weighted_recall_dict,
        weighted_precision_dict,
        weighted_f1score_dict,
        entropy_dict,
        missing_genes_dict,
        missing_genes_count_dict,
        missing_genes_proportion_dict,
        effective_size_precision_weighted_dict,
        effective_size_JI_weighted_dict,
    ):

    round_precision = 1
    data_dict = [
        (
            refog_key,
            len(refog),
            len(V_prime[refog_key]),
            np.round(100.0 * weighted_recall_dict[refog_key], round_precision),
            np.round(100.0 * weighted_precision_dict[refog_key], round_precision),
            np.round(100.0 * weighted_f1score_dict[refog_key], round_precision),
            np.round(entropy_dict[refog_key], 2),
            missing_genes_count_dict[refog_key],
            np.round(100.0 * missing_genes_proportion_dict[refog_key], round_precision),
            int(effective_size_precision_weighted_dict[refog_key]),
            int(effective_size_JI_weighted_dict[refog_key]),
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
        "Entropy",
        "TotalMissingGenes",
        "MissingGenes (%)",
        "EffectiveSize (JI_weighted)",
        "Missing_Genes",
    ]

    df = pd.DataFrame.from_records(data_dict, columns=colnames)
    # print(matched_df[["RefOG_intersection_PredOG", "min_FN", "FP"]])
    df.sort_values(
        by=["avg_Recall (%)", "avg_Precision (%)", "Entropy"],
        inplace=True,
        ascending=[False, False, True],
    )

    return df
