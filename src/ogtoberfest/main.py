import sys
from typing import Dict, List, Optional, Callable, Set
import re
import pathlib
import numpy as np
import pandas as pd

from ogtoberfest import preprocess, process_args, utils, files
from ogtoberfest import orthogroups_analyser as opa
from ogtoberfest import scorefuncs as sf

def preprocess_file(
        manager: process_args.Manager,
        file: pathlib.Path,
        funcs_dict: Dict[str, Callable], 
        method_func_name: Optional[str] = None,
    ):
    if method_func_name is not None:
        method_func = funcs_dict[method_func_name]
        if method_func_name in ["hieranoid", "broccoli"]:
            method_func(
                file,
                manager.options.output_path / file.name,
                manager.options.database_path,
            )
        else:
            method_func(
                file,
                manager.options.output_path / file.name,
            )


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


def general_stats(ref_ogs: Dict[str, Set[str]], 
                  pred_ogs: Dict[str, Set[str]]):

    V_prime = sf.V_raw_to_V_prime(ref_ogs, pred_ogs)
    V = sf.V_prime_to_V(V_prime)

    N = np.sum([len(refog) for refog in ref_ogs.values()])
    M = np.sum([len(predog) for predog in V.values()])
    M_raw = np.sum([len(predog) for predog in pred_ogs.values()])

    refog_species_dict = opa.get_refog_species_dict(ref_ogs)
    predog_species_dict = opa.get_predog_species_dict(V_prime)

    print("\nCalculating benchmarks:")
    print("%d  Number of RefOGs" % len(ref_ogs))
    print("%d  Number of genes in RefOGs" % N)
    print("%d  Number of PredOGs" % len(pred_ogs))
    print("%d  Number of genes in PredOGs (overlap)" % M)
    print("%d  Number of genes in PredOGs" % M_raw)
    print()

    general_stats_dict = {
        "nrefogs": len(ref_ogs),
        "n_genes_refogs": N,
        # "n_species_refogs": ,
        "npredogs": len(pred_ogs),
        # "n_species_predogs": 
        "n_overlaped_genes_predogs": M,
        "n_genes_predogs": M_raw
    }

    return general_stats_dict, refog_species_dict, predog_species_dict


def combine_scores(scores_list, score_names, avg_method = "mean"):

    score_dict = dict(zip(score_names, scores_list))

    combined_scores = 0.0
    for name, score in score_dict.items():
        if "recall" in name.lower() or "precision" in name.lower():
            combined_scores += score / 100
        elif "entropy" in name.lower():
            combined_scores += 1 - score
        else:
            combined_scores += 1 - score / 100

    combined_scores /= len(scores_list)
    return combined_scores


def get_scores(ref_ogs: Dict[str, Set[str]], 
               pred_ogs: Dict[str, Set[str]],  
               global_additional_scores: Optional[List[str]] = None,
               local_additional_scores: Optional[List[str]] = None,
               ) -> List[float]:

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

    total_num_species_refog = \
        len(set(utils.flatten_list_of_list([*refog_species_dict.values()])))

    total_num_species_predog_raw = len(
        set(utils.flatten_list_of_list([*predog_species_dict_raw.values()]))
    )

    total_num_species_predog_overlap = \
        len(
        set(utils.flatten_list_of_list([*predog_species_dict_overlap.values()]))
    )

    print("\nCalculating benchmarks:")
    print("%d  Number of RefOGs" % len(ref_ogs))
    print("%d  Number of species in RefOGs" % total_num_species_refog)
    print("%d  Number of genes in RefOGs" % N)
    print("%d  Number of PredOGs" % len(pred_ogs))
    print("%d  Number of species in PredOGs" % total_num_species_predog_raw)
    print("%d  Number of genes in PredOGs" % M_raw)
    print("%d  Number of species in PredOGs (overlap)" % total_num_species_predog_overlap)
    print("%d  Number of genes in PredOGs (overlap)" % M)

    print()

    # print(os.path.basename(ogs_filename))
    # gene_pair_TP, gene_pair_FP, gene_pair_FN, \
    # gene_pair_recall, gene_pair_precision, gene_pair_f1score = \
    # sf.calculate_benchmarks_pairwise(ref_ogs,
    #                                 pred_ogs,
    #                                 ref_ogs_uncertain)

    # print("%d  Correct gene pairs" % gene_pair_TP)
    # print("%d  False Positives gene pairs" % gene_pair_FP)
    # print("%d  False Negatives gene pairs\n" % gene_pair_FN)

    # print("%0.1f%%  Gene Pair Recall" % (100. * gene_pair_recall))
    # print("%0.1f%%  Gene Pair Precision" % (100. * gene_pair_precision))
    # print("%0.1f%%  Gene Pair F1-score\n" % (100. * gene_pair_f1score))

    missing_predogs_list = sf.missing_predogs(V_prime)
    missing_predogs = len(missing_predogs_list) / len(ref_ogs)
    print("%0.1f%%  Missing PrefOGs" % (100. * missing_predogs))

    total_missing_genes, total_missing_genes_set, \
    missing_genes_dict, missing_genes_count_dict, \
    missing_genes_proportion_dict = sf.check_missing_genes(ref_ogs, V_prime)
    missing_genes = total_missing_genes / N
    print("%0.1f%%  Missing Genes" % (100.0 * missing_genes))

    fussion_refog_set, fussion_predog_dict = sf.fussion(V_prime, V)
    fussion_refog_score = len(fussion_refog_set) / len(ref_ogs)
    if len(V) != 0:
        fussion_predog_score = 100. * len(fussion_predog_dict) / len(V)
    else:
        fussion_predog_score = 0.0
    print("%0.1f%%  Fussion (RefOG)" % (100.0 * fussion_refog_score))
    # print("%0.1f%%  Fussion (PredOG)" % fussion_predog_score)

    fission_refog_set, fission_predog_set = sf.fission(ref_ogs, V_prime)
    fission_refog_score = len(fission_refog_set) / len(ref_ogs)
    if len(V) != 0:
        fission_predog_score = 100. * len(fission_predog_set) / len(V)
    else:
        fission_predog_score = 0.0
    print("%0.1f%%  Fission (RefOG)" % (100.0 * fission_refog_score))
    # print("%0.1f%%  Fission (PredOG)" % fission_predog_score)
    # print()

    # micro_recall, micro_precision, \
    # micro_f1score, micro_fowlkes_mallows_index, \
    # micro_jaccard_index, micro_dissimilarity, micro_distance = sf.micro_scores(ref_ogs, V_prime, N)

    # print("%0.1f%%  micro Recall" % (100. * micro_recall))
    # print("%0.1f%%  micro Precision" % (100. * micro_precision))
    # print("%0.1f%%  micro F1-score" % (100. * micro_f1score))
    # print("%0.1f%%  micro Fowlkes Mallows Index" % (100. * micro_fowlkes_mallows_index))
    # print("%0.1f%%  micro Jaccard Index" % (100. * micro_jaccard_index))
    # print("%0.1f%%  micro Disimilarity" % (100. * micro_dissimilarity))
    # print("%0.2f   micro Distance" % (micro_distance))
    # print()

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

    # print("%0.1f%%  macro Recall" % (100. * macro_recall))
    # print("%0.1f%%  macro Precision" % (100. * macro_precision))
    # print("%0.1f%%  macro F1-score" % (100. * macro_f1score))
    # print("%0.1f%%  macro Fowlkes Mallows Index" % (100. * macro_fowlkes_mallows))
    # print("%0.1f%%  macro Jaccard Index" % (100. * macro_JI))
    # print("%0.1f%%  macro Disimilarity" % (100. * macro_dissimilarity))
    # print("%0.2f   macro Distance" % (macro_distance))
    # print("%0.2f   macro Disimilarity Distance" % (macro_dissimilarity_distance))
    # print()

    print("%0.1f%%  Weighted Avg Recall" % (100. * total_weighted_recall))
    print("%0.1f%%  Weighted Avg Precision" % (100. * total_weighted_precision))
    print("%0.1f%%  Weighted Avg F1-score" % (100. * total_weighted_f1score))
    # print("%0.1f%%  Weighted Avg Fowlkes Mallows Index" % (100. * total_weighted_fowlkes_mallows))
    # print("%0.1f%%  Weighted Avg Jaccard Index" % (100. * total_weighted_JI))
    # print("%0.1f%%  Weighted Avg Disimilarity" % (100. * total_weighted_dissimilarity))
    # print("%0.2f   Weighted Avg Distance" % (total_weighted_distance))
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

    num_decimal = 1
    scores = [
        np.round(100.0 * missing_predogs, num_decimal),
        np.round(100.0 * missing_genes, num_decimal),
        np.round(100.0 * fussion_refog_score, num_decimal),
        np.round(100.0 * fission_refog_score, num_decimal),
        np.round(100.0 * total_weighted_recall, num_decimal),
        np.round(100.0 * total_weighted_precision, num_decimal),
        np.round(total_entropy, 2),
    ]

    # if isinstance(global_additional_scores, list):
    return scores


def main(args: Optional[List[str]] = None):

    if not args:
        args = sys.argv[1:]

    args = utils.check_cmd_args(args)

    task = args.pop(0)

    manager = process_args.create_options(args, task=task)
    filewriter = files.FileWriter(manager.options.output_path)
    method_name_maps = utils.get_func_name_map()
    if task == "preprocess":
        funcs_dict = utils.get_func_name(preprocess)
        if manager.options.input_path.is_dir():
            for file in manager.options.input_path.iterdir():
                print(f"Preprocessing {file}")
                method = re.split("_|\.", file.name)[0]
                method = method.lower()
                method_func_name = method_name_maps.get(method)

                if method_func_name == "hieranoid" \
                    and manager.options.database_path is None:
                    print("Hieranoid needs to provide a database!")
                    continue
                preprocess_file(
                    manager,
                    file,
                    funcs_dict, 
                    method_func_name,
                )    

        elif manager.options.input_path.is_file():
            method = re.split("_|\.", manager.options.input_path.name)[0]
            method = method.lower()
            method_func_name = method_name_maps.get(method)
            if method_func_name == "hieranoid":
                if manager.options.database_path is None:
                    print("Hieranoid needs to provide a database!")
                    sys.exit(1)

            preprocess_file(
                manager,
                manager.options.input_path,
                funcs_dict,
                method_func_name,
            )

    elif task == "benchmark":

        global_score_colnames = [
            "Methods",
            "Missing PredOGs",
            "Missing Genes",
            "Fussion (RefOG)",
            "Fission (RefOG)",
            "Weighted Avg Recall",
            "Weighted Avg Precision",
            "Entropy",
            "Combined scores"
        ]

        global_additional_scores = None
        if hasattr(manager.options, "additional_global_scores"):
            global_additional_scores = manager.options.additional_global_scores

        ogreader = files.OGReader(manager.options.refog_path)
        if "OrthoBench" in manager.options.input_path.parent.name:
            print("\nReading RefOGs from: %s" % manager.options.refog_path)
            exp_genes = opa.get_expected_genes(
                    manager.options.database_path,
                    manager.options.outgroups,
                    manager.options.additional_species,
                )
            refogs_dict = ogreader.read_orthobench_refogs()
            predogs_dict = {}
            if manager.options.input_path.is_dir():
                global_scores_dict = {}
                for file in manager.options.input_path.iterdir():
                    print("\nReading predicted orthogroups from: %s" % file.name)
                    method = re.split("_|\.", file.name)[0]
                    method = method.lower()
                    method_func_name = method_name_maps.get(method)
                    predogs_dict[method_func_name] = ogreader.read_orthobench_predogs(
                        file
                    )

                    opa.check_orthobench_orthogroups(predogs_dict[method_func_name], exp_genes)

                    global_scores = get_scores(
                        refogs_dict, 
                        predogs_dict[method_func_name],
                        global_additional_scores,
                        )
                    combined_score = combine_scores(global_scores, global_score_colnames[1:-1])
                    
                    global_scores.append(np.round(combined_score, 3))
                    global_scores_dict[file.name.rsplit(".", 1)[0]] = global_scores

                    print("*" * 50)

                global_score_filename = manager.options.input_path.parent.name + "_global_scores.tsv"
                filewriter.save_global_scores(
                    global_score_colnames, global_score_filename, global_scores_dict
                )
            elif manager.options.input_path.is_file():
                print("\nReading predicted orthogroups from: %s" % manager.options.input_path.name)
                method = re.split("_|\.", manager.options.input_path.name)[0]
                method = method.lower()
                method_func_name = method_name_maps.get(method)
                predogs_dict[method_func_name] = \
                    ogreader.read_orthobench_predogs(manager.options.input_path)

    print()


if __name__ == "__main__":
    main()
