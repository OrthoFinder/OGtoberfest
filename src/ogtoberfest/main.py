import sys
from typing import Dict, List, Optional, Callable, Set
import re
import pathlib
import numpy as np
import pandas as pd


from ogtoberfest import preprocess, process_args, utils, files
from ogtoberfest import orthogroups_analyser as opa
from ogtoberfest import scorefuncs as sf
from ogtoberfest import score_analyser as sa

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

        if len(manager.options.additional_global_scores) != 0:
            if isinstance(manager.options.additional_global_scores, str):
                manager.options.global_scores.append(manager.options.additional_global_scores)
            elif isinstance(manager.options.additional_global_scores, list):
                manager.options.global_scores.extend(manager.options.additional_global_scores)
        
        global_score_colnames = ["Methods"] \
            + manager.options.global_scores \
            + [manager.options.combined_global_score] 

        if "OrthoBench" in manager.options.input_path.parent.name:
            print("\nReading RefOGs from: %s" % manager.options.refog_path)
            exp_genes = opa.get_expected_genes(
                    manager.options.database_path,
                    manager.options.outgroups,
                    manager.options.additional_species,
                )
            ogreader = files.OrthoBenchOGReader(manager.options.refog_path, 
                                      manager.options.uncertian_refog_path)
            refogs_dict = ogreader.read_orthobench_refogs()
            # uncertain_refogs_dict = ogreader.read_uncertain_orthobench_refogs()
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

                    global_stats_dict, local_stats_dict, global_score_dict, local_score_dict, predogs_info = sa.get_scores(
                        refogs_dict,
                        predogs_dict[method_func_name],
                        manager.options.precision,
                    )
                    
                    global_scores = [
                        score for score_name, score in global_score_dict.items()
                        if score_name in manager.options.global_scores
                    ]

                    if manager.options.combined_global_score == "Simple Avg Score":
                        combined_score = sa.combine_scores(global_scores, 
                                                           global_score_colnames[1:-1], 
                                                           precision=manager.options.precision,)
                        global_scores.append(np.round(combined_score, 3))

                    global_scores_dict[file.name.rsplit(".", 1)[0]] = global_scores

                    print("*" * 50)
                if manager.options.combined_global_score == "Avg Rank Score":
                    ranked_global_scores_dict, global_scores_rank_dict = sa.rand_score(
                        global_scores_dict, 
                        global_score_colnames[1:-1], 
                        precision=manager.options.precision, 
                        rank_method=manager.options.rank_method
                    )

                global_score_filename = manager.options.input_path.parent.name + "_global_scores.tsv"
                filewriter.save_global_scores(
                    global_score_colnames, global_score_filename, ranked_global_scores_dict
                )

                global_score_rank_filename = (
                    manager.options.input_path.parent.name + "_global_scores_rank.tsv"
                )
                filewriter.save_global_scores(
                    global_score_colnames[:-1], global_score_rank_filename, global_scores_rank_dict
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
