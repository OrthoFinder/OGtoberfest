import pathlib
from typing import Dict, List, Optional, Callable, Set, Tuple
import inspect

import numpy as np
from concurrent.futures import ThreadPoolExecutor, ProcessPoolExecutor, as_completed
from functools import partial
import traceback   

from ogtoberfest.utils import util
from ogtoberfest.run import process_args

from ogtoberfest.orthogroups import orthogroups_analyser as opa
from ogtoberfest.orthogroups import score_analyser as sa


runtime_dict = {
    "OF3_Linear": 38.525,
    "OF3_DB": 17.325,
    "OF3_Align_ST": 28.30833333,
    "OF2_DB": 16.40833333,
    "OF2_Align_ST": 100.5916667,
    "SP_def": 21.875,
    "SP_fast": 114.8583333,
    "SP_sens": 13.73333333,
    "Broccoli": 38.325,
    "ProteinOrtho": 5.466666667,
    "FastOMA": 106.4583333,
    "Hieranoid": 58.38333333,
    "OrthoMCL": 2000
}


def preprocess_file(
        file: pathlib.Path,
        method_func_name: str,
        method: str,
        speciesInfoObj: util.SpeciesInfo,
        manager: process_args.Manager,
        funcs_dict: Dict[str, Callable], 
        database: Optional[str] = None
    ):
    print(f"Processing {file.name}...")
    if method_func_name is not None:
        method_func = funcs_dict[method_func_name]
        if method is None:
            output_path = manager.options.output_path / file.name

            if method_func_name in ["hieranoid", "fastoma", "broccoli", "orthohmm"]:
                method_func(
                    file,
                    output_path,
                    manager.options.database_path,
                    speciesInfoObj.sequence2id_dict,
                    database,
                )
            else:
                method_func(
                    file,
                    output_path,
                    speciesInfoObj.sequence2id_dict,
                    database,
                )
        else:
            output_path = manager.options.output_path / method 
            output_path.mkdir(parents=True, exist_ok=True)
            output_path = output_path / f"{method}_all_genes_pairs_combined.txt"
            if method_func_name in ["hieranoid", "fastoma", "broccoli", "orthohmm"]:
                method_func(
                    file,
                    output_path,
                    manager.options.database_path,
                    speciesInfoObj.sequence2id_dict,
                )
            else:
                method_func(
                    file,
                    output_path,
                    speciesInfoObj.sequence2id_dict,
                )


def preprocess_files_parallel(
    file_method_list: List[Tuple[pathlib.Path, str, str]],
    speciesInfoObj,
    manager,
    funcs_dict,
    database=None,
    nthreads: int = 1,
):
    worker = partial(
        preprocess_file,
        speciesInfoObj=speciesInfoObj,
        manager=manager,
        funcs_dict=funcs_dict,
        database=database,
    )

    failures = {}       
    successes = 0

    with ThreadPoolExecutor(max_workers=nthreads) as pool:
        future_to_file = {pool.submit(worker, *t): t for t in file_method_list}

        for fut in as_completed(future_to_file):
            f = future_to_file[fut]
            filename = f[0].name 
            try:
                fut.result()
                successes += 1
            except Exception:
                failures[filename] = traceback.format_exc()

    if failures:
        print(f"✅ {successes} files processed successfully.")
        print(f"❌ {len(failures)} files failed:")
        for f, tb in failures.items():
            print(f"\n--- {f} ---\n{tb}")
        failed_files = ', '.join(str(p) for p in failures)
        raise RuntimeError(f"{len(failures)} files failed: {failed_files}")
    
    print(f"🎉 All {successes} files processed!!")

def compute_scores(
        file,
        predogs_dict,
        manager,
        filewriter,
        global_score_colnames,
        refogs_dict,
        refogs_ngenes_dict,
        metric,
        nthreads,
    ):  

    # print("\nProcessing orthogroups from: %s" % file.name)

    global_stats_dict, local_stats_dict, global_score_dict, \
        local_properties_dict, local_score_dict, other_info = sa.get_scores(
        refogs_dict,
        # predogs_dict[method_func_name],
        predogs_dict,
        manager.options.precision,
        nthreads,
        manager.options.use_id
    )

    effective_size_dist = sa.distribution_analyser(
        refogs_ngenes_dict, 
        local_properties_dict["Effective Size"], 
        metric
    )

    global_score_dict["Effective Size"] = np.round(effective_size_dist, manager.options.precision)
    method_runtime = runtime_dict.get(file.name.split(".")[0])
    global_score_dict["Runtime"] = np.round(method_runtime, manager.options.precision) if method_runtime is not None else np.nan

    global_scores = [
        global_score_dict.get(score_name)
        for score_name in manager.options.global_scores
    ]

    local_score_dict = {score: local_score_dict.get(score) for score in manager.options.local_scores}
    # local_scores_dict[file.name.split(".")[0]] = local_score_dict
    # complete_predog_dict[file.name.split(".")[0]] = other_info["Overlapped Genes"]
    complete_predog_dict = opa.add_missing_genes(
        other_info["Missing Genes"], 
        other_info["Overlapped Genes"]
    )

    local_filename = file.name.rsplit(".", 1)[0] + ".tsv"
## ---------------------------- save local scores ------------------------------------
    # local_score_filename = file.name.rsplit(".", 1)[0] + ".tsv"
    # filewriter.save_local_scores(
    #     manager.options.refog_stats_path, 
    #     local_score_filename, 
    #     local_score_dict,
    #     manager.options.precision,
    # )

    filewriter.save_local_info(
            local_filename, 
            local_score_dict,
            manager.options.precision,
            info_type="score",
    )

## --------------------------- Save local properties -------------------------------------
    filewriter.save_local_info(
            local_filename, 
            local_properties_dict,
            manager.options.precision,
            info_type="properties",
    )

## --------------------- Average scores ---------------------------------
    if manager.options.combined_global_score in ["Avg Score", "RMS Score"]:
        combined_score = sa.combine_scores(
            global_scores, 
            global_score_colnames[1:-1], 
            precision=manager.options.precision,
            avg_method=manager.options.combined_global_score,
            weights=manager.options.global_scor_weights
        )
        global_scores.append(np.round(combined_score, 3))

    # global_scores_dict[file.name.rsplit(".", 1)[0]] = global_scores
    if nthreads == 1:
        print("*" * 50)
## ----------------------------- Save missing genes ---------------------------------

    ## Save missing genes
    # missing_genes_filename = file.name.rsplit(".", 1)[0] + ".tsv"
    filewriter.save_missing_genes(
        local_filename,
        other_info["Missing Species"],
        other_info["Missing Genes"]
    )

##  ---------------------------- Save Fusion and Fission genes -------------------------------
    filewriter.save_fusion_genes(
        local_filename, 
        other_info["Fusion Genes"]
    )
    filewriter.save_fission_genes(
        local_filename, 
        other_info["Fission Genes"],
    )

## ----------------------------- Save overlap PredOGs ------------------------------------
    filewriter.save_overlap_predogs(
        local_filename, 
        other_info["Overlap PredOGs"]
    )

    return global_scores, local_score_dict, complete_predog_dict

def compute_scores_parallel(
    file_method_list,
    manager,
    filewriter,
    global_score_colnames,
    refogs_dict,
    refogs_ngenes_dict,
    metric,
    nthreads: int = 1,
):
    worker = partial(
        compute_scores,
        manager=manager,
        filewriter=filewriter,
        global_score_colnames=global_score_colnames,
        refogs_dict=refogs_dict,
        refogs_ngenes_dict=refogs_ngenes_dict,
        metric=metric,
        nthreads=nthreads,
    )

    global_scores_dict   = {}
    local_scores_dict    = {}
    complete_predog_dict = {}
    failures, successes  = {}, 0

    with ProcessPoolExecutor(max_workers=nthreads) as pool:
        future_to_file = {pool.submit(worker, *t): t for t in file_method_list}

        for fut in as_completed(future_to_file):
            f = future_to_file[fut]
            filename = f[0].name 
            print("Done processing orthogroups from: %s" % filename)

            method = filename.split(".")[0]
            g_scores, l_scores, comp_pred = fut.result()
            
            global_scores_dict[method] = g_scores
            local_scores_dict[method] = l_scores
            complete_predog_dict[method] = comp_pred
            try:
                fut.result()
                successes += 1
            except Exception:
                failures[filename] = traceback.format_exc()

    if failures:
        print(f"✅ {successes} succeeded.")
        print(f"❌ {len(failures)} failed:")
        for f, tb in failures.items():
            print(f"\n--- {f} ---\n{tb}")
        raise RuntimeError(f"{len(failures)} jobs failed")

    print(f"🎉 All {successes} files processed!!")
    return global_scores_dict, local_scores_dict, complete_predog_dict