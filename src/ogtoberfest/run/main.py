import sys
from typing import Dict, List, Optional, Callable, Set
import re
import pathlib
import asyncio
import inspect

import numpy as np
import pandas as pd

from ogtoberfest.orthogroups import orthogroups_preprocess, orthogroups_files
from ogtoberfest.utils import util, files
from ogtoberfest.run import process_args
from ogtoberfest.orthologues import orthologues_preprocess
from ogtoberfest.orthogroups import orthogroups_analyser as opa
from ogtoberfest.orthogroups import scorefuncs as sf
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
        speciesInfoObj: util.SpeciesInfo,
        manager: process_args.Manager,
        file: pathlib.Path,
        funcs_dict: Dict[str, Callable], 
        method_func_name: Optional[str] = None,
        method: Optional[str] = None,
        database: Optional[str] = None
    ):
    if method_func_name is not None:
        method_func = funcs_dict[method_func_name]
        if method is None:
            output_path = manager.options.output_path / file.name
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
                database,
            )
        else:
            method_func(
                file,
                output_path,
                speciesInfoObj.sequence2id_dict,
                database,
            )

def compute_scores(
        file,
        manager,
        method_func_name,
        refogs_dict,
        predogs_dict,
        local_scores_dict,
        complete_predog_dict,
        filewriter,
        global_score_colnames,
        global_scores_dict,
        refogs_ngenes_dict,
        metric,
):  

    global_stats_dict, local_stats_dict, global_score_dict, \
        local_properties_dict, local_score_dict, other_info = sa.get_scores(
        refogs_dict,
        predogs_dict[method_func_name],
        manager.options.precision,
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
    local_scores_dict[file.name.split(".")[0]] = local_score_dict
    # complete_predog_dict[file.name.split(".")[0]] = other_info["Overlapped Genes"]
    complete_predog_dict[file.name.split(".")[0]] = opa.add_missing_genes(
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

    global_scores_dict[file.name.rsplit(".", 1)[0]] = global_scores

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

    return global_scores_dict, local_scores_dict, other_info, complete_predog_dict

def main(args: Optional[List[str]] = None):

    if not args:
        args = sys.argv[1:]

    args = util.check_cmd_args(args)

    main_task = args.pop(0)

    manager = process_args.create_options(args, task=main_task)
    method_name_maps = util.get_func_name_map()
    
    filehandler = files.FileHandler(
        manager.options.database_path,
        manager.options.wd_base,
    )
    
    speciesInfoObj = filehandler.process_fasta_files()
    input_path_isfile = False
    method_file_list = []
    if manager.options.input_path.is_dir():
        for dir_file in manager.options.input_path.iterdir():
            method = re.split("_|\.", dir_file.name)[0]
            method = method.lower()
            method_func_name = method_name_maps.get(method)
        
            if dir_file.is_file():
                method_file_list.append((method, method_func_name, dir_file))

            else:
                for file in dir_file.iterdir():
                    method_file_list.append((method, method_func_name, file))

    else:
        input_path_isfile = True
        method = re.split("_|\.", manager.options.input_path.name)[0]
        method = method.lower()
        method_func_name = method_name_maps.get(method)
        method_file_list.append((method, method_func_name, manager.options.input_path))
    
    database = "OrthoBench" if "orthobench" == manager.options.input_path.parent.name.lower() else None
    use_id = False if database == "OrthoBench" else manager.options.use_id
    ## ----------------------------- Preprocessing ---------------------------
    if "preprocess" in  main_task:
        if main_task == "orthogroups_preprocess":
            funcs_dict = util.get_func_name(orthogroups_preprocess)
        elif main_task == "orthologues_preprocess":
            funcs_dict = util.get_func_name(orthologues_preprocess)

        if manager.options.input_path.is_dir():
            for dir_file in manager.options.input_path.iterdir():
                method = re.split("_|\.", dir_file.name)[0]
                method = method.lower()
                method_func_name = method_name_maps.get(method)

                if method_func_name in ["hieranoid", "fastoma", "broccoli", "orthohmm"] \
                    and manager.options.database_path is None:
                    print(f"{method_func_name.titile()} needs to provide a database!")
                    continue

                if main_task == "orthogroups_preprocess":
                    method = None

                if dir_file.is_file():
                    print(f"Preprocessing {dir_file}")

                    preprocess_file(
                        speciesInfoObj,
                        manager,
                        dir_file,
                        funcs_dict, 
                        method_func_name,
                        method,
                        database=database
                    ) 
                else:
                    for file in dir_file.iterdir():
                        print(f"Preprocessing {file}")
                        
                        preprocess_file(
                            speciesInfoObj,
                            manager,
                            file,
                            funcs_dict, 
                            method_func_name,
                            method,
                            database=database
                        ) 


        elif manager.options.input_path.is_file():
            method = re.split("_|\.", manager.options.input_path.name)[0]
            method = method.lower()
            method_func_name = method_name_maps.get(method)
            if method_func_name in ["hieranoid", "fastoma", "broccoli"]:
                if manager.options.database_path is None:
                    print(f"{method_func_name.titile()} needs to provide a database!")
                    sys.exit(1)

            preprocess_file(
                speciesInfoObj,
                manager,
                manager.options.input_path,
                funcs_dict,
                method_func_name,
                database=database
            )
        print(f"All the files have been preprocessed!")
        sys.exit(0)
    
    ## ---------------------- Benchmarking --------------------------

    if main_task == "orthogroups_benchmark":
        # filehandler = files.FileHandler(manager.options.output_path)
        filewriter = orthogroups_files.FileWriter(
            manager.options.output_path,
            input_path_isfile,
            speciesInfoObj.id2sequence_dict,
            use_id
        )
        # filewriter = files.FileWriter(
        #     filehandler.global_score_path,
        #     filehandler.local_score_path, 
        #     filehandler.other_info_path,
        #     filehandler.vi_score_path,
        #     filehandler.dist_path,
        #     )

        if manager.options.additional_global_scores is not None:
            if isinstance(manager.options.additional_global_scores, str):
                manager.options.global_scores.append(manager.options.additional_global_scores)
            elif isinstance(manager.options.additional_global_scores, list):
                manager.options.global_scores.extend(manager.options.additional_global_scores)
        
        global_score_colnames = ["Methods"] \
            + manager.options.global_scores \
            + [manager.options.combined_global_score] 

        if manager.options.additional_local_scores is not None:
            if isinstance(manager.options.additional_local_scores, str):
                manager.options.local_scores.append(manager.options.additional_local_scores)
            elif isinstance(manager.options.additional_global_scores, list):
                manager.options.local_scores.extend(manager.options.additional_local_scores)
        
        print("\nReading RefOGs from: %s" % manager.options.refog_path)

        if "OrthoBench" in manager.options.input_path.parent.name:
            exp_genes = opa.get_expected_genes(
                    manager.options.database_path,
                    manager.options.outgroups,
                    manager.options.additional_species,
                )
            ogreader = orthogroups_files.OGReader(
                manager.options.refog_path, 
                manager.options.uncertian_refog_path
            )
            refogs_dict = ogreader.read_orthobench_refogs()
            refogs_nspecies_dict, refogs_ngenes_dict = \
                ogreader.read_ogs_stats(manager.options.refog_stats_path)
            # uncertain_refogs_dict = ogreader.read_uncertain_orthobench_refogs()
            predogs_dict = {}
            # if manager.options.input_path.is_dir():
            global_scores_dict = {}
            local_scores_dict = {}
            complete_predog_dict = {}
            
            for method, method_func_name, file in method_file_list:

            # for file in manager.options.input_path.iterdir():
                print("\nReading predicted orthogroups from: %s" % file.name)
                # method = re.split("_|\.", file.name)[0]
                # method = method.lower()
                # method_func_name = method_name_maps.get(method)
                predogs_dict[method_func_name] = ogreader.read_orthobench_predogs(
                    file
                )

                opa.check_orthobench_orthogroups(predogs_dict[method_func_name], exp_genes)

                global_scores_dict, local_scores_dict, other_info, complete_predog_dict = \
                    compute_scores(
                        file,
                        manager,
                        method_func_name,
                        refogs_dict,
                        predogs_dict,
                        local_scores_dict,
                        complete_predog_dict,
                        filewriter,
                        global_score_colnames,
                        global_scores_dict,
                        refogs_ngenes_dict,
                        manager.options.metric,
                    )
            
            # elif manager.options.input_path.is_file():
            #     print("\nReading predicted orthogroups from: %s" % manager.options.input_path.name)
            #     method = re.split("_|\.", manager.options.input_path.name)[0]
            #     method = method.lower()
            #     method_func_name = method_name_maps.get(method)
            #     predogs_dict[method_func_name] = \
            #         ogreader.read_orthobench_predogs(manager.options.input_path)
                
        else: #if "sim" in manager.options.input_path.parent.name.lower():
            ogreader = orthogroups_files.OGReader(
                manager.options.refog_path, 
                manager.options.uncertian_refog_path,
                speciesInfoObj.sequence2id_dict
            )
            refogs_dict = ogreader.read_ogs(manager.options.refog_path, use_id=use_id)
            refogs_nspecies_dict, refogs_ngenes_dict = \
                ogreader.read_ogs_stats(manager.options.refog_stats_path)

            predogs_dict = {}
            # if manager.options.input_path.is_dir():
            global_scores_dict = {}
            local_scores_dict = {}
            complete_predog_dict = {}
            for method, method_func_name, file in method_file_list:
            # for file in manager.options.input_path.iterdir():
                print("\nReading predicted orthogroups from: %s" % file.name)
                method = re.split("_|\.", file.name)[0]
                method = method.lower()
                method_func_name = method_name_maps.get(method)
                predogs_dict[method_func_name] = ogreader.read_ogs(file, use_id=False)
                global_scores_dict, local_scores_dict, other_info, complete_predog_dict = \
                    compute_scores(
                        file,
                        manager,
                        method_func_name,
                        refogs_dict,
                        predogs_dict,
                        local_scores_dict,
                        complete_predog_dict,
                        filewriter,
                        global_score_colnames,
                        global_scores_dict,
                        refogs_ngenes_dict,
                        manager.options.metric,
                    )
            

        ## =============================== Score analysis =============================
        if manager.options.input_path.is_dir():
        ## ------------------------- Correlation analysis -----------------------------
            sa.corr_vi_analysis(
                complete_predog_dict, 
                refogs_ngenes_dict,
                filewriter,
                precision=manager.options.precision,
            )

            sa.corr_dist_analysis(
                local_scores_dict,
                manager.options.local_scores,
                filewriter,
                metric=manager.options.metric,
                precision=manager.options.precision,
            )
    
        ## --------------------- Rank score -----------------------------------
            if manager.options.combined_global_score == "Rank Score":

                global_scores_dict, global_scores_rank_dict = sa.rank_score(
                    global_scores_dict, 
                    global_score_colnames[1:-1], 
                    precision=manager.options.precision, 
                    rank_method=manager.options.rank_method
                )

                global_score_rank_filename = (
                    manager.options.input_path.parent.name + "_global_scores_rank.tsv"
                )
                filewriter.save_global_scores(
                    global_score_colnames[:-1], global_score_rank_filename, global_scores_rank_dict
                )
        ##  --------------------------- Z-score ---------------------------------- 
            elif manager.options.combined_global_score == "Z-score":
                    global_scores_dict = sa.z_score(
                        global_scores_dict, 
                        global_score_colnames[1:-1], 
                        precision=manager.options.precision, 
                    )

            global_score_filename = manager.options.input_path.parent.name + "_global_scores.tsv"
            filewriter.save_global_scores(
                global_score_colnames, global_score_filename, global_scores_dict
            )

    elif main_task == "orthologues_benchmark":
        pass

    print()


if __name__ == "__main__":
    main()

