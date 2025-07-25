import sys
from typing import Dict, List, Optional, Callable, Set
import re
import pathlib
import inspect
from ogtoberfest.orthogroups import orthogroups_preprocess, orthogroups_files
from ogtoberfest.orthogroups import parallel_processing as orthogroups_pp
from ogtoberfest.utils import util, files
from ogtoberfest.run import process_args

from ogtoberfest.orthologues import orthologues_preprocess
from ogtoberfest.orthogroups import orthogroups_analyser as opa
from ogtoberfest.orthogroups import score_analyser as sa


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
            
        file_method_list = []
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
                    file_method_list.append((dir_file, method_func_name, method))
                else:
                    for file in dir_file.iterdir():
                        file_method_list.append((file, method_func_name, method))

        elif manager.options.input_path.is_file():
            method = re.split("_|\.", manager.options.input_path.name)[0]
            method = method.lower()
            method_func_name = method_name_maps.get(method)
            if method_func_name in ["hieranoid", "fastoma", "broccoli"]:
                if manager.options.database_path is None:
                    print(f"{method_func_name.titile()} needs to provide a database!")
                    sys.exit(1)
            file_method_list.append((manager.options.input_path, method_func_name, method))

        orthogroups_pp.preprocess_files_parallel(
            file_method_list,
            speciesInfoObj,
            manager,
            funcs_dict, 
            database=database,
            nthreads=manager.options.nthreads
        )         
        sys.exit(0)
    
    ## ---------------------- Benchmarking --------------------------

    if main_task == "orthogroups_benchmark":
        # filehandler = files.FileHandler(manager.options.output_path)
        filewriter = orthogroups_files.FileWriter(
            manager.options.output_path,
            input_path_isfile,
            speciesInfoObj.id2sequence_dict,
            speciesInfoObj.id2species_dict,
            use_id
        )

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
        
        print("\nReading RefOGs from: %s" % manager.options.refog_path, end=2*"\n")

        file_method_list = []
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
            
            for method, method_func_name, file in method_file_list:

                predogs_dict = ogreader.read_orthobench_predogs(file)
                print(f"*** Check orthogroups for {method} ***")
                opa.check_orthobench_orthogroups(predogs_dict, exp_genes)
                print()

                file_method_list.append((
                    file,
                    predogs_dict
                ))    
                
        else: #if "sim" in manager.options.input_path.parent.name.lower():
            ogreader = orthogroups_files.OGReader(
                manager.options.refog_path, 
                manager.options.uncertian_refog_path,
                speciesInfoObj.sequence2id_dict
            )
            refogs_dict = ogreader.read_ogs(manager.options.refog_path, use_id=use_id)
            refogs_nspecies_dict, refogs_ngenes_dict = \
                ogreader.read_ogs_stats(manager.options.refog_stats_path)

            for method, method_func_name, file in method_file_list:
            # for file in manager.options.input_path.iterdir():
                # print("\nReading predicted orthogroups from: %s" % file.name)
                method = re.split("_|\.", file.name)[0]
                method = method.lower()
                method_func_name = method_name_maps.get(method)
                predogs_dict = ogreader.read_ogs(file, use_id=False)
                file_method_list.append((
                    file,
                    predogs_dict,
                ))  
                  
        global_scores_dict, local_scores_dict, complete_predog_dict = \
            orthogroups_pp.compute_scores_parallel(
            file_method_list,
            manager,
            filewriter,
            global_score_colnames,
            refogs_dict,
            refogs_ngenes_dict,
            manager.options.metric,
            manager.options.nthreads,
        )

        
        ## =============================== Score analysis =============================
        if manager.options.input_path.is_dir():
        ## ------------------------- Correlation analysis -----------------------------
            # sa.corr_vi_analysis(
            #     complete_predog_dict, 
            #     refogs_ngenes_dict,
            #     filewriter,
            #     precision=manager.options.precision,
            # )

            # sa.corr_dist_analysis(
            #     local_scores_dict,
            #     manager.options.local_scores,
            #     filewriter,
            #     metric=manager.options.metric,
            #     precision=manager.options.precision,
            # )
    
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

