import pathlib
import os 
import re
import sys
from typing import Optional, List, Dict, Any, Set
import numpy as np

class OGReader:

    def __init__(self, 
                 refog_path, 
                 uncerain_refog_path: Optional[pathlib.Path]=None
    ):
        self.refog_path = refog_path
        self.uncertain_refog_path = uncerain_refog_path


    def read_ogs(self, og_path)-> Dict[str, Set[str]]:
        refogs = {}
        try:
            with open(og_path, "r") as reader:
                for line in reader:
                    line = line.strip().split(": ")
                    refog_key, genes = line[0], line[-1]
                    genes = [g.strip() for g in genes.split(", ") if g.strip() != ""]
                    refogs[refog_key] = set(genes)

        except:
            print("ERROR: RefOG file not found in: %s" % og_path)
            sys.exit(1)
        
        return refogs
    
    def read_orthobench_refogs(self) -> Dict[str, Set[str]]:
        n_expected = 1945
        refogs = self.read_ogs(self.refog_path)
        n_genes = np.sum([len(refog) for refog in refogs.values()])

        if "uncertain" in self.refog_path.name:
            assert (
                n_expected == n_genes
            ), "ERROR: There are genes missing from the RefOG files. Found %d, expected %d" % (
                n_genes,
                n_expected,
            )

        return refogs

    # def read_orthobench_refogs(self) -> Dict[str, Set[str]]:
    #     n_expected = 1945
    #     n_genes = 0
    #     refogs = {}
    #     refogs_ngenes_dict = {}
    #     refogs_nspecies_dict = {}
    #     try:
    #         with open(self.refog_path, "r") as reader:
    #             for line in reader:
    #                 if "Orthogroups" in line:
    #                     continue
    #                 line = line.strip().split("\t")
    #                 refog_key, nspecies, ngenes, genes = line[0], line[1], line[2], line[-1]
    #                 genes = [g.strip() for g in genes.split(", ") if g.strip() != ""]
    #                 refogs[refog_key] = set(genes)
    #                 refogs_ngenes_dict[refog_key] = int(ngenes)
    #                 refogs_nspecies_dict[refog_key] = int(nspecies)
    #                 n_genes += len(refogs[refog_key])

    #     except:
    #         print("ERROR: RefOG file not found in: %s" % self.refog_path)
    #         sys.exit(1)

    #     if "uncertain" in self.refog_path.name:
    #         assert (
    #             n_expected == n_genes
    #         ), "ERROR: There are genes missing from the RefOG files. Found %d, expected %d" % (
    #             n_genes,
    #             n_expected,
    #         )

    #     return refogs, refogs_nspecies_dict, refogs_ngenes_dict
    
    
    def read_ogs_stats(self, og_stats_path) -> Dict[str, Set[str]]:
        refogs_ngenes_dict = {}
        refogs_nspecies_dict = {}
        try:
            with open(og_stats_path , "r") as reader:
                for n, line in enumerate(reader):
                    line = line.strip().split("\t")
                    if n == 0:
                        colnames = {
                            col: i for i, col in enumerate(line)
                        }
                    else:
                        nspecies_index = colnames["num_species"]
                        ngenes_index = colnames["num_genes"]
                        refog_key, nspecies, ngenes = line[0], line[nspecies_index], line[ngenes_index],
                        refogs_ngenes_dict[refog_key] = int(ngenes)
                        refogs_nspecies_dict[refog_key] = int(nspecies)

        except:
            print("ERROR: RefOG stats file not found in: %s" % og_stats_path)
            sys.exit(1)

        return refogs_nspecies_dict, refogs_ngenes_dict

    def read_uncertain_orthobench_refogs(self) -> Dict[str, Set[str]]:
        uncertain_refogs = {}

        with open(self.uncertain_refog_path, "r") as reader:
            for line in reader:
                if "Orthogroups" in line:
                    continue
                line = line.strip().split("\t")
                refog_key, genes = line[0], line[-1]
                genes = [g.strip() for g in genes.split(", ") if g.strip() != ""]
                uncertain_refogs[refog_key] = set(genes)

        return uncertain_refogs

    def read_orthobench_predogs(self, predog_path):
        ogs = {}
        # "WBGene00\d+\.1|ENSCAFP\d+|ENSCINP\d+|ENSDARP\d+|FBpp0\d+|ENSGALP\d+|ENSP000\d+|ENSMODP\d+|ENSMUSP\d+|ENSPTRP\d+|ENSRNOP\d+|ENSTNIP\d+"
        gene_pattern_list = [
            "Caenorhabditis_elegans.WBGene00\d+\.1",
            "Canis_familiaris.ENSCAFP\d+",
            "Ciona_intestinalis.ENSCINP\d+",
            "Danio_rerio.ENSDARP\d+",
            "Drosophila_melanogaster.FBpp0\d+",
            "Gallus_gallus.ENSGALP\d+",
            "Homo_sapiens.ENSP000\d+",
            "Monodelphis_domestica.ENSMODP\d+",
            "Mus_musculus.ENSMUSP\d+",
            "Pan_troglodytes.ENSPTRP\d+",
            "Rattus_norvegicus.ENSRNOP\d+",
            "Tetraodon_nigroviridis.ENSTNIP\d+",
        ]

        gene_pat = re.compile("|".join(gene_pattern_list))
        with open(predog_path, "r") as infile:
            for i, l in enumerate(infile):
                genes = re.findall(gene_pat, l)
                if len(genes) > 0:
                    if ":" in l:
                        og_label = l.split(":")[0]
                    else:
                        og_label = l.split(None, 1)[0]
                    ogs[og_label] = set(genes)
        return ogs
    

    def read_predogs(self, predog_path):
        ogs = {}
        with open(predog_path, "r") as infile:
            for line in infile:
                line = line.strip().split(": ")
                og_label, genes = line[0].replace(":", ""), line[1]
                ogs[og_label] = set(genes.split(", "))
        return ogs

class FileHandler:
    def __init__(
            self, 
            output_path: pathlib.Path, 
            input_path_isfile: bool = False
    ):
        self.input_path_isfile = input_path_isfile
        self.global_score_path = output_path
        if self.input_path_isfile:
            self.local_score_path = self.global_score_path
            self.local_property_path = self.global_score_path
            self.missing_genes_path = self.global_score_path
            self.fusion_path = self.global_score_path
            self.fission_path = self.global_score_path

        else:
            self.local_score_path = self.global_score_path / "local_scores"
            self.local_property_path = self.global_score_path / "local_properties"
            self.missing_genes_path = self.global_score_path / "missing_genes"
            self.fusion_path = self.global_score_path / "fusion"
            self.fission_path = self.global_score_path / "fission"
            self.vi_score_path = self.global_score_path / "local_vi_scores"
            self.dist_path = self.global_score_path / "local_score_distance"

             
            if not os.path.exists(self.local_score_path):
                os.makedirs(self.local_score_path, exist_ok=True)

            if not os.path.exists(self.local_property_path):
                os.makedirs(self.local_property_path, exist_ok=True)
            
            if not os.path.exists(self.missing_genes_path):
                os.makedirs(self.missing_genes_path, exist_ok=True)

            if not os.path.exists(self.fusion_path):
                os.makedirs(self.fusion_path, exist_ok=True)

            if not os.path.exists(self.fission_path):
                os.makedirs(self.fission_path, exist_ok=True)
            
            if not os.path.exists(self.vi_score_path):
                os.makedirs(self.vi_score_path, exist_ok=True)
            
            if not os.path.exists(self.dist_path):
                os.makedirs(self.dist_path, exist_ok=True)

    def get_method_local_score_fn(self, method: str):

        if self.input_path_isfile:
            return self.self.local_score_path
        
        for file in self.self.local_score_path.iterdir():
            if method in file.name:
                return file 
        
    def get_method_ortho_info_fn(self, method: str):
        if self.input_path_isfile:
            return self.self.other_info_path
        for file in self.self.other_info_path.iterdir():
            if method in file.name:
                return file         

class FileWriter(FileHandler):

    # def __init__(self, 
    #              global_score_path: pathlib.Path, 
    #              local_score_path: pathlib.Path,
    #              other_info_path: pathlib.Path,
    #              vi_score_path: pathlib.Path,
    #              dist_path: pathlib.Path,
    #              ):
    #     self.global_score_path = global_score_path
    #     self.local_score_path = local_score_path
    #     self.other_info_path = other_info_path
    #     self.vi_score_path = vi_score_path
    #     self.dist_path = dist_path

    def __init__(self, output_path: pathlib.Path, input_path_isfile: bool = False):
        super().__init__(output_path, input_path_isfile)

    def save_global_scores(self, 
                           colnames: List[str], 
                           filename: str, 
                           scores_dict: Dict[str, Any],):
        colnames_str = "\t".join(colnames)
        score_filename = self.global_score_path / filename 
        with open(score_filename, "w") as writer:
            writer.write(colnames_str + "\n")
            for method, scores in scores_dict.items():
                scores_str_list = [str(s) for s in scores]
                scores_str = "\t".join(scores_str_list)
                line = "\t".join((method, scores_str))
                writer.write(line + "\n")

    def save_local_info(self, 
                        filename: str, 
                        info_dict: Dict[str, Any],
                        precision: int,
                        info_type: str = "score",
                        ):
        
        colnames = ["RefOGs"]
        local_score_names, local_scores = [*zip(*info_dict.items())] 
        colnames.extend(local_score_names) 
        colnames_str = "\t".join(colnames)
        
        if info_type == "score":
            filename = self.local_score_path / filename 
        else:
            filename = self.local_property_path / filename

        refogs_dict = {}
        for local_score in local_scores:
            if local_score is None:
                continue
            for refog_key, info in local_score.items():
                if refog_key not in refogs_dict:
                    refogs_dict[refog_key] = [np.round(info, precision)]
                elif refog_key in refogs_dict:
                    refogs_dict[refog_key].append(np.round(info, precision))

        info_lists = [
            [refog_key] + info
            for refog_key, info in refogs_dict.items()
        ]

        info_lists = sorted(info_lists, key=lambda x: (x[1], x[2], x[3]), reverse=True)
        with open(filename, "w") as writer:
            writer.write(colnames_str + "\n")
            for info_list in info_lists:
                info_str_list = [str(s) for s in info_list]
                info_str = "\t".join(info_str_list) 
                writer.write(info_str + "\n")



    # def save_local_scores(self, 
    #                       refog_stats_path: str,
    #                       filename: str, 
    #                       scores_dict: Dict[str, Any],
    #                       precision: int,
    #                       ):
    #     refogs_dict = {}
    #     colnames = []
    #     with open(refog_stats_path, "r") as reader:
    #         for i, line in enumerate(reader):
    #             line = line.strip().split("\t")
    #             if i == 0:
    #                 colnames = line
    #             else:
    #                 refog_key = line[0]
    #                 refogs_dict[refog_key] = line[1:]

    #     local_score_names, local_scores = [*zip(*scores_dict.items())] 
    #     colnames.extend(local_score_names) 
    #     colnames_str = "\t".join(colnames)

    #     score_filename = self.local_score_path / filename 
     
    #     for local_score in local_scores:
    #         if local_score is None:
    #             continue
    #         for refog_key, scores in local_score.items():
    #             refogs_dict[refog_key].append(np.round(scores, precision))

    #     scores_lists = [
    #         [refog_key] + scores
    #         for refog_key, scores in refogs_dict.items()
    #     ]

    #     scores_lists = sorted(scores_lists, key=lambda x: x[3])
    #     with open(score_filename, "w") as writer:
    #         writer.write(colnames_str + "\n")
    #         for score_list in scores_lists:
    #             scores_str_list = [str(s) for s in score_list]
    #             scores_str = "\t".join(scores_str_list) 
    #             writer.write(scores_str + "\n")


    def save_missing_genes(self, 
                           filename: str, 
                           missing_species_dict: Dict[str, Set[str]],
                           missing_genes_dict: Dict[str, Set[str]]):
        
        colnames = ["RefOGs", "Missing Species", "Missing Genes"]
        colnames_str = "\t".join(colnames)
        other_info_filepath = self.missing_genes_path / filename
        with open(other_info_filepath, "w") as writer:
            writer.write(colnames_str + "\n")
            for refog_key, missing_genes in missing_genes_dict.items():
                missing_species = missing_species_dict[refog_key]
                missing_species_str = ", ".join(missing_species)
                missing_genes_str = ", ".join(missing_genes)
                writer.write("\t".join((refog_key, missing_species_str, missing_genes_str)) + "\n")

    def save_fusion_genes(self, 
                           filename: str, 
                           fusion_genes_dict: Dict[str, Set[str]]):
        
        colnames = ["RefOGs", "PredOGs", "Fused Genes"]
        colnames_str = "\t".join(colnames)
        fused_genes_filepath = self.fusion_path / filename
        with open(fused_genes_filepath, "w") as writer:
            writer.write(colnames_str + "\n")
            for refog_key, tp_dict in fusion_genes_dict.items():
                if len(tp_dict) != 0:
                    for predog_key, tp in tp_dict.items():
                        line = "\t".join((refog_key, predog_key, ", ".join(tp)))
                        writer.write(line + "\n")

    def save_fission_genes(self, 
                           filename: str, 
                           fission_genes_dict: Dict[str, Set[str]]):
        
        colnames = ["RefOGs", "PredOGs", "Fission Genes"]
        colnames_str = "\t".join(colnames)
        fission_genes_filepath = self.fission_path / filename
        with open(fission_genes_filepath, "w") as writer:
            writer.write(colnames_str + "\n")
            for refog_key, tp_dict in fission_genes_dict.items():
                if len(tp_dict) != 0:
                    for predog_key, tp in tp_dict.items():
                        line = "\t".join((refog_key, predog_key, ", ".join(tp)))
                        writer.write(line + "\n")

    def save_global_VI_scores(
            self, 
            global_VI_dict, 
    ):

        global_VI_filepath = self.dist_path / "avg_VI.tsv"
        with open(global_VI_filepath, "w") as writer:
            for i, (m1, VI_dict) in enumerate(global_VI_dict.items()):
                if i == 0:
                    colnames = [*VI_dict.keys()]
                    colnames_str = "\t".join(colnames)
                    writer.write(" " + "\t" + colnames_str + "\n")

                row = [str(vi) for vi in VI_dict.values()]
                row_str = "\t".join(row)
                writer.write(m1 + "\t" + row_str + "\n")


    def save_local_VI_scores(
                self, 
                local_VI_dict
        ):
        
        for m1, m1_vi_dict in local_VI_dict.items():
            m2_list, local_VI_scores = zip(*m1_vi_dict.items())
            m2_list = ["RefOGs"]
            vi_scores_list = []
            for m2, local_VI_scores in m1_vi_dict.items():
                if m1 != m2:
                    m2_list.append(m2)
                    refog_tuple, vi_scores = zip(*local_VI_scores.items())
                    vi_scores_list.append(vi_scores)
            
            vi_scores_arr = np.array(vi_scores_list).T 
            
            VI_filename = m1 + ".tsv"
            local_VI_filepath = self.vi_score_path / VI_filename
            local_VI_colname_str = "\t".join(m2_list)
            with open(local_VI_filepath, "w") as writer:
                writer.write(local_VI_colname_str + "\n")
                for refog_key, vi_score in zip(refog_tuple, vi_scores_arr):
                    row = [str(vi) for vi in vi_score]
                    row_str = "\t".join(row)
                    writer.write(refog_key + "\t" + row_str + "\n")


    def save_dist_metrics(
        self, 
        metric: str,
        score: str,
        metric_arr, 
        methods: List[str],
    ):
        filename = "_".join(score.lower().split()) + f"_{metric}.tsv"
        metric_filepath = self.dist_path / filename
        with open(metric_filepath, "w") as writer:
            colnames_str = "\t".join(methods)
            writer.write(" " + "\t" + colnames_str + "\n")
            for method, metric_row in zip(methods, metric_arr):
                row = [str(vi) for vi in metric_row]
                row_str = "\t".join(row)
                writer.write(method + "\t" + row_str + "\n")