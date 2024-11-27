import pathlib
import os 
import re
import sys
from typing import Optional, List, Dict, Any, Set
import numpy as np

class OrthoBenchOGReader:

    def __init__(self, refog_path, uncerain_refog_path):
        self.refog_path = refog_path
        self.uncertain_refog_path = uncerain_refog_path

    def read_orthobench_refogs(self) -> Dict[str, Set[str]]:
        n_expected = 1945
        n_genes = 0
        refogs = {}
        refogs_ngenes_dict = {}
        refogs_nspecies_dict = {}
        try:
            with open(self.refog_path, "r") as reader:
                for line in reader:
                    if "Orthogroups" in line:
                        continue
                    line = line.strip().split("\t")
                    refog_key, nspecies, ngenes, genes = line[0], line[1], line[2], line[-1]
                    genes = [g.strip() for g in genes.split(", ") if g.strip() != ""]
                    refogs[refog_key] = set(genes)
                    refogs_ngenes_dict[refog_key] = int(ngenes)
                    refogs_nspecies_dict[refog_key] = int(nspecies)
                    n_genes += len(refogs[refog_key])

        except:
            print("ERROR: RefOG file not found in: %s" % self.refog_path)
            sys.exit(1)

        if "uncertain" in self.refog_path.name:
            assert (
                n_expected == n_genes
            ), "ERROR: There are genes missing from the RefOG files. Found %d, expected %d" % (
                n_genes,
                n_expected,
            )

        return refogs, refogs_nspecies_dict, refogs_ngenes_dict

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
            self.other_info_path = self.global_score_path 

        else:
            self.local_score_path = self.global_score_path / "local_scores"
            self.other_info_path = self.global_score_path / "other_info"
            self.vi_score_path = self.global_score_path / "local_vi_scores"
            self.dist_path = self.global_score_path / "local_score_distance"
             
            if not os.path.exists(self.local_score_path):
                os.makedirs(self.local_score_path, exist_ok=True)
            
            if not os.path.exists(self.other_info_path):
                os.makedirs(self.other_info_path, exist_ok=True)
            
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

class FileWriter:

    def __init__(self, 
                 global_score_path: pathlib.Path, 
                 local_score_path: pathlib.Path,
                 other_info_path: pathlib.Path,
                 vi_score_path: pathlib.Path,
                 dist_path: pathlib.Path,
                 ):
        self.global_score_path = global_score_path
        self.local_score_path = local_score_path
        self.other_info_path = other_info_path
        self.vi_score_path = vi_score_path
        self.dist_path = dist_path

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
    
    def save_local_scores(self, 
                          refog_path: str,
                          filename: str, 
                          scores_dict: Dict[str, Any],
                          precision: int,
                          ):
        refogs_dict = {}
        colnames = []
        with open(refog_path, "r") as reader:
            for line in reader:
                line = line.strip().split("\t")
                if "Orthogroups" in line:
                    colnames = line[:-1]
                    
                else:
                    refog_key = line[0]
                    refogs_dict[refog_key] = line[1:-1]

        local_score_names, local_scores = [*zip(*scores_dict.items())] 
        colnames.extend(local_score_names)           
        colnames_str = "\t".join(colnames)
        score_filename = self.local_score_path / filename 
     
        for local_score in local_scores:
            if local_score is None:
                continue
            for refog_key, scores in local_score.items():
                refogs_dict[refog_key].append(np.round(scores, precision))

        scores_lists = [
            [refog_key] + scores
            for refog_key, scores in refogs_dict.items()
        ]

        scores_lists = sorted(scores_lists, key=lambda x: x[3])
        with open(score_filename, "w") as writer:
            writer.write(colnames_str + "\n")
            for score_list in scores_lists:
                scores_str_list = [str(s) for s in score_list]
                scores_str = "\t".join(scores_str_list) 
                writer.write(scores_str + "\n")


    def save_missing_genes(self, 
                           filename: str, 
                           missing_species_dict: Dict[str, Set[str]],
                           missing_genes_dict: Dict[str, Set[str]]):
        
        colnames = ["RefOGs", "Missing Species", "Missing Genes"]
        colnames_str = "\t".join(colnames)
        other_info_filepath = self.other_info_path / filename
        with open(other_info_filepath, "w") as writer:
            writer.write(colnames_str + "\n")
            for refog_key, missing_genes in missing_genes_dict.items():
                missing_species = missing_species_dict[refog_key]
                missing_species_str = ", ".join(missing_species)
                missing_genes_str = ", ".join(missing_genes)
                writer.write("\t".join((refog_key, missing_species_str, missing_genes_str)) + "\n")


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