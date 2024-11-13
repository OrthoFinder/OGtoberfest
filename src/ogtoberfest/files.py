import pathlib
import os 
import re
import sys
from typing import Optional, List, Dict, Any, Set


class OrthoBenchOGReader:

    def __init__(self, refog_path, uncerain_refog_path):
        self.refog_path = refog_path
        self.uncertain_refog_path = uncerain_refog_path

    def read_orthobench_refogs(self) -> Dict[str, Set[str]]:
        n_expected = 1945
        n_genes = 0
        refogs = {}
        try:
            with open(self.refog_path, "r") as reader:
                for line in reader:
                    if "Orthogroups" in line:
                        continue
                    line = line.strip().split("\t")
                    refog_key, genes = line[0], line[-1]
                    genes = [g.strip() for g in genes.split(", ") if g.strip() != ""]
                    refogs[refog_key] = set(genes)
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

        return refogs

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


class FileWriter:

    def __init__(self, output_path):
        self.output_path = output_path

    def save_global_scores(self, 
                           colnames: List[str], 
                           filename: str, 
                           scores_dict: Dict[str, Any]):

        colnames_str = "\t".join(colnames)
        score_filename = self.output_path / filename 
        with open(score_filename, "w") as writer:
            writer.write(colnames_str + "\n")
            for method, scores in scores_dict.items():
                scores_str_list = [str(s) for s in scores]
                scores_str = "\t".join(scores_str_list)
                line = "\t".join((method, scores_str))
                writer.write(line + "\n")
