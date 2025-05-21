import pathlib
import os 
import re
import sys
from typing import Optional, List, Dict, Any, Set
import numpy as np
import ete3
from ogtoberfest.utils import util


class FileHandler:
    def __init__(
            self, 
            input_path: pathlib.Path,
            working_dir: pathlib.Path,
            output_path: pathlib.Path, 
            input_path_isfile: bool = False,
            method: str = "OrthoFinder3",
            msa_method: str = "mafft",
            tree_method: str = "fasttree",
    ):
        
        self.wd_current = working_dir / method
        if not os.path.exists(self.wd_current):
            os.makedirs(self.wd_current, exist_ok=True)

        self.msa_dir = self.wd_current / msa_method 
        if not os.path.exists(self.msa_dir):
            os.makedirs(self.msa_dir, exist_ok=True)

        self.tree_dir = self.wd_current / tree_method
        if not os.path.exists(self.tree_dir):
            os.makedirs(self.tree_dir, exist_ok=True)  

        self.fasta_dir = self.wd_current / "orthologues_sequences"
        if not os.path.exists(self.fasta_dir):
            os.makedirs(self.fasta_dir, exist_ok=True)  


    def get_species_ids_fn(self): 
        return self.wd_current + "SpeciesIDs.txt"
        
    def get_sequence_ids_fn(self):
        # It is always in the first of the 'extension' directories (as this is the relevant one)
        return self.wd_current + "SequenceIDs.txt"

    def get_msa_dir(self):
        return self.msa_dir 
    
    def get_tree_dir(self):
        return self.tree_dir

    def get_fasta_dir(self):
        return self.fasta_dir

    def get_species_fasta_fn(self, iSpecies, qForCreation=False):
        """
        qForCreation: A path is required at which the file should be created (don't search for it)
        """
        if len(self.wd_base) == 0: raise Exception("No wd1")
        fn = "%sSpecies%d.fa" % (self.wd_current, iSpecies)
        if qForCreation:
            return fn
        if os.path.exists(fn): 
            return fn
        raise Exception(fn + " not found")



class FileReader:

    def read_trees(self, tree_dir):
        pass 



    def ProcessesNewFasta(fastaDir, q_dna, speciesInfoObj_prev = None, speciesToUse_prev_names=[]):
        """
        Process fasta files and return a Directory object with all paths completed.
        """
        fastaExtensions = {"fa", "faa", "fasta", "fas", "pep", "fna"}
        # Check files present
        qOk = True
        if not os.path.exists(fastaDir):
            print("\nDirectory does not exist: %s" % fastaDir)
            util.Fail()
        files_in_directory = sorted([f for f in os.listdir(fastaDir) if os.path.isfile(os.path.join(fastaDir,f))])
        originalFastaFilenames = []
        excludedFiles = []

        for f in files_in_directory:
            if len(f.rsplit(".", 1)) == 2 and f.rsplit(".", 1)[1].lower() in fastaExtensions and not f.startswith("._"):
                originalFastaFilenames.append(f)
            else:
                excludedFiles.append(f)

        if len(excludedFiles) != 0:
            print("\nWARNING: Files have been ignored as they don't appear to be FASTA files:")
            for f in excludedFiles:
                print(f)
            print("OrthoFinder expects FASTA files to have one of the following extensions: %s" % (", ".join(fastaExtensions)))
        
        speciesToUse_prev_names = set(speciesToUse_prev_names)
        if len(originalFastaFilenames) + len(speciesToUse_prev_names) < 2:
            print("ERROR: At least two species are required")
            util.Fail()

        if any([fn in speciesToUse_prev_names for fn in originalFastaFilenames]):
            print("ERROR: Attempted to add a second copy of a previously included species:")
            for fn in originalFastaFilenames:
                if fn in speciesToUse_prev_names: print(fn)
            print("")
            util.Fail()

        if len(originalFastaFilenames) == 0:
            print("\nNo fasta files found in supplied directory: %s" % fastaDir)
            util.Fail()

        if speciesInfoObj_prev == None:
            # Then this is a new, clean analysis 
            speciesInfoObj = util.SpeciesInfo()
        else:
            speciesInfoObj = speciesInfoObj_prev

        iSeq = 0
        iSpecies = 0
        # If it's a previous analysis:
        if len(speciesToUse_prev_names) != 0:
            with open(files.FileHandler.GetSpeciesIDsFN(), 'r') as infile:
                for line in infile: pass
            if line.startswith("#"): line = line[1:]
            iSpecies = int(line.split(":")[0]) + 1
        speciesInfoObj.iFirstNewSpecies = iSpecies

        newSpeciesIDs = []

        with open(files.FileHandler.GetSequenceIDsFN(), 'a') as idsFile, open(files.FileHandler.GetSpeciesIDsFN(), 'a') as speciesFile:
            for fastaFilename in originalFastaFilenames:
                newSpeciesIDs.append(iSpecies)
                outputFasta = open(files.FileHandler.GetSpeciesFastaFN(iSpecies, qForCreation=True), 'w')
                fastaFilename = fastaFilename.rstrip()
                speciesFile.write("%d: %s\n" % (iSpecies, fastaFilename))
                baseFilename, extension = os.path.splitext(fastaFilename)
                mLinesToCheck = 100
                qHasAA = False
                with open(fastaDir + os.sep + fastaFilename, 'r') as fastaFile:
                    for iLine, line in enumerate(fastaFile):
                        if line.isspace(): continue
                        if len(line) > 0 and line[0] == ">":
                            newID = "%d_%d" % (iSpecies, iSeq)
                            acc = line[1:].rstrip()
                            if len(acc) == 0:
                                print("ERROR: %s contains a blank accession line on line %d" % (fastaDir + os.sep + fastaFilename, iLine+1))
                                util.Fail()
                            idsFile.write("%s: %s\n" % (newID, acc))
                            outputFasta.write(">%s\n" % newID)    
                            iSeq += 1
                        else:
                            line = line.upper()    # allow lowercase letters in sequences
                            if not qHasAA and (iLine < mLinesToCheck):
    #                            qHasAA = qHasAA or any([c in line for c in ['D','E','F','H','I','K','L','M','N','P','Q','R','S','V','W','Y']])
                                qHasAA = qHasAA or any([c in line for c in ['E','F','I','L','P','Q']]) # AAs minus nucleotide ambiguity codes
                            outputFasta.write(line)
                    outputFasta.write("\n")
                if (not qHasAA) and (not q_dna):
                    qOk = False
                    print("ERROR: %s appears to contain nucleotide sequences instead of amino acid sequences. Use '-d' option" % fastaFilename)
                iSpecies += 1
                iSeq = 0
                outputFasta.close()
            if not qOk:
                util.Fail()

        if len(originalFastaFilenames) > 0: outputFasta.close()
        speciesInfoObj.speciesToUse = speciesInfoObj.speciesToUse + newSpeciesIDs
        speciesInfoObj.nSpAll = max(speciesInfoObj.speciesToUse) + 1      # will be one of the new species
        
        return speciesInfoObj


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

    def __init__(
        self, 
        input_path: pathlib.Path,
        working_dir: pathlib.Path,
        output_path: pathlib.Path, 
        input_path_isfile: bool = False,
        method: str = "OrthoFinder3",
        msa_method: str = "mafft",
        tree_method: str = "fasttree",
    ):
        super().__init__(
            input_path, 
            working_dir, 
            output_path, 
            input_path_isfile, 
            method, 
            msa_method, 
            tree_method
        )
        self.base_og_format = "OG%07d"
        self.orthologues_dir = self.input_path / self.method
        
        
    def combined_gene_pairs(self):
        combined_al_gene_pairs_fn = self.orthologues_dir / "all_gene_pairs_combined.txt"
        with open(combined_al_gene_pairs_fn, "w") as writer:
            for orthologues_file in self.orthologues_dir.iterdir():
                if orthologues_file == combined_al_gene_pairs_fn:
                    continue 
                with open(orthologues_file) as reader:
                    for line in reader:
                        writer.write(line)


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

    def save_local_info(
            self, 
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
                else:
                    writer.write(refog_key + "\n")

    def save_fission_genes(
            self, 
            filename: str, 
            fission_genes_dict: Dict[str, Set[str]]
        ):
        
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
                else:
                    writer.write(refog_key + "\n")

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


    def save_overlap_predogs(
            self, 
            filename: str, 
            predogs_dict: Dict[str, Set[str]]
        ):

        # colnames = ["PredOGs", "Genes"]
        # colnames_str = "\t".join(colnames)
        overlap_predog_filepath = self.overlap_predog_path / filename
        with open(overlap_predog_filepath, "w") as writer:
            # writer.write(colnames_str + "\n")
            for predog_key, predogs in predogs_dict.items():
                line = predog_key + ": " + ", ".join(predogs)
                writer.write(line + "\n")
