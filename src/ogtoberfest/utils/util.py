import os
import sys
import re
import traceback
from typing import Any, Dict, List, Literal, Optional, Tuple

from inspect import getmembers, isfunction
from rich import print

CMD_MANAGER: Literal[
    "orthogroups_preprocess", 
    "orthologues_preprocess", 
    "orthogroups_benchmark",
    "orthologues_benchmark",
] = "orthogroups_preprocess"



class SpeciesInfo(object):
    def __init__(self):
        self.speciesToUse = []  #       seqsInfo.iSpeciesToUse   - which to include for this analysis
        self.nSpAll = None  #       seqsInfo.nSpAll => 0, 1, ..., nSpAll - 1 are valid species indices
        self.iFirstNewSpecies = None  #       iFirstNew   => (0, 1, ..., iFirstNew-1) are from previous and (iFirstNew, iFirstNew+1, ..., nSpecies-1) are the new species indices
        self.speciesToUseName = []
        self.sequence2id_dict = {}
        self.id2sequence_dict = {}
        self.id2species_dict = {}
        self.species2id_dict = {}
        
    def __str__(self):
        return str((self.speciesToUse, self.nSpAll, self.iFirstNewSpecies))

    def get_original_species(self):
        if self.iFirstNewSpecies is None:
            return self.speciesToUse
        else:
            return [iSp for iSp in self.speciesToUse if iSp < self.iFirstNewSpecies]

    def get_species2ids(self):

        if len(self.speciesToUseName) != 0 and len(self.speciesToUse) != 0:
            self.species2id_dict = {
                name: species_id
                for name, species_id in zip(self.speciesToUseName, self.speciesToUse)
            }
        else:
            self.species2id_dict = {}

        for name, species_id in zip(self.speciesToUseName, self.speciesToUse):
            print(name, species_id)
        return self.species2id_dict
    
    def get_ids2species(self):
        if len(self.speciesToUseName) != 0 and len(self.speciesToUse) != 0:
            self.id2species_dict = {
                species_id: name
                for species_id, name in zip(self.speciesToUse, self.speciesToUseName)
            }
        else:
            self.id2species_dict = {}
        return self.id2species_dict


    def get_species_to_use(self, species_ids_fn):
        """Returns species indices (int) to use and total number of species available"""
        species_to_use = []
        species_to_use_names = []
        nSkipped = 0
        with open(species_ids_fn, "r") as speciesF:
            for line in speciesF:
                line = line.rstrip()
                if not line:
                    continue
                if line.startswith("#"):
                    nSkipped += 1
                else:
                    iSp, spName = line.split(": ")
                    species_to_use.append(int(iSp))
                    species_to_use_names.append(spName)
        return species_to_use, len(species_to_use) + nSkipped, species_to_use_names


def curate_labels(label: str):
    COMPILE = re.compile(r"[_\s]+")
    label = COMPILE.sub(" ", label).split()
    label = "_".join(label)

    return label


def check_cmd_args(args: List[str]) -> List[str]:

    if (not args or len(args) == 0) \
        or (args[0] == "--help" or args[0] == "help" or args[0] == "-h"):
        print_help()
        sys.exit()

    if len(args) == 1:
        if args[0] == "-v" or args[0] == "--version":
            print("0.0.0")
            sys.exit()

        if args[0] in [
            "orthogroups_preprocess", 
            "orthologues_preprocess", 
            "orthogroups_benchmark",
            "orthologues_benchmark"
            ]:
            print_help(args[0])
            sys.exit()
        else:
            print(
                "Missing/Incorrect task definition argument!"
                " Please choose from 'orthogroups_preprocess', 'orthologues_preprocess' and 'benchmark'."
            )
            sys.exit(1)

    elif len(args) == 2:
        if (args[0] in [            
            "orthogroups_preprocess", 
            "orthologues_preprocess", 
            "orthogroups_benchmark",
            "orthologues_benchmark"
            ]) \
            and (args[1] == "--help" or args[1] == "-h"):
            print_help(args[0])
            sys.exit()
        elif args[0] not in [
            "orthogroups_preprocess", 
            "orthologues_preprocess", 
            "orthogroups_benchmark",
            "orthologues_benchmark"
            ]:
            print(
                "Unrecognised input arguments!"
            )
            print_help()
            sys.exit(1)

    return args


def print_help(task=Optional[CMD_MANAGER]):

    print("Thank you for using OGtoberfest!")
    pass


def fail():
    sys.stderr.flush()
    print(traceback.format_exc())
    print_help()
    sys.exit(1)


def get_dir_arg(arg: str) -> str:
    directory = os.path.abspath(arg)
    if not os.path.exists(directory):
        print("Specified directory doesn't exist: %s" % directory)
        fail()
    if not os.path.isfile(directory) and directory[-1] != os.sep:
        directory += os.sep
    return directory


def get_file_arg(arg: str) -> str:
    file_path = os.path.abspath(arg)
    directory = os.path.dirname(file_path)
    if not os.path.exists(directory):
        print("Directory points to the file doesn't exist: %s" % directory)
        fail()
    if not os.path.isfile(file_path):
        print("Specified file doesn't exist: %s" % file_path)
        fail()
    return file_path

def get_func_name(object_name):
    func_dict = dict(getmembers(object_name, isfunction))
    return func_dict

def get_func_name_map():
    func_name = {
        "orthofinder": "orthofinder",
        "orthofinder2": "orthofinder",
        "orthofinder3": "orthofinder",
        "of2": "orthofinder",
        "of3": "orthofinder",
        "of": "orthofinder",
        "sonicparanoid": "sonicparanoid",
        "sonicparanoid2": "sonicparanoid",
        "sp2": "sonicparanoid",
        "sp": "sonicparanoid",
        "sonic": "sonicparanoid",
        "broccoli": "broccoli",
        "broc": "broccoli",
        "hier": "hieranoid",
        "hieranoid": "hieranoid",
        "omcl": "orthomcl",
        "orthomcl": "orthomcl",
        "po": "proteinortho",
        "prot": "proteinortho",
        "proteinortho": "proteinortho",
        "fastoma": "fastoma",
        "fast": "fastoma",
        "fo": "fastoma",
        "so": "swiftortho",
        "swiftortho": "swiftortho",
        "orthohmm": "orthohmm",
    }
    return func_name


def flatten_list_of_list(alist: List[List[str]]):
    flat_list = []
    for row in alist:
        flat_list += row
    return flat_list
