import os
import sys
import re
import traceback
from typing import Any, Dict, List, Literal, Optional, Tuple

from inspect import getmembers, isfunction
from rich import print

CMD_MANAGER: Literal["preprocess", "benchmark"] = "preprocess"


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

        if args[0] in ["preprocess", "benchmark"]:
            print_help(args[0])
            sys.exit()
        else:
            print(
                "Missing/Incorrect task definition argument!"
                " Please choose from 'preprocess' and 'benchmark'."
            )
            sys.exit(1)

    elif len(args) == 2:
        if (args[0] in ["preprocess", "benchmark"]) \
            and (args[1] == "--help" or args[1] == "-h"):
            print_help(args[0])
            sys.exit()
        elif args[0] not in ["preprocess", "benchmark"]:
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
        "sonic": "sonicparanoid",
        "broccoli": "broccoli",
        "broc": "broccoli",
        "hier": "hieranoid",
        "hieranoid": "hieranoid",
        "omcl": "orthomcl",
        "orthomcl": "orthomcl",
    }
    return func_name


def flatten_list_of_list(alist: List[List[str]]):
    flat_list = []
    for row in alist:
        flat_list += row
    return flat_list
