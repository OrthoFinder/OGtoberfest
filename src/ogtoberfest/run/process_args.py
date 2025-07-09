from __future__ import annotations
import json
import os
import pathlib
from typing import Any, Dict, List, Literal, Optional, Tuple, TypedDict, TypeVar, Union
from dataclasses import dataclass
from ogtoberfest.utils import util

CMD_MANAGER: Literal[
    "orthogroups_preprocess", 
    "orthologues_preprocess", 
    "orthogroups_benchmark",
    "orthologues_benchmark",
    ] = "orthogroups_preprocess"
ArgValue = Any

THIS_DIR = pathlib.Path(__file__).parent
ARGS_JSON_PATH = THIS_DIR / "./input_args.json"

CWD = pathlib.Path.cwd()

# class Manager:
#     def __init__(self, task: CMD_MANAGER):
#         self.task = task
#         self.options = Options()

#     def __repr__(self):
#         attributes_dict = {"task": self.task}
#         return f"Manager({attributes_dict})"

@dataclass
class Manager:
    task: str
    options: Options

class Options:
    def __init__(self):
        self._dynamic_attributes = {}

    def __getattr__(self, name: str) -> Any:
        if name in self._dynamic_attributes:
            return self._dynamic_attributes[name]

        raise AttributeError(
            f"Error: '{name}' is not a valid attribute of '{type(self).__name__}'. "
            f"Available attributes are: {list(self._dynamic_attributes.keys())}"
        )

    def __setattr__(self, name: str, value: Any) -> None:
        if name.startswith("_"):
            super().__setattr__(name, value)
        else:
            self._dynamic_attributes[name] = value

    def __delattr__(self, name: str) -> None:
        if name in self._dynamic_attributes:
            del self._dynamic_attributes[name]
        else:
            super().__delattr__(name)

    def __hasattr__(self, name):
        return name in self._dynamic_attributes or hasattr(super(), name)

    def show_options(self) -> None:
        for key, value in self._dynamic_attributes.items():
            print(f"{key}: {value}")

    def __repr__(self):
        attributes_dict = {k: v for k, v in self.__dict__.items() if k != "_dynamic_attributes"}
        attributes_dict.update(self._dynamic_attributes)

        return f"Options({attributes_dict})"


def create_options(args: List[str], task=Optional[CMD_MANAGER]):
    # manager = Manager(task)
    manager = Manager (
        **{
        "task": task,
        "options": Options()
        }
    )
    # manager = Manager(task, Options())
    args_list = read_args_from_json(task)
    output_path_name = ""
    show_args = False
    usr_input_args = []
    while args:
        arg = args.pop(0)
 
        if arg == "--show-args":
            show_args = True

        arg_dict = read_args_from_arg_list(args_list, arg)
        if arg_dict is not None:
            usr_input_args.append(arg)
            if not args:
                print(f"Missing option for argument {arg}")
                util.fail()

            attr_name = arg_dict["name"]
            arg_value: ArgValue = args.pop(0)

            if "path" in attr_name:
                if os.path.isfile(arg_value):
                    arg_value = pathlib.Path(util.get_file_arg(arg_value)).resolve()
                elif os.path.isdir(arg_value):
                    arg_value = pathlib.Path(util.get_dir_arg(arg_value)).resolve()
                else:
                    if attr_name == "input_path":
                        print(f"Error: '{arg_value}' is not a valid path or it does not exist!")
                        util.fail()
                    elif attr_name == "output_path":
                        output_path_name = arg_value
                        continue
            if attr_name in ["outgroups", "additional_species", "input_species"]:
                arg_value = [
                    util.curate_labels(item.strip())
                    for item in arg_value.strip().split(",")
                    if len(item) != 0
                ]

            if "," in str(arg_value):
                arg_value = [
                    item.strip() for item in arg_value.strip().split(",") if len(item) != 0
                ]

            # need further action
            if len(arg_dict["limit"]) != 0:
                limit = arg_dict["limit"]
                print(limit)

            if len(arg_dict["valid_options"]) != 0:
                if isinstance(arg_value, list):
                    for val in arg_value:
                        if val not in arg_dict["valid_options"]:
                            print(arg_dict["error_message"])
                            print(f"\nSupported arg for {arg}:")
                            print(arg_dict["valid_options"])
                            util.fail()
                else:
                    if arg_value not in arg_dict["valid_options"]:
                        print(arg_dict["error_message"])
                        print(f"\nSupported option for flag {arg}:")
                        print(arg_dict["valid_options"])
                        util.fail()
            setattr(manager.options, attr_name, arg_value)

    if task in ["orthogroups_preprocess", "orthologues_preprocess"]:
        prefix = "preprocessed"
    else:
        prefix = "scores"

    handle_missing_path_args(manager, output_path_name, prefix, args_list, task)
    
    # if manager.options.input_path.parent.name == "OrthoBench":# and task == "benchmark":
    read_default_args(manager, args_list, usr_input_args)
        # if not hasattr(manager.options, "outgroups"):
        #     arg_dict = read_args_from_arg_list(args_list, "--outgroups")
        #     manager.options.outgroups = arg_dict["default"]
        # if not hasattr(manager.options, "additional_species"):
        #     arg_dict = read_args_from_arg_list(args_list, "--additional-species")
        #     manager.options.additional_species = arg_dict["default"]

        # if not hasattr(manager.options, "rank_method"):
        #     arg_dict = read_args_from_arg_list(args_list, "--rank-method")
        #     manager.options.rank_method = arg_dict["default"]

    if show_args:
        print("Here is your input arguments:")
        manager.options.show_options()
        print()

    return manager


def read_args_from_json(task=CMD_MANAGER):
    args_contents = ARGS_JSON_PATH.read_text()
    args_list = json.loads(args_contents)[task]
    return args_list


def read_args_from_arg_list(args_list, flag: str) -> Any:
    for arg_dict in args_list:
        if flag in arg_dict["flag"]:
            return arg_dict


def read_default_args(manager, args_list, usr_input_args: List[str]):
    for arg_dict in args_list:
        flags = arg_dict["flag"]
        exist_flag = [flag for flag in flags if flag in usr_input_args]
        if len(exist_flag):
            continue
        arg_name = arg_dict["name"]
        arg_value = None if len(arg_dict["default"]) == 0  else arg_dict["default"]
        if arg_name == "precision" and arg_value is not None:
            arg_value = int(arg_value)

        if "input_path" == arg_name or "output_path" == arg_name:
            continue
        elif "path" in arg_name:
            arg_value = CWD / arg_value
            arg_value = arg_value.resolve()

        setattr(manager.options, arg_name, arg_value)


def handle_missing_path_args(
        manager: Manager,
        output_path_name: str,
        prefix: str,
        args_list,
        task,
    ):

    input_path_exist = True
    if not hasattr(manager.options, "input_path"):
        input_path_exist = False
        arg_dict = read_args_from_arg_list(args_list, "--input")
        input_path = CWD / arg_dict["default"]
        setattr(manager.options, "input_path", input_path.resolve())

    output_path_isdir = False
    output_path: Optional[pathlib.Path] = None
    if not hasattr(manager.options, "output_path"):
        if not input_path_exist:
            arg_dict = read_args_from_arg_list(args_list, "--output")
            default_dir = CWD / arg_dict["default"]
            output_path = default_dir.resolve() / manager.options.input_path.name
            output_path_isdir = True

        if manager.options.input_path.is_file():
            if len(output_path_name) == 0:
                output_path_filename = \
                    "_".join((prefix, manager.options.input_path.name.split(".")[0]))
            else:
                output_path_filename = \
                    "_".join((prefix, output_path_name + ".tsv"))
            output_path = manager.options.input_path.parent / output_path_filename

        elif manager.options.input_path.is_dir():
            if len(output_path_name) == 0:
                output_path_dirname = \
                    "_".join((prefix, manager.options.input_path.name))
            else:
                output_path_filename = "_".join((prefix, output_path_name))
            output_path = manager.options.input_path.parent / output_path_dirname
            output_path_isdir = True

        setattr(manager.options, "output_path", output_path)

    if output_path_isdir and not os.path.exists(manager.options.output_path):
        os.makedirs(manager.options.output_path, exist_ok=True)
    
    setattr(manager.options, "wd_base", manager.options.input_path.parent / "WorkingDirectory")
    os.makedirs(manager.options.wd_base, exist_ok=True)

    # if not hasattr(manager.options, "database_path"):
    #     if not input_path_exist \
    #         or manager.options.input_path.parent.name == "OrthoBench":
    #         arg_dict = read_args_from_arg_list(args_list, "--database")
    #         database_path = CWD / arg_dict["default"]
    #         setattr(manager.options, "database_path", database_path.resolve())

    # if task == "benchmark":
    #     if not hasattr(manager.options, "refog_path"):
    #         if manager.options.input_path.parent.name == "OrthoBench":
    #             arg_dict = read_args_from_arg_list(args_list, "--refog")
    #             refog_path = CWD / arg_dict["default"]
    #             setattr(manager.options, "refog_path", refog_path.resolve())
