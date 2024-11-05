import json
import os
import pathlib
from typing import Any, Dict, List, Literal, Optional, Tuple

from ogtoberfest import utils

CMD_MANAGER: Literal["preprocess", "benchmark"] = "preprocess"

THIS_DIR = pathlib.Path(__file__).parent
ARGS_JSON_PATH = THIS_DIR / "./input_args.json"

CWD = pathlib.Path.cwd()

class Manager:
    def __init__(self, task: CMD_MANAGER):
        self.task = task
        self.options = Options()

    def __repr__(self):
        attributes_dict = {"task": self.task}
        return f"Manager({attributes_dict})"

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


def create_options(args: List[str], task=CMD_MANAGER) -> Manager:
    manager = Manager(task)
    output_path_name = ""
    show_args = False
    while args:
        arg = args.pop(0)

        if arg == "--show-args":
            show_args = True

        arg_dict = read_args_from_json(task, arg)  # Presumed function to fetch argument info from JSON

        if arg_dict is not None:
            if not args:
                print(f"Missing option for argument {arg}")
                utils.fail()

            attr_name = arg_dict["name"]
            arg_value = args.pop(0)

            if "path" in attr_name:
                if os.path.isfile(arg_value):
                    arg_value = pathlib.Path(utils.get_file_arg(arg_value)).resolve()
                elif os.path.isdir(arg_value):
                    arg_value = pathlib.Path(utils.get_dir_arg(arg_value)).resolve()
                else:
                    if attr_name == "input_path":
                        print(f"Error: '{arg_value}' is not a valid path or it does not exist!")
                        utils.fail()
                    elif attr_name == "output_path":
                        output_path_name = arg_value
                        continue
            if attr_name in ["outgroups", "additional_species", "input_species"]:
                arg_value = [
                    utils.curate_labels(item.strip()) 
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
                if arg_value not in arg_dict["valid_options"]:
                    print(arg_dict["error_message"])
                    print(arg_dict["valid_options"])
                    utils.fail()
            setattr(manager.options, attr_name, arg_value)

    if task == "preprocess":
        prefix = "preprocessed"
    else:
        prefix = "scores"

    handle_missing_path_args(manager, output_path_name, prefix, task)

    if manager.options.input_path.parent.name == "OrthoBench" and task == "benchmark":
        if not hasattr(manager.options, "outgroups"):
            arg_dict = read_args_from_json(task, "--outgroups")
            manager.options.outgroups = arg_dict["default"]
        if not hasattr(manager.options, "additional_species"):
            arg_dict = read_args_from_json(task, "--additional-species")
            manager.options.additional_species = arg_dict["default"]

    if show_args:
        print("Here is your input arguments:")
        manager.options.show_options()
        print()

    return manager


def read_args_from_json(task=CMD_MANAGER, flag: str = "") -> Optional[Dict[str, str]]:
    args_contents = ARGS_JSON_PATH.read_text()
    args_list = json.loads(args_contents)[task]
    for arg_dict in args_list:
        if flag in arg_dict["flag"]:
            return arg_dict
    return


def handle_missing_path_args(manager: Manager, output_path_name: str, prefix: str, task=CMD_MANAGER,):

    input_path_exist = True
    if not hasattr(manager.options, "input_path"):
        input_path_exist = False
        arg_dict = read_args_from_json(task, "--input")
        input_path = CWD / arg_dict["default"]
        setattr(manager.options, "input_path", input_path.resolve())

    output_path_isdir = False
    output_path: Optional[pathlib.Path] = None
    if not hasattr(manager.options, "output_path"):
        if not input_path_exist:
            arg_dict = read_args_from_json(task, "--output")
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

    if not hasattr(manager.options, "database_path"):
        if not input_path_exist \
            or manager.options.input_path.parent.name == "OrthoBench":
            arg_dict = read_args_from_json(task, "--database")
            database_path = CWD / arg_dict["default"]
            setattr(manager.options, "database_path", database_path.resolve())

    if task == "benchmark":
        if not hasattr(manager.options, "refog_path"):
            if manager.options.input_path.parent.name == "OrthoBench":
                arg_dict = read_args_from_json(task, "--refog")
                refog_path = CWD / arg_dict["default"]
                setattr(manager.options, "refog_path", refog_path.resolve())
