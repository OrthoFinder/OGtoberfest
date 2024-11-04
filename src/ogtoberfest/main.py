import sys
from typing import Dict, List, Optional
import re

from ogtoberfest import  preprocess, process_args, utils, benchmark

def main(args: Optional[List[str]] = None):

    if not args:
        args = sys.argv[1:]

    args = utils.check_cmd_args(args)
    # filehandler = files.FileHandler()
    task = args.pop(0)

    manager = process_args.create_options(args, task=task)
    func_name_maps = utils.get_func_name_map()
    funcs_dict = utils.get_func_name(preprocess)
    for file in manager.options.input_path.iterdir():
        method = re.split("_|\.", file.name)[0]
        method = method.lower()
        method_func_name = func_name_maps.get(method)

        if method_func_name is not None:
            method_func = funcs_dict[method_func_name]
            if method_func_name == "hieranoid":
                if manager.options.database_path is None:
                    print("Hieranoid needs to provide a database!")
                    continue
                else:
                    method_func(
                        file,
                        manager.options.output_path / file.name,
                        manager.options.database_path,
                    )
            else:
                method_func(
                    file,
                    manager.options.output_path / file.name,
                )

    print()


if __name__ == "__main__":
    main()
