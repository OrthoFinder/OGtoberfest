import pathlib
from collections import Counter
from typing import Optional, List, Set, Tuple, Dict 
from ogtoberfest import utils



def grouped_orthologues(orthologues_dir: pathlib.Path):

    for file in orthologues_dir.iterdir():
        with open(file) as reader:
            pass