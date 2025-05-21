import pathlib
import os 
import re
import sys
from typing import Optional, List, Dict, Any, Set
import numpy as np


class QFO:

    def __init__(self, msa_cmd, tree_cmd):
        self.msa_cmd = msa_cmd 
        self.tree_cmd = tree_cmd


        