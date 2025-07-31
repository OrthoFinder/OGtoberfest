import argparse
import os
import itertools
import pathlib
import re
from typing import Optional, Dict


def orthofinder(
        input_file, 
        output_file, 
        sequence2id_dict: Optional[Dict[str, str]] = None,
    ):

    with open(output_file, "a") as writer:
        with open(input_file, "r") as reader:
            for n, line in enumerate(reader):
                line = line.strip().split("\t")
                if n == 0:
                    base_species = line[2]
                else:
                    other_species = line[1]
                    other_genes = [
                        sequence2id_dict.get(other_species + "." + gene)
                        if other_species.lower() not in gene.lower()
                        else sequence2id_dict.get(gene)
                        for gene in line[-1].split(", ")
                    ]
                    base_genes = [
                        sequence2id_dict.get(base_species + "." + gene)
                        if base_species.lower() not in gene.lower()
                        else sequence2id_dict.get(gene)
                        for gene in line[2].split(", ")
                    ]

                    for gene1, gene2 in itertools.product(base_genes, other_genes):
                        writer.write(gene1 + ", " + gene2 + "\n")


def broccoli(
        input_file, 
        output_file, 
        sequence2id_dict: Optional[Dict[str, str]] = None,
    ):

    with open(output_file, "w") as writer:
        with open(input_file, "r") as reader:
            for n, line in enumerate(reader):
                gene1, gene2 = line.strip().split("\t")
                gene1_id = sequence2id_dict.get(gene1)
                gene2_id = sequence2id_dict.get(gene2)
                writer.write(gene1_id + ", " + gene2_id + "\n")


def sonicparanoid(
        input_file, 
        output_file, 
        sequence2id_dict: Optional[Dict[str, str]] = None,
    ):
    with open(output_file, "a") as writer:
        with open(input_file, "r") as reader:
                
            for n, line in enumerate(reader):
                if n == 0:
                    continue
                else:
                    line = re.split("\t", line.strip())
                    base_genes = [
                        sequence2id_dict.get(item)
                        for i, item in enumerate(line[1].split())
                        if i % 2 == 0
                    ]

                    other_genes = [
                        sequence2id_dict.get(item)
                        for i, item in enumerate(line[-1].split())
                        if i % 2 == 0
                    ]

                    for gene1, gene2 in itertools.product(base_genes, other_genes):
                        writer.write(gene1 + ", " + gene2 + "\n")

