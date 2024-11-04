import argparse
import os
import pathlib
import re

from ete3 import Tree


def hieranoid(input_file, output_file, protemes_dir: pathlib.Path):

    species_gene_dict = {}
    print(protemes_dir)
    for file in protemes_dir.iterdir():
        species = file.name.split(".", 1)[0]
        with open(file, "r") as reader:
            for line in reader:
                if ">" in line:
                    gene = line[1:].strip().split(".", 1)[-1]
                    species_gene_dict[gene] = species

    with open(output_file, "w") as writer:
        with open(input_file, "r") as reader:
            for line in reader:
                line = line.strip()
                og_name, og_tree = line.split(":", 1)
                t = Tree(og_tree, format=1)
                ogs = set()
                for node in t.traverse("postorder"):
                    if node.is_leaf():
                        # gene_name = node.name
                        # species_name = [species for species, genes in species_gene_dict.items() if gene_name in genes]
                        # ogs.add(node.name)
                        ogs.add(species_gene_dict[node.name] + "." + node.name)

                og = og_name + ": " + ", ".join(ogs)
                writer.write(og + "\n")


def orthomcl(input_file, output_file):
    with open(output_file, "w") as writer:
        with open(input_file, "r") as reader:
            single_gene_og = []
            for line in reader:
                line = line.strip()
                gname = line.split(": ", 1)[0]

                genes = [
                    item.strip().replace("|", ".").replace("adjusted_", "")
                    for item in line.split(": ", 1)[1].split(" ")
                    if "combined" not in item
                ]
                if len(genes) > 1:
                    og = gname + ": " + ", ".join(genes)
                    writer.write(og + "\n")
                else:
                    single_gene_og.extend(genes)
            if len(single_gene_og) > 1:
                single_gene_og_set = set(single_gene_og)
                og = "unassigned_genes" + ": " + ", ".join(single_gene_og_set)
            writer.write(og + "\n")


def orthofinder(input_file, output_file):
    old_version = True
    with open(output_file, "w") as writer:
        with open(input_file, "r") as reader:
            for i, line in enumerate(reader):
                line = line.strip()
                if i == 0:
                    if "Gene Tree Parent Clade" in line:
                        species_list = [
                            item.split(".", 1)[0] for item in line.split("\t", 3)[-1].split()
                        ]
                    else:
                        species_list = [
                            item.split(".", 1)[0] for item in line.split("\t", 1)[1].split()
                        ]
                        old_version = False
                else:
                    line = line.split("\t")
                    predog_key = line[0]
                    if old_version:
                        genes_list = line[3:]
                    else:
                        genes_list = line[1:]

                    genes = []
                    for i, gene in enumerate(genes_list):
                        if len(gene) == 0:
                            continue
                        if species_list[i] not in gene:
                            genes.append(
                                ", ".join(
                                    [
                                        species_list[i] + "." + g
                                        for g in re.split(" |,", gene)
                                        if len(g) != 0
                                    ]
                                )
                            )
                        else:
                            genes.append(gene)

                    og = predog_key + ": " + ", ".join(genes)
                    writer.write(og + "\n")


def sonicparanoid(input_file, output_file):
    with open(output_file, "w") as writer:
        with open(input_file, "r") as reader:
            for i, line in enumerate(reader):
                line = line.strip()
                if i == 0:
                    species_list = [item.split(".", 1)[0] for item in line.split("\t", 1)[1].split()]
                else:
                    predog_key, genes_str = line.split("\t", 1)
                    genes_list = genes_str.split("\t")
                    genes = []
                    for i, gene in enumerate(genes_list):
                        if gene == "*":
                            continue
                        if species_list[i] not in gene:
                            genes.append(
                                ", ".join(
                                    [
                                        species_list[i] + "." + g
                                        for g in re.split(" |,", gene)
                                        if len(g) != 0
                                    ]
                                )
                            )
                        else:
                            genes.append(gene)
                    og = predog_key + ": " + ", ".join(genes)
                    writer.write(og + "\n")


def broccoli(input_file, output_file):
    with open(output_file, "w") as writer:
        with open(input_file, "r") as reader:
            for i, line in enumerate(reader):
                line = line.strip()
                if i == 0:
                    species_list = [
                        item.split(".", 1)[0] for item in line.split("\t", 1)[1].split()
                    ]
                else:
                    predog_key, genes_str = line.split("\t", 1)
                    genes_list = genes_str.split("\t")
                    genes = []
                    for i, gene in enumerate(genes_list):
                        if gene == "":
                            continue
                        if species_list[i] not in gene:
                            genes.append(
                                ", ".join(
                                    [
                                        species_list[i] + "." + g
                                        for g in re.split(" |,", gene)
                                        if len(g) != 0
                                    ]
                                )
                            )
                        else:
                            genes.append(gene)
                    og = predog_key + ": " + ", ".join(genes)
                    writer.write(og + "\n")

