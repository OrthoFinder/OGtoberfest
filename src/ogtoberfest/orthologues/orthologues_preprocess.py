import argparse
import os
import itertools
import pathlib
import re
from typing import Optional
from ete3 import Tree


def hieranoid(input_file, output_file, protemes_dir: pathlib.Path):

    species_gene_dict = {}
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
                        species = species_gene_dict.get(node.name)
                        if species is not None:
                            ogs.add(species + "." + node.name)
                        else:
                            ogs.add(node.name)

                og = og_name + ": " + ", ".join(ogs)
                writer.write(og + "\n")


def orthomcl(input_file, output_file):
    with open(output_file, "w") as writer:
        with open(input_file, "r") as reader:
            single_gene_og = []
            for line in reader:
                line = line.strip()
                gname = line.split(": ", 1)[0]

                # genes = [
                #     item.strip().replace("|", ".").replace("adjusted_", "")
                #     for item in line.split(": ", 1)[1].split(" ")
                #     if "combined" not in item
                # ]

                genes = []
                for item in line.split(": ", 1)[1].split(" "):
                    if "combined" not in item:
                        gene = item.strip().replace("|", ".").replace("adjusted_", "")
                        genes.append(gene)

                if len(genes) > 1:
                    og = gname + ": " + ", ".join(genes)
                    writer.write(og + "\n")
            #     else:
            #         single_gene_og.extend(genes)
            # if len(single_gene_og) > 1:
            #     single_gene_og_set = set(single_gene_og)
            #     og = "unassigned_genes" + ": " + ", ".join(single_gene_og_set)
            # writer.write(og + "\n")



# def orthomcl(input_file, output_file):
#     with open(output_file, "w") as writer:
#         with open(input_file, "r") as reader:
             
#             predog_base_name = "my_prefix%03d"
#             n = 1000
#             for line in reader:
#                 line = line.strip()
#                 gname = line.split(": ", 1)[0]
#                 genes = []
#                 for item in line.split(": ", 1)[1].split(" "):
#                     if "combined" not in item:
#                         gene = item.strip().replace("|", ".").replace("adjusted_", "")
#                         if "Branchiostoma_lanceolatum" in gene or "Schistosoma_mansoni" in gene:
#                             continue 
#                         genes.append(item)
#                 if len(genes) > 0:
#                     predog_key = predog_base_name % n
#                     og = predog_key + ": " + " ".join(genes)
#                     writer.write(og + "\n")
#                     n += 1

def orthofinder(input_file, output_file):

    with open(output_file, "w+") as writer:
        with open(input_file, "r") as reader:
            for n, line in enumerate(reader):
                line = line.strip().split("\t")
                if n == 0:
                    base_species = line[2]
                else:
                    other_species = line[1]
                    other_genes = [
                        other_species + "." + gene
                        for gene in line[-1].split(", ")
                    ]
                    base_genes = [
                        base_species + "." + gene
                        for gene in line[2].split(", ")
                    ]

                    for gene1, gene2 in itertools.product(base_genes, other_genes):
                        writer.write(gene1 + ", " + gene2 + "\n")
                


def sonicparanoid(input_file, output_file):
    with open(output_file, "w") as writer:
        with open(input_file, "r") as reader:
            for n, line in enumerate(reader):
                line = line.strip()
                if n == 0:
                    species_list = [
                        item.split(".", 1)[0] for item in line.split("\t", 1)[1].split("\t")
                    ]
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
                            genes.append(
                                ", ".join(
                                    [
                                        g
                                        for g in re.split(" |,", gene)
                                        if len(g) != 0
                                    ]
                                ))
                    og = predog_key + ": " + ", ".join(genes)
                    writer.write(og + "\n")


def broccoli_v1(input_file, output_file):

    with open(output_file, "w") as writer:
        with open(input_file, "r") as reader:
            for n, line in enumerate(reader):
                line = line.strip()
                if n == 0:
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
                            genes.append( ", ".join(
                                    [
                                        g
                                        for g in re.split(" |,", gene)
                                        if len(g) != 0
                                    ]
                                ))
                    og = predog_key + ": " + ", ".join(genes)
                    writer.write(og + "\n")

def broccoli_v2(input_file, output_file, protemes_dir: pathlib.Path):
    species_gene_dict = {}

    for file in protemes_dir.iterdir():
        species = file.name.split(".", 1)[0]
        with open(file, "r") as reader:
            for line in reader:
                if ">" in line:
                    gene = line[1:].strip().split(".", 1)[-1]
                    species_gene_dict[gene] = species

    with open(output_file, "w") as writer:
        with open(input_file, "r") as reader:
            for i, line in enumerate(reader):
                line = line.strip()
                if i == 0:
                    continue

                predog_key, genes_str = line.split("\t", 1)
                genes_list = genes_str.split()
                genes = []
                for gene in genes_list:
                    if len(gene) == 0:
                        continue
                    species = species_gene_dict.get(gene)
                    if species is not None:
                        genes.append(species + "." + gene)
                    else:
                        genes.append(gene)
                og = predog_key + ": " + ", ".join(genes)
                writer.write(og + "\n")

def broccoli(input_file, output_file, protemes_dir: Optional[pathlib.Path] = None,):

    with open(input_file) as f:
        header = f.readline().strip('\n')

    if "protein_names" in header:
        if protemes_dir is not None:
            broccoli_v2(input_file, output_file, protemes_dir)
        else:
            print(f"{input_file.name} requires species information,"
                  " please provide the path for the database!")
            return 
    else:
        broccoli_v1(input_file, output_file)


def proteinortho(input_file, output_file):
    predog_base_name = "PredOG%07d"
    with open(output_file, "w") as writer:
        with open(input_file, "r") as reader:
            for n, line in enumerate(reader):
                line = line.strip()
                if n == 0:
                    species_list = [
                        item.split(".", 1)[0] for item in line.split("\t", 3)[-1].split()
                    ]
                else:
                    predog_key = predog_base_name % n
                    genes_list = line.split("\t")[3:]
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

def fastoma(input_file, output_file, protemes_dir: pathlib.Path):
    species_gene_dict = {}

    for file in protemes_dir.iterdir():
        species = file.name.split(".", 1)[0]
        with open(file, "r") as reader:
            for line in reader:
                if ">" in line:
                    gene = line[1:].strip().split(".", 1)[-1]
                    species_gene_dict[gene] = species

    predog_dict = {}
    with open(input_file, "r") as reader:
        for i, line in enumerate(reader):
            line = line.strip()
            if i == 0:
                continue

            predog_key, gene = line.split("\t")[:2]
            predog_key = predog_key.replace(":", "")
            gene = gene.strip()
            species = species_gene_dict.get(gene) 
            if species is not None:
                predog_gene = species + "." + gene
            else:
                predog_gene = gene

            if predog_key not in predog_dict:
                predog_dict[predog_key] = [predog_gene]
            elif predog_key in predog_dict:
                predog_dict[predog_key].append(predog_gene)

    with open(output_file, "w") as writer:      
        for predog_key, genes in predog_dict.items():
            og = predog_key + ": " + ", ".join(genes)
            writer.write(og + "\n")


def swiftortho(input_file, output_file):
    predog_base_name = "PredOG%07d"
    with open(output_file, "w") as writer:
        with open(input_file, "r") as reader:
            for n, line in enumerate(reader):
                line = line.strip()
                predog_key = predog_base_name % n
                genes = []
                for gene in line.split("\t"):
                    if len(gene) == 0:
                        continue
                    species, gene = gene.split("|")
                    if species == gene.split(".", 1)[0]:
                        genes.append(gene)
                    else:
                        genes.append(species + "." + gene)
                    
                og = predog_key + ": " + ", ".join(genes)
                writer.write(og + "\n")


def orthohmm(input_file, output_file, protemes_dir: pathlib.Path):
    species_gene_dict = {}

    for file in protemes_dir.iterdir():
        species = file.name.split(".", 1)[0]
        with open(file, "r") as reader:
            for line in reader:
                if ">" in line:
                    gene = line[1:].strip().split(".", 1)[-1]
                    species_gene_dict[gene] = species

    with open(output_file, "w") as writer:
        with open(input_file, "r") as reader:
            # unassigned_genes_list = []
            for line in reader:
                line = line.strip().split(": ")
                predog_key = line[0]

                genes = []
                for gene in line[1].split():
                    if len(gene) == 0:
                        continue
                    species = species_gene_dict.get(gene) 
                    if species is not None:
                        genes.append(species + "." + gene)
                    else:
                        genes.append(gene)
                if len(genes) > 1:
                    og = predog_key + ": " + ", ".join(genes)
                    writer.write(og + "\n")
            #     else:
            #         unassigned_genes_list .extend(genes)
            # og = "unassigned_genes" + ": " + ", ".join(unassigned_genes_list)
            # writer.write(og + "\n")


                