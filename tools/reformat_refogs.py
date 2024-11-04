import pathlib 
import os 
import numpy as np
import pandas as pd

protemes_dir = r"./orthobench_proteomes"
refog_dir = r"./RefOGs"
reformated_refog = r"./OrthoBench_RefOGs.tsv"


species_gene_dict = {}
species_gene_length = {}
accession = ""
seq = ""
for file in pathlib.Path(protemes_dir).iterdir():
    species = os.path.basename(file).rsplit(".")[0]
    with open(file, "r") as reader:
        for line in reader:
            line = line.replace("\n", "").strip()
            if line.startswith(">") or ">" in line:
                if accession:
                    species_gene_length[accession] = len(seq)

                gene = line[1:].replace('>', '')
                accession = species + "." + gene
                species_gene_dict[gene] = accession
                seq = ""
            else:
                seq += line
        if accession:
            species_gene_length[accession] = len(seq)


lca_refogs_dict = {}
refogs_dict = {}
count_lca = 0
total_refog = 0
for refog_path in pathlib.Path(refog_dir).iterdir():
    if refog_path.is_dir():
        for lca_file in refog_path.iterdir():
            lca_refog = lca_file.name.rsplit(".", 1)[0]
            lca_refogs_dict[lca_refog] = set()
            with open(lca_file) as reader1:
                for line in reader1:
                    gene = line.strip()
                    if len(gene) != 0:
                        lca_refogs_dict[lca_refog].add(species_gene_dict[gene])
            count_lca += len(lca_refogs_dict[lca_refog])
    else:
        refog_key = refog_path.name.rsplit(".", 1)[0]
        refogs_dict[refog_key] = set()
        
        with open(refog_path) as reader2:
            for line in reader2:
                gene = line.strip()
                if len(gene) != 0:
                    refogs_dict[refog_key].add(species_gene_dict[gene])
                    # species_count_refogs_dict[refog_key].add(species_gene_dict[gene])
        total_refog += len(refogs_dict[refog_key])

print(count_lca, total_refog)

updated_refogs_dict = {}

for refog_key, refog in refogs_dict.items(): 
    lca_ogs = lca_refogs_dict.get(refog_key)
    if lca_ogs is not None:
        updated_refogs_dict[refog_key] = refog.difference(lca_ogs)
    else:
        updated_refogs_dict[refog_key] = refog

colnames = ["num_species", "num_genes", 
            "min_length", "max_length", "mean_length", 
            "median_length", "Orthogroups"]

refog_stats_dict = {}

for refog_key, refog in updated_refogs_dict.items(): 
    num_species = len(
        set(
            [gene.split(".", 1)[0] for gene in refog]
        )
    )
    num_genes = len(refog)

    min_length = np.amin(
        [species_gene_length[gene] for gene in refog]
    )

    max_length = np.amax(
        [species_gene_length[gene] for gene in refog]
    )

    mean_length = np.mean(
        [species_gene_length[gene] for gene in refog]
    ).round().astype(int)

    median_length = np.median(
        [species_gene_length[gene] for gene in refog]
    ).round().astype(int)

    ogs_string = ", ".join(refog)

    refog_stats_dict[refog_key] = (
        num_species,
        num_genes,
        min_length,
        max_length,
        mean_length,
        median_length,
        ogs_string
    )
    
df = pd.DataFrame.from_dict(refog_stats_dict, orient='index', columns=colnames)
df.reset_index(inplace=True)
df.rename(columns={"index": "RefOGs"}, inplace=True)
df.to_csv(reformated_refog, sep="\t", index=False)
