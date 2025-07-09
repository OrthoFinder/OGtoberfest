orig_refog_path = r"./Sim1k/reference_orthogroups_orig.txt"
orig_refog_stats_path = r"./Sim1k/reference_orthogroups_stats_orig.txt"

refog_path = r"./Sim1k/reference_orthogroups.txt"
refog_stats_path = r"./Sim1k/reference_orthogroups_stats.txt"
refog_tree_path = r"./Sim1k/reference_orthogroups_tree.txt"


with open(refog_stats_path, "w") as writer:
    with open(orig_refog_stats_path) as reader:
        for i, line in enumerate(reader):
            if i == 0:
                line.replace("Orthogroup", "RefOGs")
            line = line.strip().split("\t", 1)
            refog_key, refog_stats = line[0], line[1].rsplit("\t", 1)[0]
            writer.write(refog_key + "\t" + refog_stats + "\n")


with open(refog_path, "w") as writer:
    with open(orig_refog_path) as reader:
        for i, line in enumerate(reader):
            if i == 0:
                continue
            line = line.strip().split("\t")
            refog_key, genes = line[0], line[1:]
            
            genes = [
                ", ".join(gene.split(","))
                for gene in genes
                if len(gene) !=0
            ]
            genes_str = ", ".join(genes)
            writer.write(refog_key + ": " + genes_str + "\n")



with open(refog_tree_path, "w") as writer:
    with open(orig_refog_stats_path) as reader:
        for i, line in enumerate(reader):
            if i == 0:
                continue
            line = line.strip().split("\t", 1)
            refog_key, tree = line[0], line[1].rsplit("\t", 1)[1]
            writer.write(refog_key + "\t" + tree + "\n")