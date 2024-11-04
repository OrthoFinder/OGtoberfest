import pathlib

input_path = r"./OrthoBench/proteomes"
output_path = r"./OrthoBench/proteomes_gene_labels"

for file in pathlib.Path(input_path).iterdir():
    output_file_path = pathlib.Path(output_path) / file.name
    with open(output_file_path, "w") as writer:
        with open(file, "r") as reader:

            species = file.name.split(".", 1)[0]
            print(species)
            for line in reader:
                if line.startswith(">"):
                    gene = ">" + species + "." + line[1:]
                    writer.write(gene)
