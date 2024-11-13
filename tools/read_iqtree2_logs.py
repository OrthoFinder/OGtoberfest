import pathlib 

log_path = r"./OrthoBench/iqtree2_logs"

invariable_site_dict = {}
gamma_shape = {}
total_tree_length = {}
for file in pathlib.Path(log_path).iterdir():
    refog_key = file.name.split(".")[0]
    with open(file) as reader:
        for line in reader:
            if "Proportion of invariable sites" in line:
                invariable_site_dict[refog_key] = line.strip().split(": ")[-1]

            if "Gamma shape alpha" in line:
                gamma_shape[refog_key] = line.strip().split(": ")[-1]

            if "Total tree length (sum of branch lengths)" in line:
                total_tree_length[refog_key] = line.strip().split(": ")[-1]

print(invariable_site_dict)
print(gamma_shape)
print(total_tree_length)
