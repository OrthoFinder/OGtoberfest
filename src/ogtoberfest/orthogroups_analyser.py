import re 
import pathlib
from collections import Counter
from typing import Optional, List, Set, Tuple, Dict 
from ogtoberfest import utils


def get_expected_genes(
        input_dir: pathlib.Path,
        outgroups: Optional[List[str]] = None,
        additional_specieces: Optional[List[str]] = None,
        input_species: Optional[List[str]] = None,
    ) -> Tuple[Set[str], Set[str]]:

    all_genes = set()
    species_of_interest = set()
    for fn in input_dir.iterdir():
        species = fn.name.split(".")[0]
        species = utils.curate_labels(species)

        if outgroups is not None:
            if species.lower() in outgroups:
                continue

        if additional_specieces is not None:
            if species.lower() in additional_specieces:
                continue

        if outgroups is None and additional_specieces is None:
            if input_species is not None:
                if species not in input_species:
                    continue

        species_of_interest.add(species)

        fpath = input_dir / fn
        with open(fpath, "r") as infile:
            for l in infile:
                if l.startswith(">"):
                    all_genes.add(l[1:].rstrip())

    # assert len(all_genes) == n_genes_total, f"species_of_interest {len(species_of_interest)} all_genes {len(all_genes)} vs. n_genes_tatal {n_genes_total}"

    return all_genes


def check_orthobench_orthogroups(ogs, exp_genes):

    expected_gene_names_base = {
        "ENSP000": "Homo sapiens",
        "ENSRNO": "Rattus norvegicus",
        "ENSCAF": "Canis familiaris",
        "ENSMUS": "Mus musculus",
        "ENSMOD": "Monodelphis domestica",
        "ENSGAL": "Gallus gallus",
        "ENSCIN": "Ciona intestinalis",
        "ENSTNI": "Tetraodon nigroviridis",
        "ENSPTR": "Pan troglodytes",
        "ENSDAR": "Danio rerio",
        "FBpp": "Drosophila melanogaster",
        "WBGene": "C elegans",
    }
    n_genes_total = 251378

    all_pred_genes = set([g for og in ogs for g in ogs[og]])
    # all_pred_genes = set([g for og in ogs for g in og])
    x = all_pred_genes.difference(exp_genes)
    if len(x) != 0:
        print(
            "ERROR: found extra genes in input file, check its formatting is correct and there are no incorrect genes"
        )
        print("Examples:")
        for g in list(x)[:10]:
            print(g)
    x = exp_genes.difference(all_pred_genes)
    if len(x) != 0:
        print("Examples of genes not in file:")
        for g in list(x)[:3]:
            print(g)

    n_genes = sum([len(ogs[og]) for og in ogs])
    # n_genes = sum([len(og) for og in ogs])
    if n_genes < 0.5 * n_genes_total:
        print("ERROR: Too many missing genes in predicted orthogroups.")
        print("Orthogroups should contain at least 50% of all genes but")
        print("orthogroups file only contained %d genes" % n_genes)
        exit()

    all_genes = [g for og in ogs for g in ogs[og]]
    # all_genes = [g for og in ogs for g in og]
    n_genes_no_dups = len(set(all_genes))
    if n_genes_no_dups != n_genes:
        print(
            "ERROR: Some genes appear in multiple orthogroups, benchmark are meaningless with such data."
        )
        print("with such data")
        print((n_genes_no_dups, n_genes))
        c = Counter(all_genes)
        for g, n in c.most_common(10):
            print("%d: %s" % (n, g))
        raise Exception()
    # checked genes from each of the expected species are present
    for g_pat, sp in expected_gene_names_base.items():
        if not any(g.split(".", 1)[1].startswith(g_pat) for g in all_genes):
            print("ERROR: No genes found from %s" % sp)
    p = 100.0 * n_genes / float(n_genes_total)
    print("%d genes found in predicted orthogroups, this is %0.1f%% of all genes" % (n_genes, p))


def filter_orthogroups(ref_ogs: Dict[str, Set], pred_ogs: Dict[str, Set]):

    ref_ogs_set = set()
    for og in ref_ogs:
        ref_ogs_set = ref_ogs_set | set(ref_ogs[og])

    for pog in pred_ogs:
        diff = pred_ogs[pog] - ref_ogs_set
        if len(diff) > 0:
            pred_ogs[pog].difference_update(diff)

    return pred_ogs
