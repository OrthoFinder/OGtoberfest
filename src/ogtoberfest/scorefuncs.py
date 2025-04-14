import copy
from itertools import product

import numpy as np
from scipy.special import comb, factorial, gammaln
from scipy.stats import expon
from sklearn.metrics import cluster
from typing import Dict, List, Optional, Tuple, Set


def sparse_contingency_table_count(U, V, n):
    n_fac_log = gammaln(n + 1)

    ar_fac_log = np.sum([gammaln(len(u_set) + 1) for _, u_set in U.items()])
    bs_fac_log = np.sum([gammaln(len(v_set) + 1) for _, v_set in V.items()])

    ar_comb = np.sum([comb(len(u_set), 2) for _, u_set in U.items()])
    bs_comb = np.sum([comb(len(v_set), 2) for _, v_set in V.items()])

    log_sigma_ab = np.log2(np.exp(n_fac_log - ar_fac_log - bs_fac_log)) + 2 * ar_comb * bs_comb / n**2

    return log_sigma_ab


def dense_contingency_table_count(U, V, n):

    R = len(U)
    S = len(V)
    RS = R * S
    w = n / (n + RS / 2)

    xrs = np.array([(1 - w) / R + w * len(u_set) / n for _, u_set in U.items()])
    yss = np.array([(1 - w) / S + w * len(v_set) / n for _, v_set in V.items()])

    epsilon = 1e-10

    mu = max((R + 1) / (R * np.sum(yss**2)) - 1 / R, epsilon)
    v = max((S + 1) / (S * np.sum(xrs**2)) - 1 / S, epsilon)

    log_factorial_muR_1 = gammaln(mu * R + 1)
    log_factorial_vS_1 = gammaln(v * S + 1)
    log_factorial_v_1 = gammaln(v + 1)
    log_factorial_R_1 = gammaln(R + 1)
    log_factorial_mu_1 = gammaln(mu + 1)
    log_factorial_S_1 = gammaln(S + 1)

    log_factorial_term = (
        log_factorial_muR_1
        - log_factorial_vS_1
        - S * (log_factorial_v_1 - log_factorial_R_1)
        - R * (log_factorial_mu_1 - log_factorial_S_1)
    )

    log_sigma_ab = (
        (R - 1) * (S - 1) * np.log2(n + RS / 2)
        + (R + v - 2) * np.sum(np.log2(yss)) / 2
        + (S + mu - 2) * np.sum(np.log2(xrs)) / 2
        + log_factorial_term / 2
    )

    return log_sigma_ab


def normalised_reduced_mutual_information(U, V, n, list_products, rmi_type="dense"):
    """
    reference:
    [1] Improved mutual information measure for clustering, classification, and community detection
    [2] Comment on "Improved mutual information measure for clustering, classification, and community detection"

    """

    n_fac_log = gammaln(n + 1)

    ar_fac_log = np.sum([gammaln(len(u_set) + 1) for _, u_set in U.items()])
    bs_fac_log = np.sum([gammaln(len(v_set) + 1) for _, v_set in V.items()])

    crs_fac_log = 0
    for (u, u_set), (v, v_set) in list_products:
        crs = len(u_set & v_set)
        if crs != 0:
            crs_fac_log += gammaln(crs + 1)

    if rmi_type == "dense":
        log_sigma_ab = dense_contingency_table_count(U, V, n)
        log_sigma_aa = dense_contingency_table_count(U, U, n)
        log_sigma_bb = dense_contingency_table_count(V, V, n)
    elif rmi_type == "sparse":
        log_sigma_ab = dense_contingency_table_count(U, V, n)
        log_sigma_aa = dense_contingency_table_count(U, U, n)
        log_sigma_bb = dense_contingency_table_count(V, V, n)

    log_numerator = 2 * (np.log2(np.exp(n_fac_log + crs_fac_log - ar_fac_log - bs_fac_log)) - log_sigma_ab)
    log_denominator = (
        np.log2(np.exp(n_fac_log - ar_fac_log)) + np.log2(np.exp(n_fac_log - bs_fac_log)) - log_sigma_aa - log_sigma_bb
    )

    NRMI = log_numerator / log_denominator

    return NRMI


def mutual_information(U, V, n, list_products):
    I = 0.0

    for (u, u_set), (v, v_set) in list_products:
        nij = len(u_set & v_set)
        ai = len(U[u])
        bj = len(V[v])
        pi = ai / n
        pj = bj / n
        pij = nij / n
        # print(u, v, nij, ai, bj)
        if nij == ai == bj == n:
            I += 1.
            return I
        
        if pij != 0:
            I += pij * np.log2(pij / (pi * pj))

    return I


def homogeneity(U, V, n, list_products):
    """
    Measures whether each cluster contains only members of a single class.
    """
    I = mutual_information(U, V, n, list_products)
    HU = np.sum([-(len(u_set) / n) * np.log2(len(u_set) / n) for _, u_set in U.items()])
    h = I / HU
    return h


def completeness(U, V, n, list_products):
    """
    Measures whether all members of a given class are assigned to the same cluster
    """
    I = mutual_information(U, V, n, list_products)
    HV = np.sum([-(len(v_set) / n) * np.log2(len(v_set) / n) for _, v_set in V.items()])
    c = I / HV
    return c


def v_measure(U, V, n, list_products):
    """
    The harmonic mean of homogeneity and completeness,
    providing a single score to assess the clustering performance

    Pros:

        1. Provide a direct assessment of the match between cluster assignments
           and the class labels in terms of both purity and class representation.
        2. The scores are bounded between 0 and 1.
        3. Have intuitive interpretation.
        4. No assumption is made on the cluster structure.

    Cons:
        1. Do not take into account how data points are distributed within each cluster.
        2. Not normalized against random grouping (unlike ARI).
            This means that depending on the number of samples, clusters and classes,
            a completely random grouping of the samples will not always yield the same values for
            homogeneity, completeness, and v-measure. For this reason, for small datasets
            (number of samples < 1000) or large number of clusters (> 10) it is safer to use ARI.

    """
    h = homogeneity(U, V, n, list_products)
    c = completeness(U, V, n, list_products)
    v = 2 * h * c / (h + c)
    return v


def expected_information(U, V, n, list_products):
    EI = 0.0
    for (u, u_set), (v, v_set) in list_products:
        nij = len(u_set & v_set)
        ai = len(U[u])
        bj = len(V[v])
        
        for nij in range(np.amax((0, ai + bj - n)), np.amin((ai, bj)) + 1):
            if nij != 0:
                first_part = (nij / n) * np.log2((n * nij) / (ai * bj))

                second_part_numerator_log = (
                    gammaln(ai + 1) + gammaln(bj + 1) + gammaln(n - ai + 1) + gammaln(n - bj + 1)
                )
                second_part_denomenator_log = (
                    gammaln(n + 1)
                    + gammaln(nij + 1)
                    + gammaln(ai - nij + 1)
                    + gammaln(bj - nij + 1)
                    + gammaln(n - ai - bj + nij + 1)
                )

                EI += first_part * np.exp(second_part_numerator_log - second_part_denomenator_log)
    return EI 


def adjusted_mutual_information(U, V, n, list_products):
    """
    AMI:
        Range: [-1, 1]
            1: Perfect agreement between the two clusterings.
            0: Agreement is no better than random chance.
            Negative values: Worse than random chance.

    AMI should be used when the reference clustering is unbalanced and there exist small clusters.
    reference:
    [1] Information Theoretic Measures for Clusterings Comparison: Is a Correction for Chance Necessary?
    [2] Adjusting for Chance Clustering Comparison Measures
    """

    MI = mutual_information(U, V, n, list_products)
    EI = expected_information(U, V, n, list_products)
    HU = np.sum([-(len(u_set) / n) * np.log2(len(u_set) / n) for _, u_set in U.items()])
    HV = np.sum([-(len(v_set) / n) * np.log2(len(v_set) / n) for _, v_set in V.items()])
    if MI == 1.0:
        return 1.0
    AMI = (MI - EI) / (np.mean([HU, HV]) - EI)
    return AMI


def variation_of_information(U, V, n, list_products, vi_type="normalised"):

    I = mutual_information(U, V, n, list_products)
    EI = expected_information(U, V, n, list_products)

    HU = np.sum([-(len(u_set) / n) * np.log2(len(u_set) / n) for _, u_set in U.items()])
    HV = np.sum([-(len(v_set) / n) * np.log2(len(v_set) / n) for _, v_set in V.items()])

    if I == 1.0 or (HU + HV == 2 * EI):
        return 0.0
    
    if vi_type == "normalised":
        NVI = (HU + HV - 2 * I) / (HU + HV - 2 * EI)
        return max(0.0, NVI)
    elif vi_type == "adjusted":
        AVI = (2 * I - 2 * EI) / (HU + HV - 2 * EI)
        return max(0.0, AVI)


# def sparse_contingency_table_count(U, V, n):

#     n_fac = factorial(n, exact=False)
#     ar_fac = np.prod([factorial(len(u_set)) for _, u_set in U.items()])
#     bs_fac = np.prod([factorial(len(v_set)) for _, v_set in V.items()])

#     ar_comb = comb([len(u_set) for _, u_set in U.items()], 2).sum()
#     bs_comb = comb([len(v_set) for _, v_set in V.items()], 2).sum()
#     log_sigma_ab = np.log2(n_fac / (ar_fac * bs_fac)) + 2 * ar_comb * bs_comb / n**2 
#     return log_sigma_ab

# def dense_contingency_table_count(U, V, n):

#     R = len(U)
#     S = len(V)
#     RS = R * S
#     w = n / (n + RS / 2)
#     xrs = np.array([(1 - w) / R + w * len(u_set) / n for _, u_set in U.items()])
#     yss = np.array([(1 - w) / S + w * len(v_set) / n for _, v_set in V.items()])
#     mu = (R + 1) / (R * np.sum(yss**2)) - 1 / R 
#     v = (S + 1) / (S * np.sum(xrs**2)) - 1 / S
#     log_sigma_ab = (R - 1) * (S -1) * np.log2(n + RS / 2) + \
#             (R + v - 2) * np.sum(np.log2(yss)) / 2 + \
#             (S + mu - 2) * np.sum(np.log2(xrs)) / 2 + \
#             np.log2(factorial(mu * R - 1) * factorial(v * S - 1) / ((factorial(v - 1) * factorial(R - 1))**S * (factorial(mu - 1) * factorial(S - 1))**R)) / 2

#     return log_sigma_ab

# def normalised_reduced_mutual_information(U, V, n, list_products, rmi_type="dense"):

#     """
#     reference:
#     [1] Improved mutual information measure for clustering, classification, and community detection
#     [2] Comment on "Improved mutual information measure for clustering, classification, and community detection"
    
#     """

#     n_fac = factorial(n, exact=False)
#     ar_fac = np.prod([factorial(len(u_set)) for _, u_set in U.items()])
#     bs_fac = np.prod([factorial(len(v_set)) for _, v_set in V.items()])
    
#     crs_fac = 1
#     for (u, u_set), (v, v_set) in list_products:
#         crs = len(u_set & v_set)
#         if crs != 0:
#             crs_fac *= factorial(crs)
#     if rmi_type == "dense":
#         log_sigma_ab = dense_contingency_table_count(U, V, n)
#         log_sigma_aa = dense_contingency_table_count(U, U, n)
#         log_sigma_bb = dense_contingency_table_count(V, V, n)

#     elif rmi_type == "sparse":
#         log_sigma_ab = dense_contingency_table_count(U, V, n)
#         log_sigma_aa = dense_contingency_table_count(U, U, n)
#         log_sigma_bb = dense_contingency_table_count(V, V, n)

#     # RMI = np.log2(n_fac * crs_fac / (ar_fac * bs_fac)) / n - log_sigma_ab / n 
#     numerator = 2 * (np.log2(n_fac * crs_fac / (ar_fac * bs_fac)) - log_sigma_ab)
#     denomenator = np.log2(factorial(n) / ar_fac) + np.log2(factorial(n) / bs_fac) - log_sigma_aa - log_sigma_bb
#     NRMI = numerator / denomenator

#     return NRMI

# def mutual_information(U, V, n, list_products):
#     I = 0.0
#     for (u, u_set), (v, v_set) in list_products:

#         nij = len(u_set & v_set)
#         ai = len(U[u])
#         bj = len(V[v])
#         pi = ai / n 
#         pj = bj / n 
#         pij = nij / n
#         if pij != 0:
#             I += pij * np.log2(pij / (pi * pj))
#     return I

# def homogeneity(U, V, n, list_products):
#     """
#     Measures whether each cluster contains only members of a single class.
#     """
#     I = mutual_information(U, V, n, list_products)
#     HU = np.sum([-(len(u_set) / n) * np.log2(len(u_set) / n) for _, u_set in U.items()])
#     h = I / HU 
#     return h 

# def completeness(U, V, n, list_products):
#     """
#     Measures whether all members of a given class are assigned to the same cluster
#     """
#     I = mutual_information(U, V, n, list_products)
#     HV = np.sum([-(len(v_set) / n) * np.log2(len(v_set) / n) for _, v_set in V.items()])
#     c = I / HV 
#     return c

# def v_measure(U, V, n, list_products):
#     """
#     The harmonic mean of homogeneity and completeness, 
#     providing a single score to assess the clustering performance

#     Pros:

#         1. Provide a direct assessment of the match between cluster assignments 
#            and the class labels in terms of both purity and class representation.
#         2. The scores are bounded between 0 and 1.
#         3. Have intuitive interpretation.
#         4. No assumption is made on the cluster structure.

#     Cons:
#         1. Do not take into account how data points are distributed within each cluster.
#         2. Not normalized against random grouping (unlike ARI). 
#             This means that depending on the number of samples, clusters and classes, 
#             a completely random grouping of the samples will not always yield the same values for 
#             homogeneity, completeness, and v-measure. For this reason, for small datasets 
#             (number of samples < 1000) or large number of clusters (> 10) it is safer to use ARI.

#     """
#     h = homogeneity(U, V, n, list_products)
#     c = completeness(U, V, n, list_products)
#     v = 2 * h * c / (h + c)
#     return v

# def adjusted_mutual_information(U, V, n, list_products):

#     """
#     AMI:
#         Range: [-1, 1]
#             1: Perfect agreement between the two clusterings.
#             0: Agreement is no better than random chance.
#             Negative values: Worse than random chance.

#     AMI should be used when the reference clustering is unbalanced and there exist small clusters.
#     reference: 
#     [1] Information Theoretic Measures for Clusterings Comparison: Is a Correction for Chance Necessary?
#     [2] Adjusting for Chance Clustering Comparison Measures
#     """
    
#     I = mutual_information(U, V, n, list_products)
#     EI = 0.0
#     for (u, u_set), (v, v_set) in list_products:
#         nij = len(u_set & v_set)
#         ai = len(U[u])
#         bj = len(V[v])
#         for nij in range(np.amax((0, ai + bj - n)), np.amin((ai, bj)) + 1):
#             if nij != 0:
#                 first_part = (nij / n) * np.log2((n * nij) / (ai * bj))
            
#                 second_part_numerator = (factorial(ai) * factorial(bj) * factorial(n - ai) * factorial(n - bj))
#                 second_part_denomenator = (factorial(n) * factorial(nij) * factorial(ai - nij) * factorial(bj - nij) * factorial(n - ai - bj + nij))

#                 EI +=  first_part * second_part_numerator / second_part_denomenator
    
#     HU = np.sum([-(len(u_set) / n) * np.log2(len(u_set) / n) for _, u_set in U.items()])
#     HV = np.sum([-(len(v_set) / n) * np.log2(len(v_set) / n) for _, v_set in V.items()])
#     AMI = (I - EI) / (np.mean([HU, HV]) - EI)

#     return AMI

# def variation_of_information(U, V, n, list_products, vi_type = "normalised"):
#     HU = np.sum([-(len(u_set) / n) * np.log2(len(u_set) / n) for _, u_set in U.items()])
#     HV = np.sum([-(len(v_set) / n) * np.log2(len(v_set) / n) for _, v_set in V.items()])
#     I = mutual_information(U, V, n, list_products)
#     EI = adjusted_mutual_information(U, V, n, list_products)
#     if vi_type == "normalised":
#         NVI  = (HU + HV - 2 * I) / (HU + HV - 2 * EI)
#         return NVI
#     elif vi_type == "adjusted":
#         AVI = (2 * I - 2 * EI) / (HU + HV - 2 * EI)
#         return AVI



def rearange(label_dict):
    item_dict = {}
    for label, item_set in label_dict.items():
        for item in item_set:
            if item not in item_dict:
                item_dict[item] = label
            else:
                print(f"duplicated item found in {label}")

    return item_dict


def rand_index(U, V, n, list_products, ri_type="adjusted"):
    """
    RI:
        Measures the similarity between the cluster assignments and the true class labels by making pairwise comparisons.

        The problem with Rand Index is that it can yield high values even for random cluster assignments,
        particularly when the number of clusters is large.
        This is because when the number of clusters increases,
        the probability of randomly assigning points with different labels to different clusters increases.
        Consequently, a specific RI value can be ambiguous,
        as it is not clear how much of the score is due to chance versus actual agreement.

    ARI:
        Taking into account the expected RI score if the cluster assignments were made randomly.

        Range: [-1, 1]
            1: Perfect agreement between the two clusterings (the clusterings are identical).
            0: The level of agreement is what would be expected by random chance.
            Negative values: Indicates worse than random clustering (the two clusterings disagree more than would be expected by chance).

        ARI should be used when the reference clustering has large equal sized clusters.

    Pros:

        1. RI scores are bounded between 0 and 1, and ARI scores are bounded between -1 and 1.
           The bounded range makes it easy to compare the scores between different algorithms.
        2. Random (uniform) cluster assignments have an ARI score close to 0 for any number of samples and clusters.
        3. The adjustment of ARI for chance makes it more reliable and interpretable.
        4. No assumption is made on the cluster structure,
           which makes these metrics useful for comparing different clustering algorithms independent of the cluster shapes.
    Cons:
        1. Do not take into account how data points are distributed within each cluster.
           For example, even if the data points within a cluster are spread out or form sub-clusters, this will not affect the RI or ARI scores.

    """

    a_list = [len(u_set & v_set) for (u, u_set), (v, v_set) in list_products]
    a = sum(comb(x, 2, exact=False) for x in a_list)  
    marginal_u = [len(u_set) for u_set in U.values()]
    marginal_v = [len(v_set) for v_set in V.values()]

    marginal_u_comb = sum(comb(x, 2, exact=False) for x in marginal_u)
    marginal_v_comb = sum(comb(x, 2, exact=False) for x in marginal_v)
    total_pairs = comb(n, 2, exact=False)

    b = marginal_u_comb - a
    c = marginal_v_comb - a
    d = total_pairs - a - b - c

    if ri_type != "adjusted":
        RI = (a + d) / total_pairs
        return RI

    expected_index = marginal_u_comb * marginal_v_comb / total_pairs
    maximum_index = (marginal_u_comb + marginal_v_comb) / 2

    if np.isclose(maximum_index, expected_index):  
        ARI = 1.0
    else:
        ARI = (a - expected_index) / (maximum_index - expected_index)

    return ARI

def kl_divergence(P, Q, n, m):

    kld = np.nansum(
        [
            (len(refog) / n) * np.log2((len(refog) / n) / (Q[refog_key] / m)) if Q[refog_key] != 0 else np.nan
            for refog_key, refog in P.items()
        ]
    )

    return kld


def check_missing_genes(U, V_prime, num_round):

    missing_genes_dict = {}
    missing_genes_count_dict = {}
    missing_genes_per_dict = {}
    total_missing_genes_set = set()
    for refog_key, refog in U.items():
        predog_set = set()
        for predog in V_prime[refog_key].values():
            predog_set = predog_set.union(predog)

        nrefog = len(refog)
        FN = refog - predog_set
        for item in FN:
            for predog in V_prime[refog_key].values():
                if item in predog:
                    print(item)
        if len(FN) != 0:
            total_missing_genes_set = total_missing_genes_set.union(FN)
            missing_genes_dict[refog_key] = FN
            missing_genes_count_dict[refog_key] = len(FN)
            missing_genes_per_dict[refog_key] = np.round(100. * len(FN) / nrefog, num_round)
        else:
            missing_genes_dict[refog_key] = set()
            missing_genes_count_dict[refog_key] = 0
            missing_genes_per_dict[refog_key] = 0.0

    total_missing_genes = len(total_missing_genes_set)

    return (
        total_missing_genes,
        total_missing_genes_set,
        missing_genes_dict,
        missing_genes_count_dict,
        missing_genes_per_dict,
    )


def check_misplaced_genes(U, V_prime, n):
    missing_genes_dict = {}
    for refog_key, refog in U.items():
        predog_set = set()
        for predog in V_prime[refog_key].values():
            predog_set = predog_set.union(predog)

        FN = refog - predog_set
        if len(FN) != 0:
            missing_genes_dict[refog_key] = FN

    misplaced_genes_dict = {}
    total_misplaced_genes = 0
    for refog_key, missing_genes in missing_genes_dict.items():
        if len(missing_genes) != 0:
            misplaced_genes_dict[refog_key] = [
                (predog_key, missing_gene)
                for missing_gene in missing_genes
                for predog_dict in V_prime.values()
                for predog_key, predog in predog_dict.items()
                if missing_gene in predog
            ]
            total_misplaced_genes += len(misplaced_genes_dict[refog_key])

    total_misplaced_genes_proportion = total_misplaced_genes / n

    return total_misplaced_genes, total_misplaced_genes_proportion, misplaced_genes_dict


def check_total_missing_genes(
    total_misplaced_genes, missing_genes_count_dict, missing_genes_dict, misplaced_genes_dict, n
):
    total_missing_genes_dict = {}
    total_missing_genes_set = set()
    for refog_key, missing_genes in missing_genes_dict.items():
        misplaced_genes = misplaced_genes_dict.get(refog_key)
        if misplaced_genes is not None:
            for _, gene in misplaced_genes:
                missing_genes.discard(gene)
        total_missing_genes_dict[refog_key] = missing_genes
        total_missing_genes_set = total_missing_genes_set.union(missing_genes)

    total_missing_genes = np.sum([*missing_genes_count_dict.values()]) - total_misplaced_genes
    total_missing_genes_proportion = total_missing_genes / n

    return total_missing_genes, total_missing_genes_proportion, total_missing_genes_dict, total_missing_genes_set


def entropy_score(U, V_prime, n):

    entropy_dict = {}
    for refog_key, refog in U.items():

        if len(V_prime[refog_key]) != 0:
            overlap_list = [len(refog & predog) for predog in V_prime[refog_key].values()]
            predog_set = set()
            for predog in V_prime[refog_key].values():
                predog_set = predog_set.union(predog)
            missing_genes_count = len(refog - predog_set)

            if missing_genes_count != 0:
                overlap_list.append(missing_genes_count)
            nrefog = len(refog)
            overlap_arr = np.array(overlap_list)
            probability = overlap_arr / nrefog
            max_entropy = np.log2(nrefog)
            entropy_dict[refog_key] = np.sum(-probability * np.log2(probability)) / max_entropy
        else:
            entropy_dict[refog_key] = 1 # np.nan

    total_entropy = np.nansum([entropy_dict[refog_key] * len(refog) / n for refog_key, refog in U.items()])

    return total_entropy, entropy_dict


def macro_and_weighted_avg_scores(U, V_prime, n, num_round):

    weighted_recall_dict = {}
    weighted_precision_dict = {}
    weighted_f1score_dict = {}
    weighted_fowlkes_mallows_dict = {}
    weighted_JI_dict = {}
    weighted_dissimilarity_dict = {}
    weighted_distance_dict = {}

    effective_size_precision_weighted_dict = {}
    effective_size_JI_weighted_dict = {}
    effective_size_JI_refog_weighted_dict = {}

    all_score_dict = {}
    all_TP_dict = {}

    for refog_key, refog in U.items():
        if len(V_prime[refog_key]) != 0:
            predog_keys, predog_tuple = zip(*V_prime[refog_key].items())
            npredog_arr = np.array([len(pog) for pog in predog_tuple])

            overlap_arr = np.array([len(refog & predog) for predog in V_prime[refog_key].values()])
            union_arr = np.array([len(refog | predog) for predog in V_prime[refog_key].values()])
            nrefog = len(refog)
            total_intersection = np.sum(overlap_arr)
            predog_weights_arr = overlap_arr / total_intersection

            predog_recall_arr = overlap_arr / nrefog
            predog_weighted_recall_arr = predog_recall_arr * predog_weights_arr
            weighted_recall_dict[refog_key] = np.sum(predog_weighted_recall_arr)

            predog_precision_arr = np.array(
                [len(refog & predog) / len(predog) for predog in V_prime[refog_key].values()]
            )
            predog_weighted_precision_arr = predog_precision_arr * predog_weights_arr
            weighted_precision_dict[refog_key] = np.sum(predog_weighted_precision_arr)

            predog_f1score_arr = (
                2 * predog_recall_arr * predog_precision_arr / (predog_recall_arr + predog_precision_arr)
            )
            predog_weighted_f1score_arr = predog_f1score_arr * predog_weights_arr
            weighted_f1score_dict[refog_key] = np.sum(predog_weighted_f1score_arr)

            predog_fowlkes_mallows_arr = np.sqrt(predog_recall_arr * predog_precision_arr)
            predog_weighted_fowlkes_mallows_arr = predog_fowlkes_mallows_arr * predog_weights_arr
            weighted_fowlkes_mallows_dict[refog_key] = np.sum(predog_weighted_fowlkes_mallows_arr)

            predog_JI_arr = overlap_arr / union_arr
            predog_weighted_JI_arr = predog_JI_arr * predog_weights_arr
            weighted_JI_dict[refog_key] = np.sum(predog_weighted_JI_arr)

            predog_dissimilarity_arr = 1 - overlap_arr / union_arr
            predog_weighted_dissimilarity_arr = predog_dissimilarity_arr * predog_weights_arr
            weighted_dissimilarity_dict[refog_key] = np.sum(predog_weighted_dissimilarity_arr)

            predog_distances_arr = expon.ppf(predog_dissimilarity_arr, scale=len(V_prime[refog_key]))
            weighted_distance_dict[refog_key] = np.sqrt(np.mean(predog_distances_arr**2))

            effective_size_precision_weighted_arr = overlap_arr * predog_weighted_precision_arr
            effective_size_precision_weighted_dict[refog_key] = \
                int(np.sum(effective_size_precision_weighted_arr))

            effective_size_JI_weighted_arr = npredog_arr * predog_JI_arr
            effective_size_JI_weighted_dict[refog_key] = \
                int(np.sum(effective_size_JI_weighted_arr))

            effective_size_JI_refog_weighted_dict[refog_key] = \
                int(nrefog * weighted_JI_dict[refog_key])

            all_scores = [
                (
                    predog_key,
                    str(npredog),
                    str(overlap),
                    str(np.round(100 * recall, num_round)),
                    str(np.round(100 * precision, num_round)),
                    str(np.round(100 * fm, num_round)),
                    str(np.round(100 * ji, num_round)),
                    str(np.round(100 * dissimalrity, num_round)),
                    str(np.round(distance, num_round)),
                )
                for predog_key, npredog, overlap, recall, precision, fm, ji, dissimalrity, distance in zip(
                    predog_keys,
                    npredog_arr,
                    overlap_arr,
                    predog_recall_arr,
                    predog_precision_arr,
                    predog_fowlkes_mallows_arr,
                    predog_JI_arr,
                    predog_dissimilarity_arr,
                    predog_distances_arr,
                )
            ]

            all_scores = sorted(all_scores, key=lambda x: x[2])
            all_score_dict[refog_key] = ", ".join([":".join(item) for item in all_scores])
            all_TP_dict[refog_key] = {
                predog_key: refog & predog for predog_key, predog in V_prime[refog_key].items()
            }
        else:
            weighted_recall_dict[refog_key] = 0.0
            weighted_precision_dict[refog_key] = 0.0
            weighted_f1score_dict[refog_key] = 0.0
            weighted_fowlkes_mallows_dict[refog_key] = 0.0
            weighted_JI_dict[refog_key] = 0.0
            weighted_dissimilarity_dict[refog_key] = 0.0
            weighted_distance_dict[refog_key] = 0.0
            effective_size_precision_weighted_dict[refog_key] = 0
            effective_size_JI_weighted_dict[refog_key] = 0
            effective_size_JI_refog_weighted_dict[refog_key] = 0
            all_score_dict[refog_key] = []
            all_TP_dict[refog_key] = {}

    total_weighted_recall = np.sum(
        [weighted_recall_dict[refog_key] * len(refog) / n for refog_key, refog in U.items()]
    )
    total_weighted_precision = np.sum(
        [weighted_precision_dict[refog_key] * len(refog) / n for refog_key, refog in U.items()]
    )
    total_weighted_f1score = np.sum(
        [weighted_f1score_dict[refog_key] * len(refog) / n for refog_key, refog in U.items()]
    )
    total_weighted_fowlkes_mallows = np.sum(
        [weighted_fowlkes_mallows_dict[refog_key] * len(refog) / n for refog_key, refog in U.items()]
    )
    total_weighted_JI = np.sum([weighted_JI_dict[refog_key] * len(refog) / n for refog_key, refog in U.items()])
    total_weighted_dissimilarity = np.sum(
        [weighted_dissimilarity_dict[refog_key] * len(refog) / n for refog_key, refog in U.items()]
    )
    total_weighted_distance = np.sum(
        [weighted_distance_dict[refog_key] * len(refog) / n for refog_key, refog in U.items()]
    )

    macro_recall = np.mean([*weighted_recall_dict.values()])
    macro_precision = np.mean([*weighted_precision_dict.values()])
    macro_f1score = np.mean([*weighted_f1score_dict.values()])
    macro_fowlkes_mallows = np.mean([*weighted_fowlkes_mallows_dict.values()])
    macro_JI = np.mean([*weighted_JI_dict.values()])
    macro_dissimilarity = np.mean([*weighted_dissimilarity_dict.values()])
    macro_distance = np.mean([*weighted_distance_dict.values()])
    macro_dissimilarity_distance = expon.ppf(macro_dissimilarity, scale=len(U))

    output = (
        total_weighted_recall,
        macro_recall,
        weighted_recall_dict,
        total_weighted_precision,
        macro_precision,
        weighted_precision_dict,
        total_weighted_f1score,
        macro_f1score,
        weighted_f1score_dict,
        total_weighted_fowlkes_mallows,
        macro_fowlkes_mallows,
        weighted_fowlkes_mallows_dict,
        total_weighted_JI,
        macro_JI,
        weighted_JI_dict,
        total_weighted_dissimilarity,
        macro_dissimilarity,
        weighted_dissimilarity_dict,
        total_weighted_distance,
        macro_distance,
        weighted_distance_dict,
        macro_dissimilarity_distance,
        all_score_dict,
        all_TP_dict,
        effective_size_precision_weighted_dict,
        effective_size_JI_weighted_dict,
        effective_size_JI_refog_weighted_dict,
    )

    return output


def missing_predogs(V_prime):
    missing_predogs_list = [refog_key for refog_key, predog_dict in V_prime.items() if len(predog_dict) == 0]

    return missing_predogs_list


def V_raw_to_V_prime(U, V_raw):
    """
    V_raw: Dict[predog_key], predog_set]
        consists of all the predogs, it can contain genes that don't belong to any of the refogs.

    V_prime: Dict[refog_key, Dict[predog_key], predog_set]]
        a set of predogs that share some genes with the refogs.
        For each refog, it can have more than one corresponding predogs.
    """

    V_prime = {}
    for refog_key, refog in U.items():
        V_prime[refog_key] = {}
        for predog_key, predog in V_raw.items():
            if predog_key == "unassigned_genes":
                continue
            overlap = refog & predog
            if len(overlap) != 0:
                if predog_key not in V_prime[refog_key]:
                    V_prime[refog_key][predog_key] = predog
    return V_prime


def V_prime_to_V(V_prime):
    """ "
    V: Dict[predog_key], predog_set]
        Have the same format of V_raw, but it only contains the predog that have some overlaps with the refogs.
    """

    V = {
        predog_key: predog 
        for predog_dict in V_prime.values() for predog_key, predog in predog_dict.items()
    }
    return V


def V_to_V_complete(V, total_missing_genes_set):

    V_complete = copy.deepcopy(V)
    if len(total_missing_genes_set) != 0:
        V_complete["missing_genes"] = total_missing_genes_set

    return V_complete


def prime_information(U, V_complete, n):

    list_products = [*product(U.items(), V_complete.items())]
    AMI = adjusted_mutual_information(U, V_complete, n, list_products)
    AVI = variation_of_information(U, V_complete, n, list_products)
    ARI = rand_index(U, V_complete, n, list_products)
    # NRMI = normalised_reduced_mutual_information(U, V_complete, n, list_products)

    return AMI, AVI, ARI


def TP_FP_FN_union_dict(U, V_prime, n, return_type="count"):

    TP_dict, TP_count_dict = {}, {}
    FN_dict, FN_count_dict = {}, {}
    FP_dict, FP_count_dict = {}, {}
    union_dict, union_count_dict = {}, {}
    for refog_key, refog in U.items():
        TP_dict[refog_key] = {predog_key: refog & predog for predog_key, predog in V_prime[refog_key].items()}
        TP_count_dict[refog_key] = {predog_key: len(TP) for predog_key, TP in TP_dict[refog_key].items()}

        FN_dict[refog_key] = {predog_key: refog - predog for predog_key, predog in V_prime[refog_key].items()}
        FN_count_dict[refog_key] = {predog_key: len(FN) for predog_key, FN in FN_dict[refog_key].items()}

        FP_dict[refog_key] = {predog_key: predog - refog for predog_key, predog in V_prime[refog_key].items()}
        FP_count_dict[refog_key] = {predog_key: len(FP) for predog_key, FP in FP_dict[refog_key].items()}

        union_dict[refog_key] = {predog_key: refog | predog for predog_key, predog in V_prime[refog_key].items()}
        union_count_dict[refog_key] = {
            predog_key: len(union_set) for predog_key, union_set in union_dict[refog_key].items()
        }
    if return_type == "count":
        return TP_count_dict, FN_count_dict, FP_count_dict, union_count_dict
    elif return_type == "both":
        return TP_dict, FN_dict, FP_dict, union_dict, TP_count_dict, FN_count_dict, FP_count_dict, union_count_dict
    elif return_type == "set":
        return TP_dict, FN_dict, FP_dict, union_dict


def micro_scores(U, V_prime, n):

    TP_count_dict, FN_count_dict, FP_count_dict, union_count_dict = TP_FP_FN_union_dict(
        U, V_prime, n, return_type="count"
    )

    total_TP = np.sum([np.sum([*TP_count_dict[refog_key].values()]) for refog_key in TP_count_dict])
    total_FN = np.sum([np.sum([*FN_count_dict[refog_key].values()]) for refog_key in FN_count_dict])
    total_FP = np.sum([np.sum([*FP_count_dict[refog_key].values()]) for refog_key in FP_count_dict])
    total_union = np.sum([np.sum([*union_count_dict[refog_key].values()]) for refog_key in union_count_dict])

    if total_TP != 0 or total_FN != 0:
        micro_recall = total_TP / (total_TP + total_FN)
    else:
        micro_recall = 0.0

    if total_TP != 0 or total_FP != 0:
        micro_precision = total_TP / (total_TP + total_FP)
    else:
        micro_precision = 0.0

    if micro_recall != 0.0 or micro_precision != 0.0:
        micro_f1score = 2 * micro_recall * micro_precision / (micro_recall + micro_precision)
    else:
        micro_f1score = 0.0
    if total_union != 0:
        micro_jaccard_index = total_TP / total_union
    else:
        micro_jaccard_index = 0.0
    micro_fowlkes_mallows_index = np.sqrt(micro_recall * micro_precision)
    micro_dissimilarity = 1 - micro_jaccard_index
    micro_distance = expon.ppf(micro_dissimilarity, scale=len(U))

    return (
        micro_recall,
        micro_precision,
        micro_f1score,
        micro_fowlkes_mallows_index,
        micro_jaccard_index,
        micro_dissimilarity,
        micro_distance,
    )


def fusion(ref_ogs, V_prime, V):
    fusion_predog_dict = {}
    fusion_refog_set = set()

    for predog_key, predog in V.items():
        refog_key_set_multi = set()

        # for refog_key, refog in ref_ogs.items():
        #     if len(refog & predog) != 0:
        #         refog_key_set_multi.add(refog_key)

        # if len(refog_key_set_multi) > 1:
        #     fusion_predog_dict[predog_key] = refog_key_set_multi
        #     fusion_refog_set = fusion_refog_set.union(refog_key_set_multi)
        # elif len(refog_key_set_multi) == 1:
        #     refog_key = refog_key_set_multi.pop()
        #     refog = ref_ogs[refog_key]
        #     if len(predog - refog) > 0 and len(refog - predog) == 0:
        #         refog_key_set_multi.add(refog_key)

        #     fusion_predog_dict[predog_key] = refog_key_set_multi
        #     fusion_refog_set = fusion_refog_set.union(refog_key_set_multi)


        refog_key_set_multi = set([
            refog_key for refog_key, predog_dict in V_prime.items() 
            if predog_key in predog_dict
        ])

        if len(refog_key_set_multi) > 1:
            fusion_predog_dict[predog_key] = refog_key_set_multi
            fusion_refog_set = fusion_refog_set.union(refog_key_set_multi)

    # for refog_key, refog in ref_ogs.items():
    #     for predog_key, predog in V.items():
    #         if len(predog & refog) != 0:
    #             if len(predog - refog) > 0 and len(refog - predog) == 0:
    #                 fusion_refog_set.add(refog_key)
                
    #                 if predog_key not in fusion_predog_dict:
    #                     fusion_predog_dict[predog_key] = set([refog_key])
    #                 else:
    #                     fusion_predog_dict[predog_key].add(refog_key)

    
    # for refog_key, predog_dict in V_prime.items():
    #     if len(predog_dict) == 1:
    #         predog_key = [*predog_dict.keys()][0]
    #         if len(predog_dict[predog_key] - ref_ogs[refog_key]) > 0 \
    #             and len(ref_ogs[refog_key] - predog_dict[predog_key]) == 0:
    #             fusion_refog_set.add(refog_key)
            
    #             if predog_key not in fusion_predog_dict:
    #                 fusion_predog_dict[predog_key] = set([refog_key])
        

    fusion_refog_bool_dict = {
        refog_key: 1 if refog_key in fusion_refog_set else 0
        for refog_key in ref_ogs.keys() 
        
    }
    
    return fusion_refog_set, fusion_predog_dict, fusion_refog_bool_dict


def fission(U, V_prime):
    fission_predog_set = set()
    fission_refog_set = set()
    fission_refog_bool_dict = {}
    for refog_key in U:
        if len(V_prime[refog_key]) > 1:
            fission_refog_set.add(refog_key)
            fission_predog_set = fission_predog_set.union(set([*V_prime[refog_key].keys()]))
            fission_refog_bool_dict[refog_key] = 1
        else:
            fission_refog_bool_dict[refog_key] = 0

    return fission_refog_set, fission_predog_set, fission_refog_bool_dict


def calculate_benchmarks_pairwise(
        ref_ogs: Dict[str, Set[str]],
        pred_ogs: Dict[str, Set[str]],
        # uncert_genes: Optional[Dict] = None,
        q_even: bool = True,
        # q_remove_uncertain: bool = True,
    ):
    referenceOGs = ref_ogs
    predictedOGs = pred_ogs
    totalFP = 0.0
    totalFN = 0.0
    totalTP = 0.0
    totalGroundTruth = 0.0
    n_exact = 0
    n_splits = []
    # so as not to count uncertain genes either way remove them from the
    # expected and remove them from any predicted OG (as though they never existed!)
    totalGenes = 0.0

    for refog_key in referenceOGs:
        refOg = referenceOGs[refog_key]
        # if uncert_genes is not None:
        #     uncert = uncert_genes[refog_key]

        thisFP = 0.0
        thisFN = 0.0
        thisTP = 0.0
        this_split = 0
        # if q_remove_uncertain and uncert_genes is not None:
        #     refOg = refOg.difference(uncert)

        nRefOG = len(refOg)
        totalGenes += nRefOG
        not_present = set(refOg)
        intersection = 0

        for predog_key in predictedOGs:
            predOg = predictedOGs[predog_key]
            overlap = len(refOg.intersection(predOg))

            intersection += overlap

            if overlap > 0:
                # if q_remove_uncertain and uncert_genes is not None:
                #     predOg = predOg.difference(
                #         uncert
                #     )  # I.e. only discount genes that are uncertain w.r.t. this RefOG
                overlap = len(refOg.intersection(predOg))

            if overlap > 0:
                this_split += 1
                not_present = not_present.difference(predOg)

                thisTP += overlap * (overlap - 1) / 2  # n-Ch-2
                thisFP += overlap * (len(predOg) - overlap)
                thisFN += (nRefOG - overlap) * overlap

        # print(intersection_dict[refognum])

        # Are FNs more from splintered OGs or missing genes?
        # print("%f\t%f" % (thisFN/2./(nRefOG-1), len(not_present)*(nRefOG-1)/2./(nRefOG-1)))
        # finally, count all the FN pairs from those not in any predicted OG
        thisFN += len(not_present) * (nRefOG - 1)
        # don't count 'orphan genes' as splits, it's more informative only to count
        # clusters that this orthogroup has been split into. Recall already counts
        #  the missing genes, this should give a distinct measure.
        # this_split += len(not_present)
        # All FN have been counted twice
        assert thisFN % 2 == 0
        n_splits.append(this_split)
        # print(this_split)      # Orthogroup fragments
        # print(len(not_present))  # Unclustered genes
        thisFN /= 2
        # sanity check
        nPairs1 = thisTP + thisFN
        nPairs2 = nRefOG * (nRefOG - 1) / 2
        if nPairs1 != nPairs2:
            print("ERROR: %d != %d" % (nPairs1, nPairs2))

        totalGroundTruth += nPairs1
        if thisFN == 0 and thisFP == 0:
            n_exact += 1
        if q_even:
            N = float(len(refOg) - 1)
            totalFN += thisFN / N
            totalFP += thisFP / N
            totalTP += thisTP / N
        else:
            totalFN += thisFN
            totalFP += thisFP
            totalTP += thisTP

    TP, FP, FN = (totalTP, totalFP, totalFN)
    if TP != 0 or FP != 0:
        pres = TP / (TP + FP)
    else:
        pres = 0.0
    if TP != 0 or FN != 0:
        recall = TP / (TP + FN)
    else:
        recall = 0.0
    if pres != 0 or recall != 0:
        f = 2 * pres * recall / (pres + recall)
    else:
        f = 0.0

    return TP, FP, FN, recall, pres, f

def check_missing_species(refog_species_dict, predog_species_dict, num_round):
    missing_species_dict = {}
    missing_species_count_dict = {}
    missing_species_percent_dict = {}
    
    for refog_key, refog_species in refog_species_dict.items():
        refog_species_num = len(refog_species)
        refog_species_copy = refog_species.copy()
        for predog_key, predog_species in predog_species_dict[refog_key].items():
            refog_species_copy.difference_update(predog_species)

        missing_species_dict[refog_key] = refog_species_copy
        missing_species_count_dict[refog_key] = len(refog_species_copy)
        missing_species_percent_dict[refog_key] = \
            np.round(len(refog_species_copy) / refog_species_num, num_round)
    
    return  missing_species_dict, missing_species_count_dict, missing_species_percent_dict

            
  


if __name__ == "__main__":

    # U = {
    #     1: {"gene1", "gene4", "gene5"}, 
    # }

    # V_raw = {
    #     6: {"gene1", "gene4", "gene5"},
    # }


    # U = {
    #     1: {"gene1", "gene4", "gene5"}, 
    #     2: {"gene2", "gene7", "gene3"}, 
    #     3: {"gene6", "gene8", "gene9", "gene10"}
    # }

    # V_raw = {
    #     6: {"gene1", "gene4"},
    #     8: {"gene8", "gene3"},
    #     7: {"gene2", "gene6", "gene9"},
    #     9: {"gene10", "gene7", "gene5"},
    #     # 12: {"gene22", "gene16"},
    #     # 4: {10: {"gene9", "gene5"}}
    # }


    U = {'61': {'Danio_rerio.ENSDARP00000086805', 'Rattus_norvegicus.ENSRNOP00000007398', 'Mus_musculus.ENSMUSP00000147115', 'Canis_familiaris.ENSCAFP00000046753', 'Homo_sapiens.ENSP00000353590', 'Pan_troglodytes.ENSPTRP00000064140', 'Tetraodon_nigroviridis.ENSTNIP00000006941', 'Danio_rerio.ENSDARP00000124366', 'Caenorhabditis_elegans.WBGene00003777.1', 'Caenorhabditis_elegans.WBGene00003776.1', 'Monodelphis_domestica.ENSMODP00000051065', 'Mus_musculus.ENSMUSP00000090661', 'Homo_sapiens.ENSP00000379616', 'Drosophila_melanogaster.FBpp0072306', 'Rattus_norvegicus.ENSRNOP00000035440', 'Danio_rerio.ENSDARP00000125918', 'Pan_troglodytes.ENSPTRP00000078901', 'Pan_troglodytes.ENSPTRP00000065124', 'Ciona_intestinalis.ENSCINP00000035221', 'Danio_rerio.ENSDARP00000134054', 'Rattus_norvegicus.ENSRNOP00000062744', 'Canis_familiaris.ENSCAFP00000060526', 'Tetraodon_nigroviridis.ENSTNIP00000010900', 'Gallus_gallus.ENSGALP00000043966', 'Homo_sapiens.ENSP00000478109', 'Gallus_gallus.ENSGALP00000053892', 'Monodelphis_domestica.ENSMODP00000010173', 'Pan_troglodytes.ENSPTRP00000069176', 'Canis_familiaris.ENSCAFP00000057173', 'Danio_rerio.ENSDARP00000041141', 'Canis_familiaris.ENSCAFP00000061753', 'Tetraodon_nigroviridis.ENSTNIP00000006335', 'Homo_sapiens.ENSP00000216181', 'Rattus_norvegicus.ENSRNOP00000073929', 'Tetraodon_nigroviridis.ENSTNIP00000017126', 'Danio_rerio.ENSDARP00000114445', 'Monodelphis_domestica.ENSMODP00000001289', 'Mus_musculus.ENSMUSP00000156021', 'Gallus_gallus.ENSGALP00000046782', 'Rattus_norvegicus.ENSRNOP00000072636', 'Tetraodon_nigroviridis.ENSTNIP00000023109', 'Tetraodon_nigroviridis.ENSTNIP00000001018', 'Mus_musculus.ENSMUSP00000016771', 'Homo_sapiens.ENSP00000493594'}, 
         'missing_genes': {'Tetraodon_nigroviridis.ENSTNIP00000006941', 'Tetraodon_nigroviridis.ENSTNIP00000004844', 'Danio_rerio.ENSDARP00000155954', 'Tetraodon_nigroviridis.ENSTNIP00000007040'}}
    V_raw = {'61': {'Ciona_intestinalis.ENSCINP00000035221', 'Drosophila_melanogaster.FBpp0072306', 'Danio_rerio.ENSDARP00000124366', 'Homo_sapiens.ENSP00000353590', 'Tetraodon_nigroviridis.ENSTNIP00000006941', 'Danio_rerio.ENSDARP00000134054', 'Rattus_norvegicus.ENSRNOP00000072636', 'Caenorhabditis_elegans.WBGene00003777.1', 'Tetraodon_nigroviridis.ENSTNIP00000017126', 'Pan_troglodytes.ENSPTRP00000069176', 'Danio_rerio.ENSDARP00000114445', 'Tetraodon_nigroviridis.ENSTNIP00000001018', 'Canis_familiaris.ENSCAFP00000061753', 'Gallus_gallus.ENSGALP00000043966', 'Monodelphis_domestica.ENSMODP00000001289', 'Homo_sapiens.ENSP00000478109', 'Mus_musculus.ENSMUSP00000156021', 'Canis_familiaris.ENSCAFP00000057173', 'Homo_sapiens.ENSP00000216181', 'Rattus_norvegicus.ENSRNOP00000073929', 'Rattus_norvegicus.ENSRNOP00000007398', 'Monodelphis_domestica.ENSMODP00000051065', 'Homo_sapiens.ENSP00000493594', 'Mus_musculus.ENSMUSP00000147115', 'Danio_rerio.ENSDARP00000086805', 'Gallus_gallus.ENSGALP00000053892', 'Monodelphis_domestica.ENSMODP00000010173', 'Pan_troglodytes.ENSPTRP00000065124', 'Rattus_norvegicus.ENSRNOP00000062744', 'Canis_familiaris.ENSCAFP00000046753', 'Mus_musculus.ENSMUSP00000016771', 'Rattus_norvegicus.ENSRNOP00000035440', 'Homo_sapiens.ENSP00000379616', 'Caenorhabditis_elegans.WBGene00003776.1', 'Canis_familiaris.ENSCAFP00000060526', 'Gallus_gallus.ENSGALP00000046782', 'Pan_troglodytes.ENSPTRP00000078901', 'Tetraodon_nigroviridis.ENSTNIP00000023109', 'Tetraodon_nigroviridis.ENSTNIP00000010900', 'Danio_rerio.ENSDARP00000125918', 'Tetraodon_nigroviridis.ENSTNIP00000006335', 'Pan_troglodytes.ENSPTRP00000064140', 'Danio_rerio.ENSDARP00000041141', 'Mus_musculus.ENSMUSP00000090661'}, 
             'missing_genes': {'Tetraodon_nigroviridis.ENSTNIP00000006941', 'Tetraodon_nigroviridis.ENSTNIP00000004844', 'Tetraodon_nigroviridis.ENSTNIP00000007040', 'Danio_rerio.ENSDARP00000155954'}}
    # V_prime = {
    #     1: {6: {"gene1"}, 8: {"gene8", "gene3", "gene4"}},
    #     2: {7: {"gene2", "gene6"}, 9: { "gene10", "gene7"}, 8: {"gene8", "gene3", "gene4"}},
    #     3: {9: { "gene10", "gene7"}, 7: {"gene2", "gene6"}, 8: {"gene8", "gene3", "gene4"}},
    #     # 4: {10: {"gene9", "gene5"}}
    # }

    V_prime = V_raw_to_V_prime(U, V_raw)
    V = V_prime_to_V(V_prime)
    # print(V)

    N = np.sum([len(refog) for refog in U.values()])
    M = np.sum([len(predog) for predog in V.values()])
    M_raw = np.sum([len(predog) for predog in V_raw.values()])

    missing_predogs_list = missing_predogs(V_prime)
    # print("%0.1f%%  Missing prefOGs (%%)" % (100.0 * len(missing_predogs_list) / len(U)))

    (
        total_missing_genes,
        total_missing_genes_set,
        missing_genes_dict,
        missing_genes_count_dict,
        missing_genes_proportion_dict,
    ) = check_missing_genes(U, V_prime, 3)

    # print("%0.1f%%  Missing Genes (%%)" % (100.0 * total_missing_genes / N))

    # fusion_refog_set, fusion_predog_dict = fusion(V_prime, V)
    # fusion_refog_score = 100.0 * len(fusion_refog_set) / len(U)
    # fusion_predog_score = 100.0 * len(fusion_predog_dict) / len(V)
    # print("%0.1f%%  Fussion (refOG) (%%)" % fusion_refog_score)
    # print("%0.1f%%  Fussion (predOG) (%%)" % fusion_predog_score)
    # print()

    # (
    #     micro_recall,
    #     micro_precision,
    #     micro_f1score,
    #     micro_fowlkes_mallows_index,
    #     micro_jaccard_index,
    #     micro_dissimilarity,
    #     micro_distance,
    # ) = micro_scores(U, V_prime, N)

    # print("%0.1f%%  micro Recall (%%)" % (100.0 * micro_recall))
    # print("%0.1f%%  micro Precision (%%)" % (100.0 * micro_precision))
    # print("%0.1f%%  micro F1-score (%%)" % (100.0 * micro_f1score))
    # print("%0.1f%%  micro Fowlkes Mallows Index (%%)" % (100.0 * micro_fowlkes_mallows_index))
    # print("%0.1f%%  micro Jaccard Index (%%)" % (100.0 * micro_jaccard_index))
    # print("%0.1f%%  micro Disimilarity (%%)" % (100.0 * micro_dissimilarity))
    # print("%0.2f   micro Distance" % (micro_distance))
    # print()

    # (
    #     total_weighted_recall,
    #     macro_recall,
    #     weighted_recall_dict,
    #     total_weighted_precision,
    #     macro_precision,
    #     weighted_precision_dict,
    #     total_weighted_f1score,
    #     macro_f1score,
    #     weighted_f1score_dict,
    #     total_weighted_fowlkes_mallows,
    #     macro_fowlkes_mallows,
    #     weighted_fowlkes_mallows_dict,
    #     total_weighted_JI,
    #     macro_JI,
    #     weighted_JI_dict,
    #     total_weighted_dissimilarity,
    #     macro_dissimilarity,
    #     weighted_dissimilarity_dict,
    #     total_weighted_distance,
    #     macro_distance,
    #     weighted_distance_dict,
    #     macro_dissimilarity_distance,
    #     all_score_dict,
    #     all_TP_dict,
    #     effective_size_precision_weighted_dict,
    #     effective_size_JI_weighted_dict,
    #     effective_size_JI_refog_weighted_dict,
    # ) = macro_and_weighted_avg_scores(U, V_prime, N)

    # print("%0.1f%%  macro Recall (%%)" % (100.0 * macro_recall))
    # print("%0.1f%%  macro Precision (%%)" % (100.0 * macro_precision))
    # print("%0.1f%%  macro F1-score (%%)" % (100.0 * macro_f1score))
    # print("%0.1f%%  macro Fowlkes Mallows Index (%%)" % (100.0 * macro_fowlkes_mallows))
    # print("%0.1f%%  macro Jaccard Index (%%)" % (100.0 * macro_JI))
    # print("%0.1f%%  macro Disimilarity (%%)" % (100.0 * macro_dissimilarity))
    # print("%0.2f   macro Distance" % (macro_distance))
    # print("%0.2f   macro Disimilarity Distance" % (macro_dissimilarity_distance))
    # print()

    # print("%0.1f%%  Weighted Avg Recall (%%)" % (100.0 * total_weighted_recall))
    # print("%0.1f%%  Weighted Avg Precision (%%)" % (100.0 * total_weighted_precision))
    # print("%0.1f%%  Weighted Avg F1-score (%%)" % (100.0 * total_weighted_f1score))
    # print("%0.1f%%  Weighted Avg Fowlkes Mallows Index (%%)" % (100.0 * total_weighted_fowlkes_mallows))
    # print("%0.1f%%  Weighted Avg Jaccard Index (%%)" % (100.0 * total_weighted_JI))
    # print("%0.1f%%  Weighted Avg Disimilarity (%%)" % (100.0 * total_weighted_dissimilarity))
    # print("%0.2f   Weighted Avg Distance" % (total_weighted_distance))
    # print()

    # total_entropy, entropy_dict = entropy_score(U, V_prime, N)
    # print("%0.2f  Entropy" % (total_entropy))

    # precision_weighted_KLD = kl_divergence(U, effective_size_precision_weighted_dict, N, M)
    # JI_weighted_KLD = kl_divergence(U, effective_size_JI_weighted_dict, N, M)
    # JI_refog_weighted_KLD = kl_divergence(U, effective_size_JI_refog_weighted_dict, N, M)

    # print("%0.2f  Precision Weighted KLD" % (precision_weighted_KLD))
    # print("%0.2f  JI Weighted KLD" % (JI_weighted_KLD))
    # print("%0.2f  JI refOG Weighted KLD" % (JI_refog_weighted_KLD))
    # print()

    V_complete = V_to_V_complete(V, total_missing_genes_set)
    AMI, AVI, ARI = prime_information(U, V_complete, N)

    # print("%0.2f  Reduced Mutual Information (global)" % (NRMI))
    print("%0.5f  Adjusted Mutual Information (global)" % (AMI))
    print("%0.5f  Adjusted Variation of Information (global)" % (AVI))
    print("%0.5f  Ajusted Rand Index (global)" % (ARI))
    print()

    # print(missing_genes_dict, missing_genes_count_dict, missing_genes_proportion_dict)

    # print(np.sum([*missing_genes_count_dict.values()]))
    # total_misplaced_genes, total_misplaced_genes_proportion, misplaced_genes_dict = check_misplaced_genes(U, V_prime, N)
    # print(misplaced_genes_dict)
    # total_missing_genes, total_missing_genes_proportion, total_missing_genes_dict, total_missing_genes_set = \
    #     check_total_missing_genes(total_misplaced_genes,
    #                               missing_genes_count_dict,
    #                               missing_genes_dict,
    #                               misplaced_genes_dict,
    #                               N)
    # print(total_misplaced_genes)
    # print(total_missing_genes, total_missing_genes_proportion)
    # print(total_missing_genes_dict)
    # print(total_missing_genes_set)

    # prime_information(U, V_prime, total_missing_genes_set, N)
    # U = {
    #     1: {"gene1", "gene4", "gene5"},
    #     2: {"gene2", "gene7"},
    #     3: {"gene6", "gene8", "gene3", "gene9", "gene10"}
    # }

    # V = {
    #     1: {"gene1", "gene4"},
    #     2: {"gene2", "gene7", "gene8"},
    #     3: {"gene6", "gene3", "gene10"},
    #     # 1: {"gene1", "gene4", "gene5"},
    #     # 2: {"gene2", "gene7", "gene8"},
    #     # 3: {"gene6", "gene9", "gene10"},
    #     4: {"gene9", "gene5"}
    # }

    # N = 10
    # list_products = [*product(U.items(), V.items())]
    # NRMI = normalised_reduced_mutual_information(U, V, N, list_products)
    # print(NRMI)

    # AMI = adjusted_mutual_information(U, V, N, list_products)
    # print(AMI)

    U_item_dict = rearange(U)
    V_item_dict = rearange(V)
    label_true_list = []
    label_pred_list = []
    for item, label_true in U_item_dict.items():
        label_true_list.append(label_true)
        label_pred_list.append(V_item_dict[item])

    print(cluster.adjusted_mutual_info_score(label_true_list, label_pred_list))

    print(cluster.rand_score(label_true_list, label_pred_list))
    # print(rand_index(U, V, N, list_products, ri_type="normal"))
    # print(rand_index(U, V, N, list_products))
    print(cluster.adjusted_rand_score(label_true_list, label_pred_list))
