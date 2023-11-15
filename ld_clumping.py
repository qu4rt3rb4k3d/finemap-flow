import copy
import random
import numpy

def clump_disjoint(pvalues, ld_matrix, ld_threshold):
    num_variants = len(pvalues)
    variants = [i for i in range(num_variants)]
    clumps = []
    while variants:
        lead_variant = max(variants, key=lambda v: pvalues[v])
        ld_row = ld_matrix[lead_variant, :].tolist()
        clump_idxs = [i for i in variants if ld_row[i] >= ld_threshold]
        variants = [v for v in variants if v not in clump_idxs]
        clumps.append(clump_idxs)
    return clumps

def clump_disjoint_leads_first(pvalues, ld_matrix, ld_threshold):
    num_variants = len(pvalues)
    leads = [max(clump, key=lambda idx: pvalues[idx]) for clump in clump_disjoint(pvalues, ld_matrix, ld_threshold)]
    num_leads = len(leads)
    clumps = [[] for i in range(num_leads)]
    for i in range(num_variants):
        lds = [ld_matrix[i, j] for j in leads]
        idx = lds.index(max(lds))
        clumps[idx].append(i)
    return clumps

def merge_small_clumps(clumps, min_size, pvalues, ld_matrix):
    leads = [max(clump, key=lambda idx: pvalues[idx]) for clump in clumps]
    merged_clumps = copy.deepcopy(clumps)
    for clump in clumps:
        if len(clump) < min_size:
            for i in clump:
                lds = [ld_matrix[i, j] if j != i and len(lead_clump) >= min_size else -1 for j, lead_clump in zip(leads, clumps)]
                if all([ld <= 0 for ld in lds]):
                    idxs = [i for i, x in enumerate(lds) if x >= 0]
                    new_idx = random.choice(idxs)
                else:
                    new_idx = lds.index(max(lds))
                merged_clumps[new_idx].append(i)
    return [clump for clump in merged_clumps if len(clump) >= min_size]
