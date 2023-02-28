import numpy as np


def consensus(x):
    """
    Select dominant variants from each alignment.
    At ties, select randomly
    """

    # zeroing non-dominant variants, keeping ties of dominant variants
    d = x.copy()
    non_dominant = d < d.max(axis=2)[:, :, np.newaxis]
    d[non_dominant] = 0

    # from ties of dominant variants randomly select variant
    p_pos = np.array([0, 1, 2, 3])
    ties = d[(d > 0).sum(axis=2) > 1, :]
    for i, p in enumerate(ties):
        rand_select = np.array([0, 0, 0, 0])
        tie_idx = p_pos[(p > 0)]
        rand_select_idx = np.random.choice(tie_idx, 1)[0]
        rand_select[rand_select_idx] = p[rand_select_idx]
        ties[i] = rand_select
    d[(d > 0).sum(axis=2) > 1, :] = ties

    return d
