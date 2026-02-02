import numpy as np

def compute_exsent_metrics(signal, lambda_, m, alpha, *, normalize_mode="zscore"):
    """
    Compute ExSEnt metrics (H_D, H_A, H_DA) using Sample Entropy (SampEn) on:
      D: durations between successive extrema segments
      A: cumulative amplitudes across those segments
      joint: interleaved, normalized (D, A) sequence

    Parameters
    ----------
    signal : array_like, shape (N,)
        1D time series.
    lambda_ : float
        Threshold parameter passed to extract_DA (e.g., % of IQR).
    m : int
        Embedding dimension for SampEn on univariate series (D and A).
    alpha : float
        Scaling factor for tolerance r (r = alpha * std).
    normalize_mode : {"zscore", "range"}, default="zscore"
        Normalization applied to D and A before building the joint sequence.
        "zscore" matches MATLAB normalize(x) default behavior (mean 0, std 1).
        "range" performs min-max scaling to [0, 1].

    Returns
    -------
    HD : float
        SampEn of D at embedding m with r_D = alpha * std(D).
    HA : float
        SampEn of A at embedding m with r_A = alpha * std(A).
    H_joint : float
        SampEn of interleaved normalized (D, A) at embedding 2*m with r_joint.
    M : int
        Number of detected extrema segments (len(D)).
    range_D : float
        max(D) - min(D)
    range_A : float
        max(A) - min(A)
    r_D : float
        Tolerance used for D.
    r_A : float
        Tolerance used for A.
    segment_ids : array_like
        Segment indices/labels returned by extract_DA.

    Notes
    -----
    Dependencies expected in scope:
      extract_DA(signal, lambda_) -> (D_vals, A_vals, segment_ids)
      sample_entropy(x, m, r) -> float
    """

    x = np.asarray(signal, dtype=float).ravel()

    # --- dependency: extrema segmentation ---
    D_vals, A_vals, segment_ids = extract_DA(x, lambda_)

    D_vals = np.asarray(D_vals, dtype=float).ravel()
    A_vals = np.asarray(A_vals, dtype=float).ravel()

    if D_vals.size == 0:
        raise ValueError("No valid segments found in the raw signal with the current threshold.")

    # --- tolerances for univariate SampEn ---
    r_D = alpha * np.std(D_vals, ddof=1)
    r_A = alpha * np.std(A_vals, ddof=1)

    # --- univariate entropies ---
    HD = float(sample_entropy(D_vals, m, r_D))
    HA = float(sample_entropy(A_vals, m, r_A))

    # --- normalization for joint embedding ---
    if normalize_mode == "zscore":
        Dn = (D_vals - D_vals.mean()) / np.std(D_vals, ddof=1)
        An = (A_vals - A_vals.mean()) / np.std(A_vals, ddof=1)
        # MATLAB code used r_joint = alpha assuming std≈1 after normalization
        r_joint = float(alpha)
    elif normalize_mode == "range":
        Dn = (D_vals - D_vals.min()) / (D_vals.max() - D_vals.min() + 1e-12)
        An = (A_vals - A_vals.min()) / (A_vals.max() - A_vals.min() + 1e-12)
        # for range scaling, joint std is not fixed; use alpha * std(joint) if desired
        r_joint = float(alpha * np.std(np.r_[Dn, An], ddof=1))
    else:
        raise ValueError("normalize_mode must be 'zscore' or 'range'.")

    # Interleave as [D1, A1, D2, A2, ...]
    joint_data_norm = np.empty(Dn.size + An.size, dtype=float)
    joint_data_norm[0::2] = Dn
    joint_data_norm[1::2] = An

    mAD = 2 * int(m)
    H_joint = float(sample_entropy(joint_data_norm, mAD, r_joint))

    M = int(D_vals.size)
    range_D = float(np.max(D_vals) - np.min(D_vals))
    range_A = float(np.max(A_vals) - np.min(A_vals))

    return HD, HA, H_joint, M, range_D, range_A, r_D, r_A, segment_ids
