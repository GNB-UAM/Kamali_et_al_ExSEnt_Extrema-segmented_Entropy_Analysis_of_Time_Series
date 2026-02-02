import numpy as np

def sample_entropy(data, m, r):
    """
    Sample Entropy (SampEn) for a 1D series, matching the MATLAB implementation:
      - Chebyshev distance (max absolute difference) between templates
      - counts unordered pairs (i < j)
      - SampEn = -log(A/B), where:
          B = number of similar pairs of length m
          A = number of similar pairs of length m+1

    Parameters
    ----------
    data : array_like, shape (N,)
        Input 1D series.
    m : int
        Embedding dimension.
    r : float
        Tolerance.

    Returns
    -------
    samp_ent : float
        Sample entropy. Returns NaN if N <= m+1. Returns inf if A==0 or B==0.
    """
    x = np.asarray(data, dtype=float).ravel()
    N = x.size
    m = int(m)
    r = float(r)

    if N <= m + 1:
        return np.nan

    # Templates of length m: X shape = (N-m+1, m)
    Xm = np.lib.stride_tricks.sliding_window_view(x, window_shape=m)
    # Templates of length m+1: X1 shape = (N-m, m+1)
    Xm1 = np.lib.stride_tricks.sliding_window_view(x, window_shape=m + 1)

    def count_similar_pairs(X, tol):
        nT = X.shape[0]
        cnt = 0
        for i in range(nT):
            # Chebyshev distance to all j>i
            d = np.max(np.abs(X[i+1:] - X[i]), axis=1)
            cnt += int(np.sum(d <= tol))
        return cnt

    B = count_similar_pairs(Xm, r)
    A = count_similar_pairs(Xm1, r)

    if B == 0 or A == 0:
        return np.inf

    return float(-np.log(A / B))
