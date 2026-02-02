import numpy as np

def extract_DA(time_series, lambda_=0.001):
    """
    Extract duration (D) and cumulative amplitude (A) pairs from a 1D time series
    by segmenting based on sign changes in the raw derivative, with an adaptive
    threshold based on IQR(derivative).

    Parameters
    ----------
    time_series : array_like, shape (N,)
        1D time series.
    lambda_ : float, default=0.001
        Scaling factor for thresholding (threshold = lambda_ * IQR(diff(x))).

    Returns
    -------
    duration_sequence : np.ndarray, shape (M,)
        Duration (number of points) of each segment.
    amp_sequence : np.ndarray, shape (M,)
        Segment amplitude computed as x[i] - x[start_idx] at each cutpoint i.
    segments_idx : np.ndarray, shape (M,)
        Segment end indices i (1-based in MATLAB; here returned 0-based indices).
    """
    x = np.asarray(time_series, dtype=float).ravel()
    if x.size < 2:
        return np.array([], dtype=int), np.array([], dtype=float), np.array([], dtype=int)

    # Step 1: raw derivative (MATLAB diff)
    derivative = np.diff(x)

    # Step 2: adaptive threshold = lambda * IQR(derivative)
    q75, q25 = np.percentile(derivative, [75, 25])
    threshold = float(lambda_ * (q75 - q25))

    # Step 3: segment based on sign changes and jump size condition
    duration_list = []
    amp_list = []
    idx_list = []

    start_idx = 0  # Python 0-based; MATLAB start_idx=1
    # MATLAB: for i = 2:length(derivative)
    for i in range(1, derivative.size):
        if (np.sign(derivative[i]) != np.sign(derivative[i - 1])
            and (abs(derivative[i] - derivative[i - 1]) > threshold)):

            # MATLAB duration = i - start_idx; with MATLAB indices
            # Here, "i" is derivative index; corresponding signal index is i (0-based)
            duration = i - start_idx
            amplitude = x[i] - x[start_idx]

            duration_list.append(duration)
            amp_list.append(amplitude)
            idx_list.append(i)  # signal index where segment ends (0-based)

            start_idx = i

    return (np.asarray(duration_list, dtype=int),
            np.asarray(amp_list, dtype=float),
            np.asarray(idx_list, dtype=int))
