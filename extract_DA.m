% extract_DA extracts duration (D) and cumulative amplitude (A) pairs from a
% 1D time series by segmenting based on sign changes in the derivative.

function [duration_sequence, amp_sequence,segments_idx] = extract_DA(time_series, lambda)

% INPUT:
%   time_series - 1D array representing the time series data.
%   lambda      - Scaling factor for thresholding.
%
% OUTPUT:
%   duration_sequence   - Array containing the duration (number of points) of each segment.
%   amp_sequence     - Array containing the cumulative amplitude (sum of values) of each segment.
if nargin<2
    lambda=0.001;
end
    % Step 1: Compute the raw derivative (no smoothing)
    derivative = diff(time_series);

    % Step 2: Compute an adaptive threshold using robust quantiles
    % This threshold is based on the difference between the 75th and 25th percentiles of diff(x)
    threshold = lambda * iqr(derivative);
    
    % Step 3: Segment the signal based on sign changes or low derivative values.
    % No minimum segment length is enforced.
    duration_sequence = [];
    amp_sequence = [];
    segments_idx=[];
    start_idx = 1;
    n = length(time_series);
    
    for i = 2:length(derivative)
        if sign(derivative(i)) ~= sign(derivative(i -1)) && (abs(derivative(i) - derivative(i-1))>threshold)
            % End the current segment
            duration = i - start_idx;
            amplitude = time_series(i) - time_series(start_idx);
            duration_sequence = [duration_sequence, duration];
            amp_sequence = [amp_sequence, amplitude];
            segments_idx = [segments_idx, i ];
            start_idx = i;  % Start a new segment
        end
    end


 
end
