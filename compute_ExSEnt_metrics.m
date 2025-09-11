function [HD, HA, H_joint,M,range_D,range_A,r_D,r_A, segment_ids] = compute_Sampentropies(signal,lambda,m,alpha)

% COMPUTE_SAMPENTROPIES
%   Computes Sample Entropy (SampEn) on:
%     (1) durations between successive extrema segments, D
%     (2) cumulative amplitudes across those segments, A
%     (3) the joint (D,A) sequence
%
% INPUTS
%   signal   : 1D time series (row or column vector)
%   lambda   : threshold parameter used inside extract_DA (as % of IQR) to
%              enforce noise tolerance in extrema segmentation
%   m        : embedding dimension for SampEn on the univariate series
%   alpha    : scaling factor for the tolerance r (i.e., r = alpha * std)
%
% OUTPUTS
%   HD         : SampEn of durations D at embedding m (scalar)
%   HA         : SampEn of cumulative amplitudes A at embedding m (scalar)
%   H_joint    : SampEn of the concatenated, normalized joint series (D,A)
%                at embedding 2*m (scalar)
%   M          : number of detected extrema segments
%   range_D    : range(D) = max(D) - min(D)
%   range_A    : range(A) = max(A) - min(A)
%   r_D, r_A   : SampEn tolerances used for D and A respectively
%   segment_ids: indices/labels of segments returned by extract_DA
%
% DEPENDENCIES (must exist on path)
%   extract_DA.m      → [D_vals, A_vals, segment_ids] = extract_DA(signal, lambda)
%   sample_entropy.m  → H = sample_entropy(x, m, r)
%
% NOTES
%   • r choice: standard practice is r = alpha * std(x).
%   • Joint entropy: we embed a single 1D sequence built by interleaving
%     normalized D and A, then use mAD = 2*m to match total state size.
%   • Normalization: To make r_joint interpretable across datasets, we put
%     D and A on comparable scales before interleaving.
%     IMPORTANT: MATLAB’s normalize(x) by default does z-score (mean 0, std 1).
%     If you intend [0,1] range normalization, use normalize(x,'range').

    % Extract durations (D_vals) and cumulative amplitudes (A_vals)
    [D_vals, A_vals, segment_ids] = extract_DA(signal, lambda);
    
    % Check if valid segments were found
    if isempty(D_vals)
         error('No valid segments found in the raw signal with the current threshold.');
    end
    
    % Set parameters for sample entropy computation
    r_D = alpha * std(D_vals); % tolerance for durations (D)
    r_A = alpha * std(A_vals); % tolerance for amplitudes (A)
    
    % Compute sample entropy for durations (D)
    i=1;
    % for m=2:20
    HD(i) =  sample_entropy(D_vals, m, r_D);%sample_entropy(D_vals, m, r_D);
    
    % Compute sample entropy for cumulative amplitudes (A)
    HA(i) = sample_entropy(A_vals, m, r_A);%sampen


    % Compute oint sample entropy for the normalized paired (D, A)
    D_vals_norm=normalize(D_vals);  % Normalizes to [0,1]
    A_vals_norm=normalize(A_vals);  % Normalizes to [0,1]
    joint_data_norm = reshape([D_vals_norm(:) A_vals_norm(:)].',[],1);
    mAD = 2* m;
    
    % r_joint = alpha * std(joint_data(:)); % overall tolerance for joint data
    r_joint = alpha; % Tolerance for joint data since (std=1 for normalized data)
    H_joint(i) = sample_entropy(joint_data_norm, mAD, r_joint);
    % i=i+1;end
    % set(gca,'fontsize',12)
    % figure;subplot(311);plot(2:20,HD,'r-s','LineWidth',2);ylabel('H_D');hold on;subplot(312);plot(2:20,HA,'b-v','LineWidth',2);ylabel('H_A');
    % subplot(313);plot(2:20,H_joint,'k-d','LineWidth',2);ylabel('H_{DA}');xlabel('m');

    M=length(D_vals); %number of segments
    range_D = range(D_vals);%range of duration values
    range_A = range(A_vals);%range of amplitudes
    

end

