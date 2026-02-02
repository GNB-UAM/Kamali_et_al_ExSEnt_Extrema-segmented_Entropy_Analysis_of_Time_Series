% Bifurcation of the Rulkov map over sigma, color-coded by joint ExEnt (H_DA)
% + HD, HA, H_DA rug plot
% + Sample Entropy vs sigma
% clear; clc;

%% Control parameter (sigma) and fixed parameters
sigma_range = -1.5:0.0025:1.5;      % control parameter (adjust step as needed)
num_s = numel(sigma_range);
alphaR = 4.3;                       % fast subsystem gain
mu     = 0.001;                  %  timescale ---> slow dynamics

%% Simulation parameters
N           = 20000;              % total iterations
N_transient = 5000;               % discard first samples
N_ss        = 1000;                  % number of steady-state points kept per sigma
x0 = -1.5;
y0 =-2;            % initial conditions

%% ExEnt / SampEn parameters (match your logistic script)
lambda = 0.01;      % segmentation tolerance
m      = 2;         % embedding dimension
alpha  = 0.2;       % SampEn tolerance factor (r = alpha * std)

%% Preallocate rug arrays
HD_values       = zeros(num_s,1);
HA_values       = zeros(num_s,1);
H_joint_values  = zeros(num_s,1);
% Optional: store M if you want it in rug (kept for completeness)
M_values        = zeros(num_s,1);

%% Preallocate bifurcation cloud (exact size)
bifurcation_sigma  = zeros(num_s * N_ss, 1);
bifurcation_x      = zeros(num_s * N_ss, 1);
bifurcation_Hjoint = zeros(num_s * N_ss, 1);
result_idx = 1;

%% Sample entropy per sigma
samp_en = zeros(num_s,1);

%% Sweep sigma
for i = 1:num_s
    sigma = sigma_range(i);

    % --- iterate Rulkov map ---
    [x, y] = rulkov_iteration(alphaR, mu, sigma, N, x0, y0);

    % steady-state segment
    x_ss = x(N_transient:end);

    % --- ExEnt metrics (HD, HA, H_DA) ---
    [HD, HA, H_joint, M, ~, ~, ~, ~, ~] = compute_ExSEnt_metrics(x_ss, lambda, m, alpha);

    % take last N_ss points for bifurcation cloud (like logistic script)
    seg = x_ss(end-N_ss+1:end);

    
    % --- Sample Entropy ---
    r  = alpha * std(x_ss);
    SE = sample_entropy(x_ss, m, r);

    % store rug values
    HD_values(i)      = HD;
    HA_values(i)      = HA;
    H_joint_values(i) = H_joint;
    M_values(i)       = M;
    samp_en(i)        = SE;

    % store bifurcation cloud (replicate H_joint for this stripe)
    idx = result_idx : result_idx + N_ss - 1;
    bifurcation_sigma(idx)  = sigma;
    bifurcation_x(idx)      = seg(:);
    bifurcation_Hjoint(idx) = H_joint;
    result_idx = result_idx + N_ss;

    % progress
    fprintf('sigma = %+5.3f done! %d more to go!\n', sigma, num_s - i);
end

%% Normalize H_joint to [0,1] for color
Hmin = min(H_joint_values); 
Hmax = max(H_joint_values);
Hscaled = (bifurcation_Hjoint - Hmin) ./ (Hmax - Hmin + eps);

%% Plot: integrated axes (bifurcation + rug)
figure('Color','w');

% Main axes: bifurcation diagram (color = H_DA)
ax1 = axes('Position',[0.1, 0.35, 0.85, 0.6]);
scatter(ax1, bifurcation_sigma, bifurcation_x, 4, Hscaled, 'filled');
set(ax1, 'XTickLabel', []);
ylabel(ax1, 'x(n)');
title(ax1, sprintf('Rulkov Bifurcation (Color-coded by H_{DA}), \\alpha=%.2f, \\mu=%.3g, N=%d', alphaR, mu, N));
colormap(ax1, jet); colorbar(ax1);
set(gca,'FontSize',22,'FontWeight','bold','LineWidth',1.5);
xlim(ax1, [min(sigma_range) max(sigma_range)]);

% Secondary axes: rug plot (HD, HA, H_DA)
ax2 = axes('Position',[0.1, 0.1, 0.85, 0.2]);
hold(ax2, 'on');
plot(ax2, sigma_range, HD_values,      'Color', [.1, .35, .65], 'LineWidth', 1.7);
plot(ax2, sigma_range, HA_values,      'Color', [.85, .2, .25], 'LineWidth', 1.7);
plot(ax2, sigma_range, H_joint_values, 'Color', [.65, .2, .75], 'LineWidth', 1.7);
xlabel(ax2, 'Control Parameter \sigma');
ylabel(ax2, 'ExEnt');
legend(ax2, {'H_D','H_A','H_{DA}'}, 'Location', 'best');
set(gca,'FontSize',22,'FontWeight','bold','LineWidth',1.5);
xlim(ax2, [min(sigma_range) max(sigma_range)]);
linkaxes([ax1, ax2], 'x');
% print(gcf,'Rulkovmap_bifur_colorcoded.png','-dpng','-r300')

%% Sample Entropy figure 
figure;
subplot(312)
plot(sigma_range, samp_en, 'Color', [0.7, 0.3, 0.1], 'LineWidth', 1.7);
xlabel('\sigma'); ylabel('SampEn');
set(gca,'FontSize',24,'FontWeight','bold','LineWidth',1.5);
xlim([min(sigma_range) max(sigma_range)]);
% print(gcf,'Rulkovmap_SampEn.png','-dpng','-r300')

%% Save workspace
save('Rulkov_ExEnt_bifurcation_workspcae_2.mat')

