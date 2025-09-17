% Plot the bifurcation diagram of Rossler system, color coded with joint ExEnt entropy 
% and evaluate HD, HA, and H_DA for a range of c values
% Sara Kamali June 2025
clc;

%% 1. Rossler Signals
% System's conditions
fs = 10;           % Sampling frequency (Hz)
a = 0.2; 
b = 0.2;
c = 2:0.002:8;  % Parameter range for c
k =1;
init_state = [.1,.1,.1];
rossler_conditions = zeros(length(c), 3);
for i = 1:length(c)
    rossler_conditions(i, :) = [a, b, c(i)];
end

N = 10000;
tspan = 0:1/fs:(N-1)/fs;
% To discard transient portion (first 20% of points)
transient = 500;
% Parameters for computations of InSEnt
lambda = 0.01; % Noise tolerance for segmentation
m=2; % Embedding dimension
alpha=0.2;% Noise tolerance for SampEnt

% Preallocate arrays for the bifurcation diagram and joint entropy values
num_c = length(rossler_conditions);
max_peaks_est = 2000; % Estimate maximum number of peaks per condition
bifurcation_c = zeros(num_c * max_peaks_est, 1);  % Preallocate generously
bifurcation_x = zeros(num_c * max_peaks_est, 1);
bifurcation_Hjoint = zeros(num_c * max_peaks_est, 1);


% Prepare arrays for HD, HA, and H_joint (for the rug plot)
HD_values = zeros(num_c, 1);
HA_values = zeros(num_c, 1);
H_joint_values = zeros(num_c, 1);
M_values = zeros(num_c, 1);

% Preallocate arrays for other entropy measures
samp_en = zeros(num_c, 1);

% Parfor loop for parallel processing
temp_results = cell(1,num_c); %Preallocate cell array
parfor cond = 1:num_c
    params = rossler_conditions(cond, :);
    rossler_eq = @(t, state) k * (...
        [-state(2) - state(3);
         state(1) + params(1)*state(2);
         params(2) + state(3)*(state(1) - params(3))] );
    [X] = ode4(rossler_eq, tspan, init_state);

    signal_ss = X(transient+1:end,1);    % Discard transient portion 

    
    % Extract local maxima using findpeaks on the steady state signal
    [pks, ~] = findpeaks(signal_ss);  % Adjust prominence as needed
    
    % Compute sample entropy measures on the steady state signal
    [HD, HA, H_joint,M,~,~,~,~] = compute_Sampentropies(signal_ss, lambda,m,alpha);


    % Store HD, HA, and H_joint for the rug plot
    HD_values(cond) = HD;
    HA_values(cond) = HA;
    H_joint_values(cond) = H_joint;
    M_values(cond)=M;


    % Compute other entropy metrics on the steady state signal
    r = alpha * std(signal_ss); 
    samp_en(cond) = sample_entropy(signal_ss, m, r);
    
    % Store bifurcation data in a temporary cell array
    temp_data = cell(1,1);
    temp_data{1,1} = [repmat(params(3), length(pks), 1), pks, repmat(H_joint, length(pks), 1)];

    % Store results in a temporary variable that will be merged later.
    temp_results{cond} = temp_data; %Corrected line
    if ismember(c,2:.5:8)
    fprintf('Rossler entropies for c = %.3f done! %d more to go!\n', params(3), cond);
    end

end

% Merge the temporary results from the parfor loop
result_idx = 1;
for cond = 1:num_c
    temp_data = temp_results{cond}{1,1};
    n_pks = size(temp_data,1);
    bifurcation_c(result_idx:result_idx+n_pks-1) = temp_data(:,1);
    bifurcation_x(result_idx:result_idx+n_pks-1) = temp_data(:,2);
    bifurcation_Hjoint(result_idx:result_idx+n_pks-1) = temp_data(:,3);
    result_idx = result_idx + n_pks;
end

% Trim unused array space
bifurcation_c = bifurcation_c(1:result_idx-1);
bifurcation_x = bifurcation_x(1:result_idx-1);
bifurcation_Hjoint = bifurcation_Hjoint(1:result_idx-1);

%% Plotting with Integrated Axes
figure;

% ------- Normalise H_joint to the range [0 1] -------
Hmin = min(bifurcation_Hjoint);
Hmax = max(bifurcation_Hjoint);
Hnorm = (bifurcation_Hjoint - Hmin) ./ (Hmax - Hmin + eps);   % eps avoids ÷0

% Main axes: Bifurcation Diagram (colour-coded by Joint Entropy)
ax1 = axes('Position',[0.1,0.35,0.85,0.6]);
scatter(ax1, bifurcation_c, bifurcation_x, 5, Hnorm, 'filled');  

set(ax1, 'XTickLabel', []);  % Remove x-axis labels
ylabel(ax1, 'Local Maxima of x(t)');
title(ax1, sprintf('Bifurcation Diagram (Color-coded by Joint Entropy), a = %.1f, b = %.1f, N = %d, fs = %d', a, b, N, fs), 'fontweight', 'bold');
% grid(ax1, 'on');
colormap(ax1, jet);
cl = colorbar(ax1);

cl.Label.String = 'Joint Entropy (Normalized)';
xlim(ax1, [min(c) max(c)]);
set(gca, 'fontsize', 24, 'FontWeight', 'bold');

% Secondary axes: Plot  HD and HA
ax2 = axes('Position', [0.1, 0.1, 0.85, 0.2]);
hold(ax2, 'on');
plot(ax2, c, HD_values, 'Color', [0.1, 0.35, 0.65],'LineWidth',2);
plot(ax2, c, HA_values, 'Color', [0.85, 0.2, 0.25],'LineWidth',2);
plot(ax2, c, H_joint_values, 'Color', [0.65, 0.2, 0.75],'LineWidth',2);

xlabel(ax2, 'c');
ylabel(ax2, 'ExEnt');
legend(ax2, {'H_D', 'H_A', 'H_{DA}'}, 'Location', 'best');
% grid(ax2, 'on');
hold(ax2, 'off');
xlim(ax2, [min(c) max(c)]);
xlabel('c')
set(gca, 'fontsize', 30, 'FontWeight', 'bold','LineWidth',1.5);
linkaxes([ax1, ax2], 'x');
savefig('Rossler_ExEnt_bifur_colorcoded.fig')
% print(gcf,'Rossler_ExEnt_bifur_colorcoded.png','-dpng','-r300')
close
%%
figure;
subplot(311)
plot( c, samp_en,  'Color', [0.7, 0.3, 0.1],'linewidth',1.7);hold on
ylabel('SampEn')
ylim([0 0.6])
set(gca, 'fontsize', 30, 'FontWeight', 'bold','LineWidth',1.5);

savefig('Rossler_sampen.fig')
close
save('Rossler_Bifurcation_workspace.mat')