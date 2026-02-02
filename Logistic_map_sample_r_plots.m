% Create Logistic map time series and phase space plots for 3 sample r values
clear; clc;

% Global styles
set(groot, ...
    'DefaultAxesFontSize', 15, ...
    'DefaultAxesFontWeight', 'bold', ...
    'DefaultAxesLineWidth', 1.5, ...
    'DefaultLineLineWidth', 2);   % Axis colours in white

% System parameters for the logistic map
r_values =[3.66, 3.742, 3.97];   % Logistic-map parameters (chaotic regime)
x0 = 0.5;      % Initial condition
N = 50000;  % Length of time series

% Sample-entropy parameters
m = 2;
lambda  = 0.01;
alpha = 0.2;

% Pre-allocate array
x = zeros(N, 1);
% Create two subplots; one for time series, one for phase space
figure
for i = 1:numel(r_values)
    r      = r_values(i);
    x(1)   = x0;
    % Generate logistic-map time series
    for n = 1:N-1
        x(n+1) = r * x(n) * (1 - x(n));
    end
    % ---------- subplot 1 : Time-series subplot ----------
    subplot(3, 3,[2+3*(i-1),3+3*(i-1)]);

    idx = 1000:1400;   % Window to visualise
    if i==1
    plot(idx, x(idx), 'k', 'LineWidth', 2,'Marker','x','MarkerEdgeColor','m','MarkerSize',6);
    else
        plot(idx, x(idx), 'k', 'LineWidth', 2);
    end
    xlim([idx(1) idx(end)]);

    % Compute sample-entropy measures (discard transient 1:500)
    % [HD, HA, HDA, ~, ~, ~, ~, ~, ~] = compute_Sampentropies(x(500:end), lambda, m, alpha);
    % 
    % title(sprintf('r = %.3f | H_A = %.3f, H_D = %.3f, H_{DA} = %.3f', ...
    %       r, HA, HD, HDA));
    % set(gca, 'XColor', 'w', 'YColor', 'w');

    % subplot 2 : Phase-space subplot ----------

    subplot(3, 3, 1+3*(i-1));

    plot(x(idx-1), x(idx), '.','color',[0.7,0.2,0.7], 'MarkerSize', 10);
    xlabel('x(i-1)'); ylabel('x(i)');
    title(sprintf('Phase space | r = %.3f', r));
    axis square
    
    set(gca, 'XColor', 'w', 'YColor', 'w');
end

