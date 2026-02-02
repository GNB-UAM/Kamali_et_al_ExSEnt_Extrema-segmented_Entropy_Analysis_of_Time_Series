% Plot Rossler time series and 2D phase space for sample C values
% and evaluate ExSEnt metrics (HD, HA, and H_DA) 
% Sara Kamali, sara.kamali@uam.es, sara.kamali@gmail.com, UAM, GNB lab, June 2025

clear; clc;

%% System & plotting parameters
fs = 10; % sampling frequency 
a = 0.2;
b = 0.2;
cSet = [4,  6.51];  % Sample c-values (chaotic regime) You can change it to try different C values!
k = 1;  % overall scaling
init_state = [0.1 0.1 0.1];    % [x0 y0 z0]

N       = 10000;
tspan   = 0:1/fs:(N-1)/fs;
transient = 500;    % discard transient mode

% ExEnt parameters
lambda = 0.01;
m= 2;
alpha  = 0.2;

% Global plot style
set(groot, ...
    'DefaultAxesFontSize', 15, ...
    'DefaultAxesFontWeight', 'bold', ...
    'DefaultAxesLineWidth', 1.5, ...
    'DefaultLineLineWidth', 2);


h=figure;h.WindowState='Maximized';clf;  
% Loop over the c-values
for kCond = 1:numel(cSet)
    c = cSet(kCond);

    %  the Rossler System
    rossler_eq = @(t, s) k * [ ...
        -s(2) - s(3); ...
         s(1) + a * s(2); ...
         b + s(3) * (s(1) - c) ];
    X = ode4(rossler_eq, tspan, init_state);   % [x y z]

    x_ss = X(transient:end, 1);  % steady-state x-series
    y_ss = X(transient:end, 2);

    % ---------- entropy metrics ----------
    [HD, HA, HDA, ~, ~, ~, ~, ~, ~] = ...
        compute_ExSEnt_metrics(x_ss', lambda, m, alpha);

    % Time series ----------
    subplot(2, 2, kCond);
    idxWin = 100:2000;
    plot(idxWin, x_ss(idxWin), 'k');
    title(sprintf('c = %.2f | H_D = %.3f, H_A = %.3f, H_{DA} = %.3f', ...
          c, HD, HA, HDA));
    xlim([idxWin(1) idxWin(end)]);
    set(gca, 'XColor', 'w', 'YColor', 'w');    % hide axes

    % Phase space ----------
    subplot(2, 2, 2+kCond);
    plot(x_ss, y_ss, '.','color',[.7,.2,.7], 'MarkerSize', 6);
    xlabel('x'); ylabel('y');
    title(sprintf('Phase space | c = %.2f', c));
    axis square;                               % <- make axes square
    set(gca, 'FontSize', 15, 'FontWeight', 'bold');
end
