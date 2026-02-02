% Compute ExEnt over logistic map. Color-code the bifurcation diagram based on H_joint and add
% the HD and HA to the plot

clear; clc;

% Define parameter range for chaos control parameter r
r_range = 2.9:0.002:4;   %  Finer step for smoother bifurcation plot

num_r = length(r_range);

% % Preallocate arrays to plot (HD and HA vs. r)
HD_values = zeros(num_r, 1);
HA_values = zeros(num_r, 1);
H_joint_values = zeros(num_r, 1);
% % Preallocate arrays to plot (Sample entropy)
samp_en = zeros(num_r, 1);


% Prepare arrays for the bifurcation diagram and joint entropy values
bifurcation_r = [];
bifurcation_x = [];
bifurcation_Hjoint = [];

% Simulation parameters
N = 10000;         % number of iterations.
N_transient = 500;  % Number of transient iterations to discard.
N_ss = 500;      % Number of steady-state points to use for bifurcation.
x0 = 0.5;         % initial condition

for i = 1:num_r
    r = r_range(i);

    % Generate logistic map time series: x(n+1) = r*x(n)*(1-x(n))
    x = zeros(N, 1);
    x(1) = x0;
    for n = 1:(N-1)
        x(n+1) = r * x(n) * (1 - x(n));
    end

    % Discard transient portion (first N_transient iterations)
    x_ss = x(N_transient:end);

    % Collect the LAST N_ss points for the bifurcation diagram.
    bifurcation_r = [bifurcation_r; repmat(r, N_ss, 1)];
    bifurcation_x = [bifurcation_x; x_ss(end-N_ss+1:end)];


    % Define noise threshold
    lambda = .01;
    m=2;
    alpha=0.2;
    % Compute DEnt measures (HD, HA, and joint H)
    [HD, HA, H_joint,~,~,~,~,~, ~] = compute_ExSEnt_metrics(x_ss,lambda,m,alpha);


    % Replicate H_joint_norm for each point in the current steady-state
    bifurcation_Hjoint = [bifurcation_Hjoint; repmat(H_joint, N_ss, 1)];

    % Save HD and HA for the rug plot
    HD_values(i) = HD;
    HA_values(i) = HA;
    H_joint_values(i) = H_joint;

    % Compute other metrics for comparision
     rtol = alpha * std(x_ss); % tolerance for durations
    samp_en(i)=  sample_entropy(x_ss, m, rtol);
    fprintf('r=%.3f done! %d more to go!\n',r,length(r_range)-i)
end

%% Plot: Create integrated axes
figure;

% Create main axes for the bifurcation diagram
ax1 = axes('Position',[0.1, 0.35, 0.85, 0.6]);
% scatter(ax1, bifurcation_r, bifurcation_x, 4, bifurcation_Hjoint, 'filled');
H_min = min(H_joint_values);
H_max = max(H_joint_values);
bifurcation_Hscaled = (bifurcation_Hjoint - H_min) / (H_max - H_min + eps);
scatter(ax1, bifurcation_r, bifurcation_x, 4, bifurcation_Hscaled, 'filled');

xticklabels('')
ylabel(ax1, 'x(n)');
title(ax1, 'Bifurcation Diagram (Color-coded by Joint Entropy) for the Logistic Map');
% grid(ax1, 'on');
colormap(ax1, jet);
cl = colorbar(ax1);
cl.Label.String = 'Joint Entropy';  % Add label to colorbar
ylim(ax1, [0 1]);
set(gca,'fontsize',26)
% Create secondary axes for the rug plot (embedded below the main plot)
ax2 = axes('Position',[0.1, 0.1, 0.85, 0.2]);
hold(ax2, 'on');
plot(ax2, r_range, HD_values, 'Color', [.1, .35, .65],'linewidth',1.7);hold on
plot(ax2, r_range, HA_values,  'Color', [.85,.2,.25],'linewidth',1.7);hold on
plot(ax2, r_range, H_joint_values, 'Color', [.65, .2, .75],'linewidth',1.7);
xlabel(ax2, 'Control Parameter r');
ylabel(ax2, 'ExEnt');
legend(ax2, {'H_D', 'H_A','H_{DA}'}, 'Location', 'best');
% grid(ax2, 'on');
hold(ax2, 'off');
 xlim(ax2, [2.9 4-.005]);
set(gca,'fontsize',30,'linewidth',1.7)
% Link the x-axes so they share the same range
linkaxes([ax1, ax2], 'x');

print(gcf,'Logisticmap_bifur_colorcoded.png','-dpng','-r300') %Save the plot

%% Plot
figure
hold 'on';
subplot(312)
plot(r_range,samp_en, 'Color', [0.7, 0.3, 0.1],'linewidth',1.7);hold on
xlabel('Control Parameter r');
ylabel('SampEnt');
xlim([2.9 4-.005]);
set(gca,'fontsize',30,'fontweight','bold','linewidth',1.7)
ylim([0 .9])
print(gcf,'Logisticmap_sampen.png','-dpng','-r300')


%% save workspace
save('Logistic_workspace.mat')

