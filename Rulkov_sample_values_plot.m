% Sample plots of the Rulkov map, time seris and phase space

%% Parameters for RM
sigma_values =[-.605,-.6025,0.79,0.85];  % Adjusted range for regime transition 
alpha_rulkov = 4.3;               % Chaotic bursting parameter (renamed to avoid conflict) 
mu = 0.001;                % speed of dynamics parameter

% Initial conditions mathing the bifurcation diagram
x0 = -1.5;
y0 =-2; 

num_sgma = length(sigma_values);
% Simulation parameters
N = 20000;        % Total iterations
transient = 5000; % Transient period (fixed value)

% Transient mode test: to see how the transition mode behaves
% N = 3000;     
% transient = 500; 

m=2;
alpha=0.2;
lambda = 0.01;
%% Loop over sigma values and plot time series
    figure
for idx =1:num_sgma

    sigma = sigma_values(idx);
    x=[];
    % Simulate system with state continuation
    [x, y] = rulkov_iteration(alpha_rulkov, mu, sigma, N, x0, y0);
    
    % Extract steady-state portion
    signal = x(transient:end);
    
    % Compute entropy measures on the steady-state signal
    [HD, HA, H_joint,M,range_D] = compute_Sampentropies(signal, lambda,m,alpha);
    sprintf('range_D=%.3f, for sigma=%.3f',range_D, sigma)
    if idx<5;subplot(2,2,idx);elseif idx==5;figure;subplot(2,2,idx-4);elseif idx==6;subplot(2,2,idx-4);end
    plot(signal,'k', 'LineWidth',1.5)
    if HD>HA && HD>H_joint
         title(sprintf('sigma=%.3f, H_A=%.3f, \\color{red}H_D=%.3f\\color{black}, H_{DA}=%.3f',sigma,HA,HD,H_joint))
    elseif HA > HD && HA>H_joint
         title(sprintf('sigma=%.3f,  \\color{red} H_A=%.3f, H_D=%.3f \\color{black} , H_{DA}=%.3f',sigma,HA,HD,H_joint))
    elseif H_joint>HA && H_joint>HD
        title(sprintf('sigma=%.3f, H_A=%.3f, H_D=%.3f,  \\color{red} H_{DA}=%.3f',sigma,HA,HD,H_joint))
    end
    title(sprintf('sigma=%.3f | H_A=%.3f, H_D=%.3f, H_{DA}=%.3f',sigma,HA,HD,H_joint))
    xlim('tight');
   
    ax = gca; % Get current axis
    ax.XColor = 'w'; % Set x-axis color to white (optional!)
    ax.YColor = 'w'; % Set y-axis color to white  (optional!)

    set(gca,'FontSize',20,'FontWeight','bold')   
end
%% Plotting phase space
figure
for idx = 1:num_sgma
    sigma = sigma_values(idx);
    
    % Simulate system with state continuation
    [x, y] = rulkov_iteration(alpha_rulkov, mu, sigma, N, x0, y0);
    
    % Extract steady-state portion
    signal = x(transient:end);
    
    if idx<5;subplot(2,2,idx);elseif idx==5;figure;subplot(2,2,idx-4);elseif idx==6;subplot(2,2,idx-4);end
    plot(signal(1:end-1),signal(2:end),'.','markersize',7)
    axis equal
    set(gca,'FontSize',20,'FontWeight','bold')   
end
