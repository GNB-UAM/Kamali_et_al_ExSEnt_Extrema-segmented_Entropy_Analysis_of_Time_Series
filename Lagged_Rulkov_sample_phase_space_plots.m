%% Rulkov map
clear; clc;

% Global styles
set(groot, ...
    'DefaultAxesFontSize', 15, ...
    'DefaultAxesFontWeight', 'bold', ...
    'DefaultAxesLineWidth', 1.5, ...
    'DefaultLineLineWidth', 1.7);

%  Rulkov parameters
sigma_values =[-.605,-.6025,.79,.85];  % 4 regimes to visualize behavior
alphaR = 4.3;               % Chaotic bursting parameter (renamed to avoid conflict with alpha for SampEn)
mu = 0.001;  % Slow dynamics parameter

% Initial condition the same as the bifurcation diagram
x0     = -1.5;
y0 = -2;
N      = 20000;              % length of time series
Ntrns =5000;    % Transition mode to get to the stable range

% Sample-entropy parameters
m=2;
lambda = 0.01;
alpha  = 0.2;
%% Create two figure windows – one for time series, one for phase space

h=figure(1);h.WindowState='maximized';

for i = 1:numel(sigma_values)-2 % This is just to plot the 2D and 3D phase space of the first two sigma values
    sigma = sigma_values(i);

    %  Generate Rulkov-map time series
    [xinit, y] = rulkov_iteration(alphaR, mu, sigma, N, x0, y0);
    x1=xinit(Ntrns:end);

    % Laged Phase-space 3D subplot
    subplot(2, 2, i);
    plot3(x1(1:end-2),x1(2:end-1),x1(3:end), '.','color',[0.7,0,0.7]); % This shows the

    xlabel('x[n-1]'),  ylabel('x[n]'),  zlabel('x[n+1]');title(sprintf('sigma=%.3f',sigma))
    axis square
    subplot(2, 2, i+2);

    % Only plot 100 data points to see the trajectory between the 3D states
    xv=1:100;
    xl=x1(xv);
    plot3(xl(1:end-2),xl(2:end-1),xl(3:end), '-o','color',[0.7,0,0.75],'MarkerSize',3,'LineWidth',0.5);title(sprintf('sigma=%.3f',sigma))

    xlabel('x[n-1]'),  ylabel('x[n]'),  zlabel('x[n+1]')
    axis square
end

%% Second set of plot

h=figure(2);h.WindowState='maximized';
Nlong=35000; %Longer simulation range

sigma = sigma_values(3);
[xinit, y] = rulkov_iteration(alphaR, mu, sigma, Nlong, x0, y0);

% 2D phase space with lag 1
len_plot=[15000, 2500];

lag=[1,3];
for i=1:2
    subplot(2, 2, i);
    x=xinit(500: 499+len_plot(i));
    plot(x(1:end-lag(i)), x(1+lag(i):end), '.', 'MarkerSize', 5,'color',[0.65,.2,.75],'LineWidth',.8);
    xlabel('x[n]'); ylabel(sprintf('x[n+%d]',lag));
    title(sprintf('sigma = %.3f,N=%d', sigma,length(x)));
    axis square

end

sigma = sigma_values(4);
[xinit, y] = rulkov_iteration(alphaR, mu, sigma, Nlong, x0, y0);

% 2D phase space with lag 1
len_plot=[2500,15000];

lag=3;
for i=1:2
    subplot(2, 2, 2+i);
    if i==1;x=xinit(500: 499+len_plot(i));else;x=xinit(Ntrns:Ntrns+len_plot(i)-1);end
    plot(x(1:end-lag), x(1+lag:end), '.', 'MarkerSize', 5,'color',[0.65,.2,.75],'LineWidth',.8);
    xlabel('x[n]'); ylabel(sprintf('x[n+%d]',lag));
    title(sprintf('sigma = %.3f,N=%d', sigma,length(x)));
    axis square

end
