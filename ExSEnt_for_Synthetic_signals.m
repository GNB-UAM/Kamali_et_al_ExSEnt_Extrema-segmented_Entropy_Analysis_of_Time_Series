% Code to compute ExEnt metrics over 3 synthetic signals. We run the code over 100 times
% runs to get reliable results. Mean and std of all the runs are displayed in a table
% Sara Kamali, sara.kamali@uam.es, sara.kamali@gmail.com, UAM, GNB lab
%--------------------------------------------------------------------------------

numSim = 100;                % number of simulations (fixed seeds)
N = 10000;                    % number of data points
fs = 200;                    % sampling frequency (Hz)
t = 0:1/fs:(N-1)/fs;         % time vector
%ExSEnt related parameters
lambda = 0.01;              % threshold scaling factor
m = 2;                      % embedding dimension
alpha = 0.2;                 % noise tolerance for SampEn

% Define 100 frequencies for the periodic signal (< 25% of fs = 50 Hz)
f_vec = linspace(1, 48, numSim);  

%% Preallocate arrays for aggregated measures and multiscale entropy
% --- Gaussian Noise ---
HD_G_m     = zeros(numSim,1); 
HA_G_m     = zeros(numSim,1); 
Hjoint_G_m = zeros(numSim,1);
SampEn_G_m = zeros(numSim,1); 

% --- Pink Noise ---
HD_P_m     = zeros(numSim,1);  
HA_P_m     = zeros(numSim,1);  
Hjoint_P_m = zeros(numSim,1); 
M_P_m = zeros(numSim,1); 
SampEn_P_m = zeros(numSim,1);  

% --- Brownian Motion ---
HD_B_m     = zeros(numSim,1);
HA_B_m     = zeros(numSim,1);  
Hjoint_B_m = zeros(numSim,1); 
SampEn_B_m = zeros(numSim,1); 

%% Loop over 100 fixed random seeds
for sim = 1:numSim
    rng(sim);  % fixed random seed for reproducibility
    
    % %% 1. Gaussian Noise
    gaussian_noise = randn(N,1);
    r = alpha * std(gaussian_noise);
    [HD_temp, HA_temp, Hjoint_temp, ~, ~, ~, ~, ~] = compute_ExSEnt_metrics(gaussian_noise, lambda, m, alpha);
    HD_G_m(sim) = HD_temp;
    HA_G_m(sim) = HA_temp;
    Hjoint_G_m(sim) = Hjoint_temp;
    SampEn_G_m(sim) = sample_entropy(gaussian_noise, m, r);

    % 2. Pink Noise
    pinkNoiseObj = dsp.ColoredNoise('Color','pink','SamplesPerFrame',N,'NumChannels',1);
    pink_noise = pinkNoiseObj();
    r = alpha * std(pink_noise);
    [HD_temp, HA_temp,  Hjoint_temp, M_temp, ~, ~, ~, ~, ~] = compute_ExSEnt_metrics(pink_noise, lambda, m, alpha);
    HD_P_m(sim) = HD_temp;
    HA_P_m(sim) = HA_temp;
    Hjoint_P_m(sim) = Hjoint_temp;
    M_P_m(sim)=M_temp;
    SampEn_P_m(sim) = sample_entropy(pink_noise, m, r);

    % 3. Brownian Motion
    brownian_motion = cumsum(randn(N,1));
    r = alpha * std(brownian_motion);
    [HD_temp, HA_temp, Hjoint_temp, rangeD_temp, rangeA_temp, ~, ~, ~] = compute_ExSEnt_metrics(brownian_motion, lambda, m, alpha);
    HD_B_m(sim) = HD_temp;
    HA_B_m(sim) = HA_temp;
    Hjoint_B_m(sim) = Hjoint_temp;
    SampEn_B_m(sim) = sample_entropy(brownian_motion, m, r);

end

%% Define signal names and measure names
signals = {'Gaussian','Pink','Brownian'};
measures = {'HD','HA','Hjoint','SampEn'};

% Function to format stats as "median ± std (CV=val)"
formatStats = @(x) sprintf('%.3f ± %.3f (CV=%.3f)', median(x), std(x), std(x)/median(x));

%%Compute statistics for aggregated measures 
Hjoint_P_m = Hjoint_P_m(~isinf(Hjoint_P_m));
results_m = cell(length(signals), length(measures));

results_m{1,1} = formatStats(HD_G_m); results_m{1,2} = formatStats(HA_G_m);
results_m{1,3} = formatStats(Hjoint_G_m); results_m{1,4} = formatStats(SampEn_G_m);
results_m{2,1} = formatStats(HD_P_m); results_m{2,2} = formatStats(HA_P_m);
results_m{2,3} = formatStats(Hjoint_P_m); results_m{2,4} = formatStats(SampEn_P_m);
results_m{3,1} = formatStats(HD_B_m); results_m{3,2} = formatStats(HA_B_m);
results_m{3,3} = formatStats(Hjoint_B_m); results_m{3,4} = formatStats(SampEn_B_m);

ResultsTable_m = table(signals', results_m(:,1), results_m(:,2), results_m(:,3), ...
    'VariableNames', {'Signal', 'HD', 'HA', 'Hjoint', 'SampEn'});
disp('Aggregated Measures for m = 2');
disp(ResultsTable_m);


