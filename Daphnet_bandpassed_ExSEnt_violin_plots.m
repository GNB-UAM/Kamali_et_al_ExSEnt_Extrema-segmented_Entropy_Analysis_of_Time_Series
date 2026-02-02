%% ExEnt Metrics Analysis with Aggregated Violin Plots
% Computes HD, HA and HDA SampEnt metrics for ankle vertical acceleration
% across subjects in the Daphnet Freezing-of-Gait dataset. For every
% subject/session the metrics are saved, and at the end a single set of
% violin plots compares Walk vs Freeze across three frequency bands.
%Sara Kamali, sara.kamali@uam.es, sara.kamali@gmail.com, UAM, GNB lab

% -------------------------------------------------------------------------

clear; clc;
set(groot,'DefaultAxesFontSize',20,'defaultaxesfontweight','bold', ...
    'DefaultAxesLineWidth',1.5,'DefaultLineLineWidth',2);

%% Paths & parameters -----------------------------------------------------
data_path   = 'add data path';
results_path = 'add path to save results';
addpath(data_path,results_path)
if ~exist(results_path,'dir'); mkdir(results_path); end

lambda = 0.01; m = 2; alpha = 0.2;   % SampEn parameters
Fs = 64;                               % Sampling rate (Hz)
N  = 256;                              % FIR filter order  → 4 s window

% Design zero-phase FIR band-pass filters
b_delta = fir1(N-1,[0.5  3]/(Fs/2),'bandpass',hamming(N));
b_alpha = fir1(N-1,[3    8]/(Fs/2),'bandpass',hamming(N));
b_beta  = fir1(N-1,[8   20]/(Fs/2),'bandpass',hamming(N));
BANDS   = {'0.5–3 Hz','3–8 Hz','8–20 Hz'};
numBands = numel(BANDS);
subjects = [1:3,5:9];
numSubj = length(subjects);
numSess = 2;  %Number of sessions

load('-mat','Daphnet_selected_windows.m')
% ------------------------------------------------------------------------
% Pre-allocate aggregated result arrays  (subjects × bands)
HD_walk   = nan(numSubj,numBands);
HA_walk   = nan(numSubj,numBands);
HDA_walk  = nan(numSubj,numBands);
num_seg_walk = nan(numSubj,numBands);
HD_freeze = nan(numSubj,numBands);
HA_freeze = nan(numSubj,numBands);
HDA_freeze= nan(numSubj,numBands);
num_seg_freeze = nan(numSubj,numBands);

%% Main loop --------------------------------------------------------------

for sub = 1:numSubj
    subj=subjects(sub);
    tmp = struct('HDw',nan(numSess,numBands), 'HAw',nan(numSess,numBands), ...
        'HDAw',nan(numSess,numBands), 'Nsw',nan(numSess,numBands), ...
        'HDf',nan(numSess,numBands), 'HAf',nan(numSess,numBands), ...
        'HDAf',nan(numSess,numBands), 'Nsf',nan(numSess,numBands));

    for sess = 1:numSess
        fname = sprintf('S%02dR0%d.txt',subj,sess);
        if ~isfile(fname)
            warning('Missing %s  →  subject %d session %d skipped.',fname,subj,sess);
            continue;
        end

        data_path  = load(fname);
        x_raw = data_path(:,3);      % Ankle vertical (mg)
        y     = data_path(:,11);     % Annotations
        time = data_path(:,1);      %time vector

        % Zero-phase band-pass filtering (columns = bands)
        x_filt = [filtfilt(b_delta,1,x_raw), ...
            filtfilt(b_alpha,1,x_raw), ...
            filtfilt(b_beta ,1,x_raw)];

        % start/stop indices of runs imported with the subjects_data struct
        % file
        idxGait   = subject_data(subj).(sprintf('normR%02d',sess));
        idxFreeze = subject_data(subj).(sprintf('freezeR%02d',sess));

        for b = 1:numBands
            HDw=[];HAw=[];HDAw=[];Nsw=[];
            HDf=[];HAf=[];HDAf=[];Nsf=[];
            % Concatenate samples for each condition – robust to empty idx* sets
            x_gait = []; x_freeze = [];
            if ~isempty(idxGait)
                x_gait   = cellfun(@(idx) x_filt(idx(1):idx(2),b), idxGait , 'UniformOutput', false);
            end
            if ~isempty(idxFreeze)
                x_freeze = cellfun(@(idx) x_filt(idx(1):idx(2),b), idxFreeze, 'UniformOutput', false);
            end
            if isempty(x_gait) || isempty(x_freeze); continue; end

            % SampEn metrics
            for seg = 1:numel(x_gait)
                [HDw(seg),HAw(seg),HDAw(seg),~,~,~,~,~,num_sw] = compute_ExSEnt_metrics(x_gait{seg}, lambda, m, alpha);
                Nsw(seg)=numel(num_sw);
            end
            for seg = 1:numel(x_freeze)
                [HDf(seg),HAf(seg),HDAf(seg),~,~,~,~,~,num_sf] = compute_ExSEnt_metrics(x_freeze{seg},lambda, m, alpha);
                Nsf(seg)=numel(num_sf);
            end

            tmp.HDw (sess,b) = nanmean(HDw);  tmp.HAw (sess,b) = nanmean(HAw);  tmp.HDAw(sess,b) = nanmean(HDAw);
            tmp.Nsw(sess,b)=nanmean(Nsw);
            tmp.HDf (sess,b) = nanmean(HDf);  tmp.HAf (sess,b) = nanmean(HAf);  tmp.HDAf(sess,b) = nanmean(HDAf);
            tmp.Nsf(sess,b)=nanmean(Nsf);
        end
    end

    % Average across sessions
    HD_walk(sub,:) = nanmean(tmp.HDw ,1);
    HA_walk(sub,:) = nanmean(tmp.HAw ,1);
    HDA_walk(sub,:) = nanmean(tmp.HDAw,1);
    num_seg_walk(sub,:) = nanmean(tmp.Nsw,1);
    HD_freeze(sub,:) = nanmean(tmp.HDf ,1);
    HA_freeze(sub,:) = nanmean(tmp.HAf ,1);
    HDA_freeze(sub,:) = nanmean(tmp.HDAf,1);
    num_seg_freeze(sub,:) = nanmean(tmp.Nsf,1);

    % Save per-subject structure
    save(fullfile(results_path,sprintf('S%02d_metrics.mat',subj)), 'tmp');
    tmp =[];
end

save('ExEnt_Daphnet_8subject.mat','HD_walk','HA_walk','HDA_walk','num_seg_walk','HD_freeze','HA_freeze','HDA_freeze','num_seg_freeze');

%% Violin plots -----------------------------------------------------------
aggFig = figure('Position', [100 100 1200 900]);
metrics = {'HD', 'HA', 'HDA'};     % rows
conds   = {'Move', 'Freeze'};      % plotting labels
% palette ----------------------------------------------------
lightPurp =[0.43 0.16 0.55];  % Walk / Move
darkPurp  =   [0.7 0.6 0.9];     % Freeze

for b = 1:numBands                              % columns
    % Gather metric columns for current band
    dataW = [HD_walk(:, b),  HA_walk(:, b),  HDA_walk(:, b) ];
    dataF = [HD_freeze(:, b), HA_freeze(:, b), HDA_freeze(:, b)];

    % Scrub Inf for plotting
    dataW(isinf(dataW)) = NaN;
    dataF(isinf(dataF)) = NaN;

    for mIdx = 1:numel(metrics)                  % rows
        subplot(numel(metrics), numBands, (mIdx-1)*numBands + b);
        vals   = [dataW(:, mIdx); dataF(:, mIdx)];
        % Option A: set order at creation
        lbls   = [repmat({'NFG'}, numSubj, 1); repmat({'FOG'}, numSubj, 1)];
        groups = categorical(lbls, {'NFG','FOG'});   % NFG on the left, FOG on the right
        % violinplot(vals, groups);
% draw violins + points (default is ShowData=true, so just omit or set true)
v = violinplot(vals, groups, 'ViolinAlpha', 0.3, 'ShowData', true);

% Color the two violins  
cmap = {darkPurp, lightPurp};

% Move / Walk
v(1).ViolinColor = cmap(1);
v(1).BoxColor    = cmap{1};
v(1).ScatterPlot.CData = repmat(cmap{1}, numel(v(1).ScatterPlot.CData), 1);  % dots

% Freeze
v(2).ViolinColor = cmap(2);
v(2).BoxColor    = cmap{2};
v(2).ScatterPlot.CData = repmat(cmap{2}, numel(v(2).ScatterPlot.CData), 1);

        % Row labels
        if b == 1
            ylabel(metrics{mIdx});
        end
        % Column titles
        if mIdx == 1
            title(BANDS{b});
        end
        set(gca,'FontSize',22)
    end
end

sgtitle('Aggregated SampEnt Metrics Across Subjects');
saveas(aggFig, fullfile(results_path, 'AllSubjects_ViolinPlots.png'));
% close(aggFig);

%% ------------------- 4. Stats: rank‑sum & Spearman ----------------------
numMetrics = numel(metrics);
resBands   = strings(numBands * numMetrics, 1);
resMetric  = strings(numBands * numMetrics, 1);
pRankSum   = NaN(numBands * numMetrics, 1);
hRankSum   = NaN(numBands * numMetrics, 1);
rhoSpearm  = NaN(numBands * numMetrics, 1);
pSpearm    = NaN(numBands * numMetrics, 1);

row = 1;
for b = 1:numBands
    walkMat   = [HD_walk(:, b),  HA_walk(:, b),  HDA_walk(:, b)];
    freezeMat = [HD_freeze(:, b), HA_freeze(:, b), HDA_freeze(:, b)];

    % Ensure Inf→NaN for stats consistency
    walkMat(isinf(walkMat))   = NaN;
    freezeMat(isinf(freezeMat)) = NaN;

    for mIdx = 1:numMetrics
        wVals = walkMat(:, mIdx);
        fVals = freezeMat(:, mIdx);

        % Rank‑sum (Wilcoxon) — ignores NaN automatically
        [p_rs, h_rs] = ranksum(wVals, fVals, 'tail', 'both');

        % Store
        resBands(row)  = BANDS{b};
        resMetric(row) = metrics{mIdx};
        pRankSum(row)  = p_rs;
        hRankSum(row)  = h_rs;
        row = row + 1;
    end
end

statsTbl = table(resBands, resMetric, pRankSum, hRankSum,  ...
    'VariableNames', {'Band', 'Metric', 'pRankSum', 'hRankSum'});

disp(statsTbl);

%% -------- Collect every window across subjects (for bands 1–2) ----------
metrics       = {'HD','HA','HDA'};          % three SampEn metrics
numBandsPlot  = 2;                          % only the first two bands
cd(results_path)
valsAll  = cell(numel(metrics), numBandsPlot);   % values
groupAll = cell(numel(metrics), numBandsPlot);   % group labels

for sub = 1:numSubj
    subj = subjects(sub);
    fname = sprintf('S%02d_metrics.mat', subj);
    if ~isfile(fname), warning('Missing %s – skipped.', fname); continue; end
    load(fname, 'tmp');                     % brings HDw/HDf… inside tmp

    for mIdx = 1:numel(metrics)
        metricFields = struct( ...
            'HD', struct('w','HDw','f','HDf'), ...
            'HA', struct('w','HAw','f','HAf'), ...
            'HDA',struct('w','HDAw','f','HDAf'));

        wField = metricFields.(metrics{mIdx}).w;
        fField = metricFields.(metrics{mIdx}).f;

        for bIdx = 1:numBandsPlot
            wVals = tmp.(wField)(:, bIdx);
            fVals = tmp.(fField)(:, bIdx);

            % scrub NaN / Inf
            wVals = wVals(~isnan(wVals) & ~isinf(wVals));
            fVals = fVals(~isnan(fVals) & ~isinf(fVals));

            % accumulate
            valsAll{mIdx,bIdx}  = [valsAll{mIdx,bIdx};  wVals;  fVals];
            gWalk   = repmat({sprintf('S%02d_W', subj)}, numel(wVals), 1);
            gFreeze = repmat({sprintf('S%02d_F', subj)}, numel(fVals), 1);
            groupAll{mIdx,bIdx} = [groupAll{mIdx,bIdx}; gWalk; gFreeze];
        end
    end
end


%% -------------------- Box-plot grid (colour-coded) ----------------------
% Walk  = light purple, Freeze = dark purple
lightPurp = [0.78 0.73 0.95];
darkPurp  = [0.43 0.16 0.55];
gap = 3;                                   % distance step; 2 boxes per subject

fig = figure('Position',[100 100 1200 800]);

for mIdx = 1:numel(metrics)                % rows
    for bIdx = 1:numBandsPlot              % cols (first two bands)
        ax = subplot(numel(metrics),numBandsPlot,(mIdx-1)*numBandsPlot + bIdx);
        hold on

        for s = 1:numSubj
            subjID  = subjects(s);
            tagWalk = sprintf('S%02d_W',subjID);
            tagFree = sprintf('S%02d_F',subjID);

            wVals = valsAll{mIdx,bIdx}( strcmp(groupAll{mIdx,bIdx},tagWalk) );
            fVals = valsAll{mIdx,bIdx}( strcmp(groupAll{mIdx,bIdx},tagFree) );

            posW = gap*(s-1)+1;
            posF = posW+1;

            boxplot(wVals, 'positions',posW, 'Colors',lightPurp, ...
                    'Widths',0.6, 'Symbol','');   % Walk (light)
            boxplot(fVals, 'positions',posF, 'Colors',darkPurp, ...
                    'Widths',0.6, 'Symbol','');   % Freeze (dark)
        end

% 1) Thicker outlines
set(findobj(gca,'Tag','Median' ),'LineWidth',2)  % median line
set(findobj(gca,'Tag','Whisker'),'LineWidth',1.9)  % whiskers
set(findobj(gca,'Tag','Box'    ),'LineWidth',1.7)  % box edges

% 2) Transparent box faces (since boxplot has no native fill)
bH = findobj(gca,'Tag','Box');                     % each box is a patch
for k = 1:numel(bH)
    x = get(bH(k),'XData');   y = get(bH(k),'YData');
    % choose colour based on box centre (odd → walk, even → freeze)
    if mean(x) == floor(mean(x))          % integer pos
        col = lightPurp;
    else
        col = darkPurp;
    end
    patch(x,y,col,'FaceAlpha',0.25,'EdgeColor','none');  % translucent fill
end


        % Cosmetic tweaks
        xlim([0 gap*numSubj]);
        xticks(gap*(0:numSubj-1)+1.5);
        if mIdx==3
            xticklabels(arrayfun(@(s)sprintf('S%02d',s),subjects,'uni',false));
            ax.XTickLabelRotation = 45;
        else
            xticklabels({''})
        end


        if bIdx==1, ylabel(metrics{mIdx}); end
        if mIdx==1, title(sprintf('%s',BANDS{bIdx})); end
        box off
        hold off

    end
end

sgtitle('All subjects – per-window SampEn (Bands 1–2)');
saveas(fig, fullfile(results_path,'Boxplots_AllSubjects_Windows_Band1-2_Color.png'));
% close(fig);
