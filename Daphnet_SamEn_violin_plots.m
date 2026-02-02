
%% Sample entropy Analysis with Aggregated Violin Plots
% Computes SampEnt metrics for ankle vertical acceleration
% across subjects in the Daphnet Freezing-of-Gait dataset. For every
% subject/session the metrics are saved, and at the end a single set of
% violin plots compares Walk vs Freeze across three frequency bands.
% -------------------------------------------------------------------------
% Sara Kamali, sara.kamali@uam.es, sara.kamali@gmail.com, UAM, GNB lab, Aug 2025 

clear; clc;
set(groot,'DefaultAxesFontSize',16,'defaultaxesfontweight','bold', ...
    'DefaultAxesLineWidth',1.5,'DefaultLineLineWidth',2);

%% Paths & parameters -
data_path   = 'add path to data';
results_path = 'add path to save results';
addpath(data_path,results_path)
cd(data_path)
if ~exist(results_path,'dir'); mkdir(results_path); end

lambda = 0.01; m = 2; alpha = 0.2;   % ExSEnt parameters
Fs = 64;  % Sampling rate (Hz)
N  = 256; % FIR filter order  → 4 s window

% Design zero-phase FIR band-pass filters
b_delta = fir1(N-1,[0.5  3]/(Fs/2),'bandpass',hamming(N));
b_alpha = fir1(N-1,[3    8]/(Fs/2),'bandpass',hamming(N));

BANDS   = {'0.5–3 Hz','3–8 Hz'};
numBands = numel(BANDS);
subjects = [1:3,5:9];
numSubj = length(subjects);
numSess = 2; %number of sessions

load('-mat','Daphnet_selected_windows.mat')
load('ExEnt_Daphnet_embedding_data.mat')
% ------------------------------------------------------------------------
% Pre-allocate aggregated result arrays  (subjects × bands)
SampEn_walk   = nan(numSubj,numBands);
SampEn_freeze = nan(numSubj,numBands);
%% Main loop --------------------------------------------------------------
cd(data_path)

for sub = 1:numSubj
    subj=subjects(sub);
    tmp = struct('SEw',nan(numSess,numBands), 'SEf',nan(numSess,numBands));

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
            filtfilt(b_alpha,1,x_raw)];

        % start/stop indices of runs imported with the subjects_data struct
        % file
        idxGait   = subject_data(subj).(sprintf('normR%02d',sess));
        idxFreeze = subject_data(subj).(sprintf('freezeR%02d',sess));

        for b = 1:numBands
            SEw=[];SEf=[];
            % Concatenate samples for each condition – robust to empty idx* sets
            x_gait = []; x_freeze = [];
            if ~isempty(idxGait)
                x_gait   = cellfun(@(idx) x_filt(idx(1):idx(2),b), idxGait  , 'UniformOutput',false);
            end
            if ~isempty(idxFreeze)
                x_freeze = cellfun(@(idx) x_filt(idx(1):idx(2),b), idxFreeze, 'UniformOutput',false);
            end
            if isempty(x_gait) || isempty(x_freeze); continue; end

            % SampEn metrics
            for seg = 1:numel(x_gait)
                rtol = alpha * std(x_gait{seg});
                SEw(seg)= sample_entropy(x_gait{seg}  ,m,rtol);
            end

            for seg = 1:numel(x_freeze)
                rtol = alpha * std(x_freeze{seg});
                SEf(seg)= sample_entropy(x_freeze{seg}  ,m,rtol);
            end

            tmp.SEw (sess,b) = nanmean(SEw);  
            tmp.SEf (sess,b) = nanmean(SEf);  
        end
    end

    % Average across sessions
    SampEn_walk  (sub,:) = nanmean(tmp.SEw ,1);
    SampEn_freeze(sub,:) = nanmean(tmp.SEf ,1);

    % Save per-subject structure
    save(fullfile(results_path,sprintf('S%02d_SampEn.mat',subj)), 'tmp');
    tmp =[];
end

save('SampEn_Daphnet_8subject.m','SampEn_walk','SampEn_freeze');

%% Violin plots (SampEn only) --------------------------------------------
aggFig = figure('Position', [100 100 1200 450]);

% Colors (Walk/NFG = light; Freeze/FOG = dark)

colWalk   = [0.2 0.75 0.4];   % light
colFreeze = [0.05 0.4 0.2];   % dark

for b = 1:numBands
    % Gather and scrub per group
    dataW = SampEn_walk(:, b);   dataW = dataW(isfinite(dataW));
    dataF = SampEn_freeze(:, b); dataF = dataF(isfinite(dataF));

    subplot(1, numBands, b); hold on;
    if isempty(dataW) && isempty(dataF)
        title(sprintf('%s — no data', BANDS{b}));
        ylabel('SampEn'); set(gca,'FontSize',22); box off; hold off; continue;
    end

    % Build inputs only for present groups
    vals = []; groups = {};
    if ~isempty(dataW)
        vals   = [vals; dataW(:)];
        groups = [groups; repmat({'NFG'}, numel(dataW), 1)];
    end
    if ~isempty(dataF)
        vals   = [vals; dataF(:)];
        groups = [groups; repmat({'FOG'}, numel(dataF), 1)];
    end
    cats = categorical(groups, {'NFG','FOG'});

% Violin (requires violinplot.m). If absent, fallback to boxplot.
if exist('violinplot','file') == 2
    v = violinplot(vals, cats, 'ViolinAlpha', 0.3, 'ShowData', true);

    hasW = ~isempty(dataW);
    hasF = ~isempty(dataF);

    % Color mapping (note: ViolinColor expects a CELL)
    if hasW && hasF
        % Order matches category order {'NFG','FOG'}
        v(1).ViolinColor = {colWalk};   v(1).BoxColor = colWalk;
        v(2).ViolinColor = {colFreeze}; v(2).BoxColor = colFreeze;

        % Optional: color the scatter points too (guarded)
        if isprop(v(1),'ScatterPlot') && ~isempty(v(1).ScatterPlot)
            try, v(1).ScatterPlot.CData = repmat(colWalk, size(v(1).ScatterPlot.CData,1), 1); end
        end
        if isprop(v(2),'ScatterPlot') && ~isempty(v(2).ScatterPlot)
            try, v(2).ScatterPlot.CData = repmat(colFreeze, size(v(2).ScatterPlot.CData,1), 1); end
        end

    elseif hasW
        v(1).ViolinColor = {colWalk};   v(1).BoxColor = colWalk;
        if isprop(v(1),'ScatterPlot') && ~isempty(v(1).ScatterPlot)
            try, v(1).ScatterPlot.CData = repmat(colWalk, size(v(1).ScatterPlot.CData,1), 1); end
        end

    else % hasF only
        v(1).ViolinColor = {colFreeze}; v(1).BoxColor = colFreeze;
        if isprop(v(1),'ScatterPlot') && ~isempty(v(1).ScatterPlot)
            try, v(1).ScatterPlot.CData = repmat(colFreeze, size(v(1).ScatterPlot.CData,1), 1); end
        end
    end

else
    boxplot(vals, cats);
    set(findobj(gca,'Tag','Box'),'LineWidth',1.7);
    set(findobj(gca,'Tag','Median'),'LineWidth',2);
end
    ylabel('SampEn'); title(BANDS{b});
    set(gca,'FontSize',22); box off; hold off;
end

sgtitle('Aggregated SampEn Across Subjects');
saveas(aggFig, fullfile(results_path, 'AllSubjects_ViolinPlots_SampEn.png'));


%% ------------------- Stats: rank-sum (NFG vs FOG) -----------------------
resBand  = strings(numBands, 1);
nNFG     = NaN(numBands, 1);
nFOG     = NaN(numBands, 1);
pRankSum = NaN(numBands, 1);
hRankSum = NaN(numBands, 1);

for b = 1:numBands
    wVals = SampEn_walk(:, b);   wVals = wVals(isfinite(wVals));
    fVals = SampEn_freeze(:, b); fVals = fVals(isfinite(fVals));

    resBand(b) = BANDS{b};
    nNFG(b) = numel(wVals); 
    nFOG(b) = numel(fVals);

    if ~isempty(wVals) && ~isempty(fVals)
        [pRankSum(b), hRankSum(b)] = ranksum(wVals, fVals, 'tail', 'both');
    end
end

statsTbl = table(resBand, nNFG, nFOG, pRankSum, hRankSum, ...
    'VariableNames', {'Band','nNFG','nFOG','pRankSum','hRankSum'});
disp(statsTbl);


%% -------- Collect every window across subjects (SampEn; bands 1–2) ------
numBandsPlot = min(2, numBands);
valsAll  = cell(1, numBandsPlot);   % values per band
groupAll = cell(1, numBandsPlot);   % 'Sxx_W' / 'Sxx_F' labels

cd(results_path)
for s = 1:numSubj
    subj = subjects(s);
    fname = sprintf('S%02d_SampEn.mat', subj);   % matches your save() above
    if ~isfile(fname), warning('Missing %s – skipped.', fname); continue; end
    L = load(fname, 'tmp');  % fields: tmp.SEw (numSess×numBands), tmp.SEf

    for bIdx = 1:numBandsPlot
        wVals = L.tmp.SEw(:, bIdx);
        fVals = L.tmp.SEf(:, bIdx);

        % scrub NaN/Inf
        wVals = wVals(isfinite(wVals));
        fVals = fVals(isfinite(fVals));

        % accumulate
        valsAll{bIdx}  = [valsAll{bIdx};  wVals;  fVals];
        groupAll{bIdx} = [groupAll{bIdx}; ...
                          repmat({sprintf('S%02d_W', subj)}, numel(wVals), 1); ...
                          repmat({sprintf('S%02d_F', subj)}, numel(fVals), 1)];
    end
end

%% -------------------- Box-plot grid (SampEn; color-coded) ---------------
colWalk   = [0.2 0.75 0.4];   % light
colFreeze = [0.05 0.4 0.2];   % dark
gap = 3;                        % spacing; two boxes per subject

fig = figure('Position',[100 100 1200 450]);

for bIdx = 1:numBandsPlot
    ax = subplot(1, numBandsPlot, bIdx); hold(ax, 'on')
    anyPlotted = false;

    for s = 1:numSubj
        subjID  = subjects(s);
        tagWalk = sprintf('S%02d_W', subjID);
        tagFree = sprintf('S%02d_F', subjID);

        wVals = valsAll{bIdx}( strcmp(groupAll{bIdx}, tagWalk) );
        fVals = valsAll{bIdx}( strcmp(groupAll{bIdx}, tagFree) );

        if isempty(wVals) && isempty(fVals), continue; end

        posW = gap*(s-1)+1;
        posF = posW+1;

        if ~isempty(wVals)
            boxplot(wVals, 'positions', posW, 'Colors', colWalk, ...
                    'Widths', 0.6, 'Symbol', '');
            anyPlotted = true;
        end
        if ~isempty(fVals)
            boxplot(fVals, 'positions', posF, 'Colors', colFreeze, ...
                    'Widths', 0.6, 'Symbol', '');
            anyPlotted = true;
        end
    end

    if anyPlotted
        % Thicker outlines
        set(findobj(ax,'Tag','Median' ),'LineWidth',2)
        set(findobj(ax,'Tag','Whisker'),'LineWidth',1.9)
        set(findobj(ax,'Tag','Box'    ),'LineWidth',1.7)

        % Translucent fills (odd x → Walk, even x → Freeze)
        bH = findobj(ax,'Tag','Box');
        for k = 1:numel(bH)
            x = get(bH(k),'XData'); y = get(bH(k),'YData');
            if mod(round(mean(x)),2)==1, col = colWalk; else, col = colFreeze; end
            patch(x, y, col, 'FaceAlpha', 0.25, 'EdgeColor', 'none', 'Parent', ax);
        end

        xlim([0 gap*numSubj]);
        xticks(gap*(0:numSubj-1)+1.5);
        xticklabels(arrayfun(@(s)sprintf('S%02d',s),subjects,'uni',false));
        ax.XTickLabelRotation = 45;
    else
        title(sprintf('%s — no per-window data', BANDS{bIdx}));
    end

    ylabel('SampEn');
    title(BANDS{bIdx});
    box(ax,'off');  hold(ax,'off')
end

sgtitle('All subjects – per-window SampEn (Bands 1–2)');
saveas(fig, fullfile(results_path,'Boxplots_AllSubjects_Windows_Band1-2_SampEn.png'));
