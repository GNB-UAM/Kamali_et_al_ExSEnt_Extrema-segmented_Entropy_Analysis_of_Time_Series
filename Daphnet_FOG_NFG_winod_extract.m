%% This code extracts the optimal windows for freeze vs normal movement for PD subjects
% The selection of the windows is done manually by visually inspecting the FOG events.

clear; clc;
set(groot,'DefaultAxesFontSize',16,'defaultaxesfontweight','bold', ...
    'DefaultAxesLineWidth',1.5, 'DefaultLineLineWidth',2);
Fs=64;
%% Paths & parameters -----------------------------------------------------
data_path    = 'add data path';
results_path = 'add path to save results';

subjects = 1:10;
numSubj = length(subjects); 
numSess = 2; %number of sessions
subject_data=[];
%%Import subjects annotation and time data --------------------------------------------------------------
cd(data_path)

for subj = 1:numSubj
    for sess = 1:numSess
        fname = sprintf('S%02dR0%d.txt',subj,sess);
        data  = load(fname);
        ann     = data(:,11);     % Annotation 1=Walk, 2=Freeze
        figure;plot(ann);title(sprintf('subject %d, R%d',subj,sess))
    end
end

subject_data(1).normR01 = {[48000,67190],[88600,106500],[111100,117700],[123900, 136400]};
subject_data(1).freezeR01 = {[70400, 88590],[106545, 111000],[117900, 123816],[136700, 139220]};
subject_data(1).normR02 = {[16060,30300]};
subject_data(1).freezeR02 = {[34260,38250]};

subject_data(2).normR01 = {[42880,54460],[61190,68350]};
subject_data(2).freezeR01= {[54469,61188]};
subject_data(2).normR02 = {[11840,23650],[27350,30770],[38985,63120],[65350,69300],[69572,76790]};
subject_data(2).freezeR02= {[23651,27349],[30773,38982],[63125,65347]};

subject_data(3).normR01 ={[59170,74690],[85685,97870],[113045,123850],[132035,139474]};
subject_data(3).freezeR01={[50168,59163],[74693,85682],[104592,113044],[123853,132033]};
subject_data(3).normR02 ={[16640,20720],[24305,33267]};
subject_data(3).freezeR02={[20725,24303]};

subject_data(5).normR01 ={[35525,44830],[66050,68980]};
subject_data(5).freezeR01={[31800,35520],[52184,62425],[63965,66044],[96843,103818]};
subject_data(5).normR02 ={[21760,27590],[42080,46090],[64000,81200]};
subject_data(5).freezeR02={[36255,39392],[46092,50560],[81209,83906],[85824,92916]};

subject_data(6).normR01 ={[3200,51330],[73035,83180],[94080,119050],[130146,153527]};
subject_data(6).freezeR01={[51333,55192],[58388,62746],[65700,67400],[119000,120950]};
subject_data(6).normR02 ={};
subject_data(6).freezeR02={};


subject_data(7).normR01 ={[46500,53660],[61680,92580],[104900,112380]};
subject_data(7).freezeR01={[40641,46477],[53660,61546],[97690,104830]};
subject_data(7).normR02 ={[16000,29080],[34990,44790]};
subject_data(7).freezeR02={[28950,30000],[31170,34985]};


subject_data(8).normR01 ={[39680,46615],[81920,86480],[90885,95320]};
subject_data(8).freezeR01={[46615,48500],[64644,72229],[113401,122722]};
subject_data(8).normR02 ={};
subject_data(8).freezeR02={};


subject_data(9).normR01 ={[64750,86339],[107250,113000],[127810,138270],[139430,148420]};
subject_data(9).freezeR01={[38430,46510],[52610,64728],[118720,127802]};
subject_data(9).normR02 ={};
subject_data(9).freezeR02={};

cd(results_path)
save('Daphnet_selected_windows.mat','subject_data')



