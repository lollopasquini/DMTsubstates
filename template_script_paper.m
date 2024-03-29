%% Template code to identify brain activity substates
% Lorenzo Pasquini; 03/28/2028; lorenzo.pasquini@ucsf.edu

clear all; close all; clc;
addpath(genpath('~/BCT')); % Brain Connectivity Toolbox, https://sites.google.com/site/bctnet/

%% Load data
% https://github.com/timmer500/DMT_Imaging

unil=1;
for i=[2 3 6 7]
    load(['LongS0',num2str(i),'DMT.mat']);
    all_dmt(:,:,unil) = BOLD_AAL;
    clear BOLD_AAL;
    load(['LongS0',num2str(i),'PCB.mat']);
    all_pcb(:,:,unil) = BOLD_AAL;
    clear BOLD_AAL;
    unil = unil+1;
end

for ii=[10:13 15 17:19 23 25]
    load(['LongS',num2str(ii),'DMT.mat']);
    all_dmt(:,:,unil) = BOLD_AAL;
    clear BOLD_AAL;
    load(['LongS',num2str(ii),'PCB.mat']);
    all_pcb(:,:,unil) = BOLD_AAL;
    clear BOLD_AAL;
    unil = unil+1;
end

%% Generate brain activity similairty matrices
% Individual brain activity similarity matrices
for ns=1:size(all_dmt,3)
    corr_all_dmt(:,:,ns) = corr(all_dmt(:,:,ns));
    corr_all_pcb(:,:,ns) = corr(all_pcb(:,:,ns));
end

DMNmPLC = mean(corr_all_dmt,3)-mean(corr_all_pcb,3); % mean subtraction matrix

%% Brain activity substate identification

rng(1); % for reproducibility
thr = .25; % ideal threshold to apply on subtraction matrix
[Ci, Q]= community_louvain(threshold_absolute(DMNmPLC, thr), 1, [], 'negative_sym');

% Time spent in substates
for sn=1:max(Ci)
    fo(1,sn)=sum(Ci(1:240)==sn)/240; % before injection
    fo(2,sn)=sum(Ci(241:end)==sn)/(length(Ci)-240); % after injection
end

% Chi-square test comparing time spent in substates pre and post injection
for nfo = 1:size(fo,2)
    clear x1 x2;
    x1=[repmat('a',240,1); repmat('b',600,1)];
    x2=[repmat(1,round(fo(1,nfo)*240),1); repmat(2,round(240-fo(1,nfo)*240),1); repmat(1,round(fo(2,nfo)*600),1); repmat(2,round(600-fo(2,nfo)*600),1)];
    [tbl(:,:,nfo), chi2stats(nfo), pi(nfo)] = crosstab(x1,x2);
end
