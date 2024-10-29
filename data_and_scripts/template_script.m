%% Template scrip manuscript
% Lorenzo Pasquini, 10/15/2024, lorenzo.pasquini@ucsf.edu

clear all; close all; clc;

% donwload toolbox from Brain Connectivity Toolbox, https://sites.google.com/site/bctnet/
addpath(genpath('~/BCT'));

% download function below from https://www.mathworks.com/matlabcentral/fileexchange/27418-fdr_bh
addpath(genpath('~/fdr_corr/'));

set(0,'defaultfigurecolor',[1 1 1]);

%% fMRI data

cd('~/DMT_func');

inj = 240;
labels = {
    'rvis1','rvis2','rvis3','rvis4','rvis5','rvis6','rvis7','rvis8','rvis9', ...
    'rsm1', 'rsm2', 'rsm3', 'rsm4', 'rsm5','rsm6',...
    'rdan1','rdan2','rdan3','rdan4','rdan5','rdan6','rdan7','rdan8',...
    'rvan1','rvan2','rvan3','rvam4','rvan5','rvan6','rvan7',...
    'rlim1','rlim2','rlim3',...
    'rfp1','rfp2','rfp3','rfp4',...
    'rdmn1','rdmn2','rdmn3','rdmn4','rdmn5','rdmn6','rdmn7','rdmn8','rdmn9','rdmn10','rdmn11','rdmn12','rdmn13',...
    'lvis1','lvis2','lvis3','lvis4','lvis5','lvis6','lvis7','lvis8','lvis9', ...
    'lsm1', 'lsm2', 'lsm3', 'lsm4', 'lsm5','lsm6',...
    'ldan1','ldan2','ldan3','ldan4','ldan5','ldan6','ldan7','ldan8',...
    'lvan1','lvan2','lvan3','lvam4','lvan5','lvan6','lvan7',...
    'llim1','llim2','llim3',...
    'lfp1','lfp2','lfp3','lfp4',...
    'ldmn1','ldmn2','ldmn3','ldmn4','ldmn5','ldmn6','ldmn7','ldmn8','ldmn9','ldmn10','ldmn11','ldmn12','ldmn13',...
    'sc1','sc2','sc3','sc4','sc5','sc6','sc7','sc8','sc9','sc10','sc11','sc12'};

% Load DMT and PCB time series
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

% Plot connectivity dynamics
figure;
for ns=1:size(all_dmt,3)
    subplot(4,4,ns);
    imagesc(corr(all_dmt(:,:,ns)));
    colorbar;
    colormap('jet');
    title(['DMT Sub',num2str(ns)]);
end

figure;
for ns=1:size(all_dmt,3)
    subplot(4,4,ns);
    imagesc(corr(all_pcb(:,:,ns)));
    colorbar;
    colormap('jet');
    title(['placebo Sub',num2str(ns)]);
end

% Generate connectivity dynamics variables
for ns=1:size(all_dmt,3)
    corr_all_dmt(:,:,ns) = corr(all_dmt(:,:,ns));
    corr_all_pcb(:,:,ns) = corr(all_pcb(:,:,ns));
end

%% Modularity

DMNmPLC = mean(corr_all_dmt,3)-mean(corr_all_pcb,3);

rng(1);
thr = .25; %
[Ci, Q]= community_louvain(threshold_absolute(DMNmPLC, thr), 1, [], 'negative_sym');

% recode vectors of relevant states
Ci0 = Ci;
Ci0(Ci~=8) = 0;
Ci0(Ci0~=0) = 1;
Ci1 = Ci;
Ci1(Ci~=33) = 0;
Ci1(Ci1~=0) = 1;
Ci2 = Ci;
Ci2(Ci~=38) = 0;
Ci2(Ci2~=0) = 1;
Ci3 = Ci;
Ci3(Ci~=81) = 0;
Ci3(Ci3~=0) = 1;
Ci4 = Ci;
Ci4(Ci~=27) = 0;
Ci4(Ci4~=0) = 1;

%% One sample t test

vec1 = [Ci0, Ci1, Ci2, Ci3, Ci4]; % regressors

for ns = 1:14 % loop over subjects
    for nroi = 1:112 % loop over ROIs
        % DMT
        tmp = squeeze(all_dmt(nroi,:,ns))';
        tmp2 = regress(tmp, [ones(size(vec1,1),1) vec1]);
        for nb = 1:size(tmp2,1)
            beta(ns,nroi,nb) = tmp2(nb);
        end
        % PCB
        tmp3 = squeeze(all_pcb(nroi,:,ns))';
        tmp4 = regress(tmp3, [ones(size(vec1,1),1) vec1]);
        for nb = 1:size(tmp4,1)
            beta2(ns,nroi,nb) = tmp4(nb);
        end
    end
end

diff_beta = beta-beta2; % DMT - PCB

%% One sample t test maps
% Generate a "~/2_DMT_states_diff/" folder

% Subcortical areas
scaal = niftiread('~/aal_2mm.nii');
scaal = double(scaal);

drep = setdiff(unique(scaal), [71:78 37:38 41:42]);
for nrp = 1:size(drep)
    scaal(scaal==drep(nrp)) = 0;
end

scaal(scaal==71) = 101;
scaal(scaal==72) = 102;
scaal(scaal==73) = 103;
scaal(scaal==74) = 104;
scaal(scaal==75) = 105;
scaal(scaal==76) = 106;
scaal(scaal==77) = 107;
scaal(scaal==78) = 108;
scaal(scaal==37) = 109;
scaal(scaal==38) = 110;
scaal(scaal==41) = 111;
scaal(scaal==42) = 112;

namest = [8, 33, 38, 81, 27];
for nst = 1:5
    % DMT-PCB
    [H,P,CI,STATS] = ttest(diff_beta(:,:,nst+1)); %27, 81, 33
    fdrP = fdr_bh(P);
    tvalue_dmt_roi2(nst,:)= fdrP.*STATS.tstat;

    % Mean t FDR
    mean_t_fdr(:,nst,1) = mean(diff_beta(:,STATS.tstat.*fdrP>0,5),2);
    mean_t_fdr(:,nst,2) = mean(diff_beta(:,STATS.tstat.*fdrP<0,5),2);

    % Cortical areas
    sc100 = niftiread('~/Schaefer2018_100Parcels_7Networks_order_FSLMNI152_2mm.nii');
    sc100 = sc100+scaal; % Cortical and subcortical areas
    sc100(sc100>112) = 0;

    tmap3 = sc100;
    info = niftiinfo('~/Schaefer2018_100Parcels_7Networks_order_FSLMNI152_2mm.nii');

    for nr=1:112
        tmap3(sc100==nr)=fdrP(nr).*STATS.tstat(nr);
    end

    % Save
    myname = sprintf(['~/2_DMT_states_diff/tmap_TWOTAIL_dmt-pcb_state',num2str(namest(nst)),'_FDRp05.nii']);
    niftiwrite(tmap3, myname, info);
end
