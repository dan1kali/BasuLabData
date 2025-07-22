%% Preprocessing for SVM Decoding 
clear
load('state_map.mat','demoTable') %Subject List
badTrialSubs=[]; % Flag for mis-matched trials

% R = Region IDs in Parcellation_Sided_v3, RegionNames= ROIs corresponding to each ID
R=["L1","R1","L2","R2","L4","R4","L6","R6","L9","R9","L11","R11","L13","R13","L14","R14"];
RegionNames=[{'L dlPFC'},{'R dlPFC'},{'L dmPFC'},{'R dmPFC'},{'L OFC'},{'R OFC'},{'L vlPFC'},{'R vlPFC'},{'L LTL'},{'R LTL'},{'L dACC'},{'R dACC'},{'L AMY'},{'R AMY'},{'L HC'},{'R HC'}];

parfor (ss=1:length(demoTable.Subject),4)
    filename=demoTable.ID{ss};
    
    %Load Data
    origdata=load(filename,'ft_data3','Parcellation_Sided_v3','TrialDet','ch_ictal','Trials_C','Trials_I');
    ft_data3=origdata.ft_data3;
    Parcellation_Sided_v3=origdata.Parcellation_Sided_v3;
    ch_ictal=origdata.ch_ictal;
    TrialDet=origdata.TrialDet;
    
    %Remove trials with NaNs
    cfg=[];
    cfg.trials = ~cell2mat(cellfun(@anynan, ft_data3.trial, 'UniformOutput', false)); %finds trials that have nans, selects all others to be included
    if any(cfg.trials==0)
        disp('Empty trials Detected');
        ft_data3 = ft_selectdata(cfg,ft_data3);
    end

    %Check number of trials matches with TrialDet
    if length(ft_data3.trial)~=length(TrialDet)
        badTrialSubs(ss) = 1;
        continue;
    end

    %Fix time axis if not consistent
    for ii=1:length(ft_data3.time)
        ft_data3.time{ii}=ft_data3.time{1};
    end

    %Bandstop line noise + 0.5Hz highpass
    cfg = [];
    cfg.hpfilter = 'yes';
    cfg.hpfreq = 0.5;
    cfg.hpfiltord = 5;
    cfg.bsfilter = 'yes';
    cfg.bsfreq = [55 65];
    cfg.bsfiltord = 4;
    %cfg.plotfiltresp='yes' %For testing/visualizing filters
    %cfg.trials=1;
    ft_data3_rs = ft_preprocessing(cfg,ft_data3);
    
    %Downsample data to match UCMC sampling rate (you don't need to do this)
    if ft_data3_rs.fsample ~=512
        cfg=[];
        cfg.resamplefs = 512;
        cfg.detrend = 'no';
        ft_data3_rs = ft_resampledata(cfg,ft_data3_rs);
    end
    
    goodChans = find(~ismember(ft_data3.label,ch_ictal)); %Select non-ictal channels
    
    %Calculate PSD
    cfg = [];
    cfg.output = 'pow';
    cfg.keeptrials = 'yes';
    cfg.method = 'wavelet';
    cfg.pad='nextpow2'; %pads data to match desired freq resolution (see: https://www.fieldtriptoolbox.org/faq/spectral/freqanalysis_paddinginsufficient/)
    cfg.foi = [1:110]; %1-110Hz w/1Hz steps
    cfg.toi = [-1:(1/512*12):2]; % -1s to +2s w.r.t. image onset, w/23.4ms steps. If you leave sampling rate at 1000Hz, would just do [-1:0.01:2] instead. 
    ft_freq_rs = ft_freqanalysis(cfg,ft_data3_rs);
    
    %Save data
    origdata.ft_freq_rs = ft_freq_rs;
    origdata.goodChans=goodChans;
    origdata=rmfield(origdata,'ft_data3');
    origdata.Trials_C = origdata.Trials_C;
    origdata.Trials_I = origdata.Trials_I;
    save(['D:\Data\SVM\Data\',demoTable.ID{ss},'.mat'],'-fromstruct',origdata);
end