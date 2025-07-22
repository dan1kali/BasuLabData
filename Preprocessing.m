%% Preprocessing

%% Filtering
% Apply zero-phase digital notch filter (Matlab function "filtfilt") 
% to broadband signals to remove the AC line frequency at 60 Hz and harmonics

cfg = [];

cfg.bsfilter = 'yes';                    % use bandstop filter (notch)
cfg.bsfreq = [55 65; 115 125; 175 185];  % notch bands around 60, 120, 180 Hz
cfg.bsfiltord = 4;

cfg.hpfilter = 'yes';
cfg.hpfreq = 0.5;
cfg.hpfiltord = 5;

ft_data3_filt = ft_preprocessing(cfg,ft_data3);


%% Artifact Rejection
% For each electrode and each task, trials with amplitudes 
% (max-min voltage from fixation onset to stimulus off) larger than 
% three stdevs above the mean amplitude across all trials were excluded





