
tic
files = {'BW42.mat', 'MG51b.mat'};
features = struct();

for i = 1:length(files)
    patientFeatures = preProcess(files{i});
    
    fNames = fieldnames(patientFeatures);
    for k = 1:numel(fNames)
        features.(fNames{k}) = patientFeatures.(fNames{k});
    end
end

save('features.mat','features');
clear
toc
% okay these shits took 15 min

%% functions

function [features] = preProcess(filename)
    load(fullfile(filename));
    [~, patient, ~] = fileparts(filename);

    %%%%%%%%%%%%%%%%% Filtering %%%%%%%%%%%%%%%%%

    cfg = [];
    cfg.bsfilter = 'yes';                    % use bandstop filter (notch)
    cfg.bsfreq = [55 65; 115 125; 175 185];  % notch bands around 60, 120, 180 Hz
    cfg.bsfiltord = 4;
    
    cfg.hpfilter = 'yes';
    cfg.hpfreq = 0.5;
    cfg.hpfiltord = 5;
    
    ft_data3_filt = ft_preprocessing(cfg,ft_data3);
       
    nTrials = numel(ft_data3.trial);
    nChannels = numel(ft_data3.label);
    

    %%%%%%%%%%%%%%%%% Artifact Trial Rejection %%%%%%%%%%%%%%%%%

    % Calculate max-min amplitudes
    amplitudes = zeros(nChannels, nTrials);
    for tr = 1:nTrials
        for ch = 1:nChannels
            signal = ft_data3.trial{tr} (ch, :); % {trial} [channels x time]
            amplitudes(ch, tr) = max(signal) - min(signal);
        end
    end
    
    % Compute thresholds: mean + 3*std per channel
    mean_amp = mean(amplitudes, 2);               % [1 x channels]
    std_amp = std(amplitudes, 0, 2);              % [1 x channels]
    thresholds = mean_amp + 3 * std_amp;          % [1 x channels]
    
    % Identify artifact trials per channel
    artifact_mask = amplitudes > thresholds;      % [trials x channels]
    bad_trials = any(artifact_mask, 1);           % [1 x trials]
    bad_trials = bad_trials';                     % convert to column [trials x 1]
    clean_trials_idx = find(~bad_trials);         % indices of good trials
    
    ft_data_clean = ft_data3;
    ft_data_clean.trial = ft_data3_filt.trial(clean_trials_idx);
    ft_data_clean.time = ft_data3.time(clean_trials_idx);
    ft_data_clean.sampleinfo = ft_data3.sampleinfo(clean_trials_idx, :);
    % New congruent/incongruent indices
    [~, loc_C] = ismember(Trials_C, clean_trials_idx);
    [~, loc_I] = ismember(Trials_I, clean_trials_idx);
    Trials_C_clean = loc_C(loc_C > 0); % remove rejected indices
    Trials_I_clean = loc_I(loc_I > 0);

    fprintf('\nRejected %d of %d trials (%.2f%%)\n', ...
        sum(bad_trials), nTrials, 100 * sum(bad_trials) / nTrials);

    %%%%%%%%%%%%%%%%% Feature Extraction %%%%%%%%%%%%%%%%%
    
    timeData = ft_data_clean.time;
    congruentData = ft_data_clean.trial(Trials_C_clean);
    incongruentData = ft_data_clean.trial(Trials_I_clean);
    selectedChannels = ft_data3.label;
    mask1 = ~ismember(selectedChannels, ch_ictal);
    mask2 = ~strcmp(Parcellation_Sided_v3, 'RNan') & ~strcmp(Parcellation_Sided_v3, 'LNan');
    selectedChannels = find(mask1 & mask2);

    [conPowerFeatures] = extractPowerFeatures(congruentData, timeData);
    [inPowerFeatures] = extractPowerFeatures(incongruentData, timeData);
    
    features = struct();
    % features.(['conMeanPower_' patient]) = conMeanPower;
    % features.(['conMaxPower_' patient]) = conMaxPower;
    features.(['conPowerFeatures_' patient]) = conPowerFeatures;
    % features.(['inMeanPower_' patient]) = inMeanPower;
    % features.(['inMaxPower_' patient]) = inMaxPower;
    features.(['inPowerFeatures_' patient]) = inPowerFeatures;
    features.(['selectedChannels_' patient]) = selectedChannels;

end


% concat_index = ((c-1)*nTrials)+ tr; % concatenate all electrodes


function [band_power_mean_max] = extractPowerFeatures(data, timeData)
    
    nTrials = length(data);
    nChannels = size(data{1}, 1);
            
    band_power_mean_max = cell(1, nTrials); % {trial} [column 1 = mean, column 2 = max]

    for tr=1:nTrials
        band_power_mean_max{tr} = zeros(nChannels, 2);
        for c=1:nChannels

            [~,t,~,~,Snorm]=computeNormalizedFreqMag(data{tr}(c,:),1000);
            t = t + timeData{1}(1);  % shift t so its zero corresponds to 0 s in original time
                        
            S_epoched = Snorm(t >= 0, :);   % time_indices after image onset
       
        % [1st column - mean, 2nd column - max]
        band_power_mean_max{tr}(c,1) = mean(S_epoched(:));  % mean over time
        band_power_mean_max{tr}(c,2) = max(S_epoched(:));   % max over time

        end
        % band_power = [mean power, max power] for all trials per channel
    end
        
    % avg_band_power = mean(band_power_mean_max, 1);  % [average mean power, average max power]
    % sem_band_power = std(band_power_mean_max, 0, 1) ./ sqrt(size(band_power_mean_max, 1));  % SEM for each

end

function [y,t,S,f,Snorm]=computeNormalizedFreqMag(x,Fs,params)

if(nargin < 4)
    params.Fs = Fs;
end
if(~isfield(params,'Fs'))
    params.Fs = Fs;
end
if(~isfield(params,'fpass'))
    params.fpass = [70 120];
end
if(~isfield(params,'tapers'))
    params.tapers = [5 7];
end
if(~isfield(params,'movingwin'))
    params.movingwin = [0.2 0.01]; % [window size, step size]
end
if ~isfield(params,'norm_option')
    params.norm_option='mean';
end
if ~isfield(params,'err')
    params.err=[1 0.05];
end

[S,t,f,~] = mtspecgramc(x,params.movingwin,params); % Chronux function
% units of spectrum? 

if strcmp(params.norm_option,'median')
    m=median(S,1);
elseif strcmp(params.norm_option,'mean')
    m=mean(S,1);
end


Snorm = S./repmat(m, size(S,1), 1);  % normalized power (this is implicit in your code)

y=mean(Snorm,2); % mean power for each trial and channel

% get max and min of S

end


%% unused old functions

function [newdata,nTrials,nChannels] = switchChannelsTrials(data)

    nTrials = length(data);
    nChannels = size(data{1}, 1);
    newdata = cell(1, nChannels);
    
    for c = 1:nChannels
        for t = 1:nTrials
            newdata{c}(t, :) =  data{t}(c, :);
        end
    end

    nChannels = length(newdata);
    nTrials = size(newdata{1}, 1);
end



%% Future Concerns

% 1) how ft_data3_filt is passed to ft_data3_clean only for .trials
% otherwise ft_data3_clean takes everything from ft_data3
% relies on time and sampleinfo matching

% but if doesn't take from ft_data3, then no initialization for
% Trials_C/Trials_I


% 2) how preprocess loads data
% instead of taking data as argument for everything, loads and names
% ft_data3 --> not dynamic


% 3) ROI and ictal don't match the labels?