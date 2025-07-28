
%% just testing conflict mod on BW42


[conflictModChan] = conflictModAnalysis(PowerTimeData, features.PowerData_BW42, features.Zscores_BW42, Trials_C, Trials_I, responseTimes);


%% 
tic
files = {'BW42.mat', 'MG51b.mat'};
features = struct();

for i = 1:length(files)
    patientFeatures = preProcessTrials(files{i});
    
    fNames = fieldnames(patientFeatures);
    for k = 1:numel(fNames)
        features.(fNames{k}) = patientFeatures.(fNames{k});
    end
end

save('features.mat','features');
toc

% okay these took 15 min
% huzzah! only 10 min now!

%% functions

function [clean_trials_idx] = preProcessTrials(filename)
    load(fullfile(filename));
    [~, patient, ~] = fileparts(filename);

    %%%%%%%%%%%%%%%%% Filtering %%%%%%%%%%%%%%%%%

    cfg = [];
    cfg.bsfilter = 'yes';                    % bandstop filter (notch)
    cfg.bsfreq = [55 65; 115 125; 175 185];  % 60, 120, 180 Hz
    cfg.bsfiltord = 4;
    
    cfg.hpfilter = 'yes';                    % high pass filter
    cfg.hpfreq = 0.5;
    cfg.hpfiltord = 5;
    
    ft_data3_filt = ft_preprocessing(cfg,ft_data3);
       
    nTrials = numel(ft_data3.trial);
    nChannels = numel(ft_data3.label);
    

    %%%%%%%%%%%%%%%%% Artifact Trial Rejection %%%%%%%%%%%%%%%%%

    % [-1 +2] window needed for artifact rejection

    % Calculate max-min amplitudes
    amplitudes = zeros(nChannels, nTrials);
    for tr = 1:nTrials
        for ch = 1:nChannels
            signal = ft_data3.trial{tr} (ch, :); % {trial} [channels x time]
            amplitudes(ch, tr) = max(signal) - min(signal);
        end
    end
    
    % Compute thresholds: mean + 3*std per channel
    mean_amp = mean(amplitudes, 2);           % [1 x channels]
    std_amp = std(amplitudes, 0, 2);          % [1 x channels]
    thresholds = mean_amp + 5 * std_amp;      % [1 x channels] - 5 stdevs
    
    % Identify artifact trials per channel
    artifact_mask = amplitudes > thresholds;      % [trials x channels]
    bad_trials = any(artifact_mask, 1) | isnan(TrialDet(:,12))';
    bad_trials = bad_trials';                     % convert to column [trials x 1]
    clean_trials_idx = find(~bad_trials);         % indices of good trials
    
    %%%%%%%%%%%%%%%%% 27 mustnot be 0 or -1 , aka must be 1 %%%%%%%%%%%%%%%%%
    


end


function [features] = PowerAndChannelReject(cleanTrialsIdx)

    % convert rest of ft_data3
    ft_data_clean = ft_data3;
    ft_data_clean.trial = ft_data3_filt.trial(cleanTrialsIdx);
    ft_data_clean.time = ft_data3.time(cleanTrialsIdx);
    ft_data_clean.sampleinfo = ft_data3.sampleinfo(cleanTrialsIdx,:);
    responseTimes = TrialDet(cleanTrialsIdx,12);
    meanResponseTime = mean(responseTimes);
    
    % New congruent/incongruent indices
    [~, loc_C] = ismember(Trials_C, cleanTrialsIdx);
    [~, loc_I] = ismember(Trials_I, cleanTrialsIdx);
    Trials_C_clean = loc_C(loc_C > 0); % remove rejected indices
    Trials_I_clean = loc_I(loc_I > 0);

    fprintf('\nRejected %d of %d trials (%.2f%%)\n', ...
        sum(bad_trials), nTrials, 100 * sum(bad_trials) / nTrials);

    %%%%%%%%%%%%%%%%% Feature Extraction %%%%%%%%%%%%%%%%%
    
    timeData = ft_data_clean.time;


    selectedChannels = ft_data3.label;
    mask1 = ~ismember(selectedChannels, ch_ictal);
    mask2 = ~strcmp(Parcellation_Sided_v3, 'RNan') & ~strcmp(Parcellation_Sided_v3, 'LNan');
    selectedChannels = find(mask1 & mask2);

    [~, ~, PowerData, PowerTimeData] = extractPowerFeatures(ft_data_clean.trial, timeData, responseTimes);

    conPowerFeatures = PowerFeatures(Trials_C_clean);
    inPowerFeatures  = PowerFeatures(Trials_I_clean);

    features = struct();
    features.(['powerData_' patient]) = PowerData;
    features.(['powerTimeData_' patient]) = PowerTimeData;
    features.(['conPowerFeatures_' patient]) = conPowerFeatures;
    features.(['inPowerFeatures_' patient]) = inPowerFeatures;
    features.(['selectedChannels_' patient]) = selectedChannels;
    features.(['zScores_' patient]) = ZScores;
end


function [conflictModChan] = conflictModAnalysis(PowerTimeData, BandPower, Zscores, Trials_C, Trials_I, responseTimes)

    %%%%%%%%%%%%%%%%% Electrode Responsiveness Analysis %%%%%%%%%%%%%%%%%

    % assume Z score already cut: [onset time, avg RT], and all time uniform

    nConTrials = length(Trials_C);
    nInTrials = length(Trials_I);
    nChannels = size(Zscores{1}, 1);

    minDuration = 0.15; % duration in s to see if congruent z score > 1 for
    dt = mean(diff(PowerTimeData{1} (1,:)));   % assume uniform sampling
    minSamples = round(minDuration / dt); % Convert time duration to number of samples
    responsiveChannels = false(nChannels, 1);  % initialize logical array per channel

    for ch = 1:nChannels

        for tr_c = 1:nConTrials
            trialIdx = Trials_C(tr_c);
    
            % Logical vector where z-score for this trial/channel > 1
            above = Zscores{trialIdx}(ch, :) > 1;
    
            % Find all runs of true values, check if long enough
            d = diff([0, above, 0]); 
            startIdx = find(d == 1);
            endIdx   = find(d == -1) - 1;
            runLengths = endIdx - startIdx + 1;
            if any(runLengths >= minSamples)
                responsiveChannels(ch) = true;
                break;  % if yes, stop checking other trials for this channel
            end
        end
    end




    %%%%%%%%%%%%%%%%% Conflict Modulation Analysis %%%%%%%%%%%%%%%%%

    % Do permutation test for
    % [0s , 0s + time window] - stimulus align
    % [RT - time window , RT] - response align

    nChannels = size(BandPower{1}, 1);
    nTime = size(BandPower{1}, 2);
    t = PowerTimeData{1} (1,:);
    meanResponseTime = mean(responseTimes);
    resT = t - meanResponseTime;

    alpha = 0.05;
    nPermutations = 5000;
    conflictModChan = false(nChannels, 1); % final output per channel
    
    % Convert cell array to trial x time matrix per channel
    for ch = 1:nChannels
        
        % Build trial x time matrices
        stimConMatrix = zeros(nConTrials, nTime);
        stimInMatrix = zeros(nInTrials, nTime);
        resConMatrix = zeros(nConTrials, nTime);
        resInMatrix = zeros(nInTrials, nTime);
        
        for i = 1:nConTrials
            stimConMatrix(i, :) = BandPower{Trials_C(i)}(ch, :);

            rt = responseTimes(Trials_C(i));
            t_resp = t - rt; % re-align t to response at t = 0
            resConMatrix(i, :) = interp1(t_resp, BandPower{Trials_C(i)}(ch,:), resT, 'linear', NaN);

        end
        
        for i = 1:nInTrials
            stimInMatrix(i, :) = BandPower{Trials_I(i)}(ch, :);

            rt = responseTimes(Trials_I(i));
            t_resp = t - rt;
            resInMatrix(i, :) = interp1(t_resp, BandPower{Trials_I(i)}(ch,:), resT, 'linear', NaN);
   
        end

        timeWindow1 = find(t >= 0 & t <= meanResponseTime);
        timeWindow2 = find(resT >= -meanResponseTime & resT <= 0); % assume all NaNs got cut
            
        % Run permutation tests at each time bin
        p1 = NaN(1, nTime); % initialize p-values vector at each time bin
        for ti = timeWindow1
            p1(ti) = permutationTest(stimConMatrix(:, ti), stimInMatrix(:, ti), nPermutations);
        end

        p2 = NaN(1, nTime);
        for ti = timeWindow2
            p2(ti) = permutationTest(resConMatrix(:, ti), resInMatrix(:, ti), nPermutations);
        end
    
        % mask of significant p values: >=0.05 -> 1, <0.05 -> 0
        h1 = false(1, nTime);
        h2 = false(1, nTime);
        h1(p1 <= alpha) = true;
        h2(p2 <= alpha) = true;
    
        % Check for 15 consecutive significant bins (10 ms step between bins, want 150 ms)
        windowSize = 15;
        h1sum15 = movsum(h1, [windowSize - 1, 0]);  % moving sum of 15 consecutive bins
        h2sum15 = movsum(h2, [windowSize - 1, 0]);
        isConflictMod = any(h1sum15 >= windowSize) && any(h2sum15 >= windowSize);
    
        conflictModChan(ch) = isConflictMod;
    end


end


function decodeFeatures(features)
% concat_index = ((c-1)*nTrials)+ tr; % concatenate all electrodes



end






function [band_power_mean_max, band_power_zscore, normalized_band_power, power_time_data] = extractPowerFeatures(data, timeData, responseTimes)
    
    nTrials = length(data);
    nChannels = size(data{1}, 1);
            
    band_power_mean_max = cell(1, nTrials); % {trial} [column 1 = mean, column 2 = max]
    band_power_zscore = cell(1, nTrials);
    normalized_band_power = cell(1, nTrials); % {trial} [channels]

    for tr=1:nTrials
        band_power_mean_max{tr} = zeros(nChannels, 2);
        
        % normalized_band_power{tr} = zeros(nChannels, 2);

        % band_power_zscore{tr} = zeros(nChannels, 391); % hardcoded
        % length(S_epoched)=103?

        for ch=1:nChannels
            
            %%%%%%%%%%%%%%%%% Band Power Calculation %%%%%%%%%%%%%%%%%
            
            [Snorm,t,~,~]=computeNormalizedFreqMag(data{tr}(ch,:),1000);
            normalized_band_power {tr}(ch,:) = Snorm;  % all band power
            t = t + timeData{1}(1);  % shift t so its zero corresponds to 0 s in original time
            power_time_data {tr}(ch,:) = t;  % all band power
            
            % Cut time window: from  stimulus onset to behavioral response

            S_epoched_per_trial = Snorm(t >= 0 & t <= responseTimes(tr), :);   % time_indices after image onset
            meanResponseTime = mean(responseTimes);
            S_epoched_avgRT = Snorm(t >= 0 & t <= meanResponseTime, :);   % time_indices after image onset
            

            % [1st column - mean, 2nd column - max]
            band_power_mean_max{tr}(ch,1) = mean(S_epoched_per_trial(:));  % mean over time
            band_power_mean_max{tr}(ch,2) = max(S_epoched_per_trial(:));   % max over time


            %%%%%%%%%%%%% Baseline Band Power Calculation %%%%%%%%%%%%%
            
            % Calculate mean band power during 500 ms baseline

            S_baseline = Snorm(t >= -0.5 & t <= 0, :);

            mu = mean(S_baseline(:));
            sigma = std(S_baseline(:));

            band_power_zscore{tr} (ch, :) = (S_epoched_avgRT - mu) / sigma;  % Store in 3D matrix
        end
    end
end


function [y,t,S,f]=computeNormalizedFreqMag(x,Fs,params)

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

if strcmp(params.norm_option,'median')
    m=median(S,1);
elseif strcmp(params.norm_option,'mean')
    m=mean(S,1);
end

y=mean(S./repmat(m,size(S,1),1),2); % mean power for each trial and channel

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

% 3) not hardcoding the 391/103 --> just didn't preset the initialization

% 4) making window size (15) in conflictModAnalysis dynamic

%% Questions

% 1) filter order?

% 2) mismatched time vectors?

% 3) for responsive testing, ok to see if zscore>1 etc in ANY trial, to
% greenlight the electrode?

% 4) is the z score calculation correct.?