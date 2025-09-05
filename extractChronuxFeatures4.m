

%% Preprocess and extract features

tic
files = {'MG118'}; 
% 'BW42', 'MG51b', 'MG79', 'MG86', 
% 'MG89', 'MG90', 'MG91', 'MG95', 
% 'MG96', 'MG99', 'MG102', 'MG104', 
% 'MG105', 'MG106', 'MG111', 'MG112',
% 'MG116', 'MG117', 'MG118', 'MG120'

features = struct();

for i = 1:length(files)
    try
        fprintf('\nProcessing patient %s\n', files{i});
        patientFeatures = preProcess(files{i});

        fNames = fieldnames(patientFeatures);
        for k = 1:numel(fNames)
            features.(fNames{k}) = patientFeatures.(fNames{k});
        end
    catch MhE
        fprintf('Error processing patient %s: %s\n', files{i}, ME.message);
        continue;
    end
end
toc

save('features_MG118.mat','features');
% save('features_squeeze.mat','features','-v7.3');


%% Testing conflict mod on one patient (must preprocess and extract features first)

patientList = {'MG118'};  % BW42, MG51b, sub16

tic
for i = 1:length(patientList)
    try
        patient = patientList{i};
        fprintf('\nRunning conflictModAnalysis for patient: %s\n', patient);
    
        patientConflictModChans = conflictModAnalysis( ...
            patient, ...
            features.(['powerTimeData_' patient]), ...
            features.(['powerData_' patient]), ...
            features.(['trialsC_' patient]), ...
            features.(['trialsI_' patient]), ...
            features.(['responseTimes_' patient]));
    
        fNames = fieldnames(patientConflictModChans);
        for k = 1:numel(fNames)
            conflictModChans.(fNames{k}) = patientConflictModChans.(fNames{k});
        end

    catch ME
        fprintf('Error processing patient %s: %s\n', patientList{i}, ME.message);
        continue;
    end

end
toc

save('conflictModChans_MG118.mat', 'conflictModChans');
% save('conflictModChans_squeeze2.mat', 'conflictModChans','-v7.3');


%% functions

function [features] = preProcess(patient)
    load(fullfile('patientData', patient));
    % [~, patient, ~] = fileparts(filename);

    %%%%%%%%%%%%%%%%% Filtering %%%%%%%%%%%%%%%%%

    % got rid of preprocessing since already built in

    ft_data3_filt = ft_data3_filt_rs;
       
    nTrials = numel(ft_data3_filt.trial);
    nChannels = numel(ft_data3_filt.label);

    %%%%%%%%%%%%%%%%% Artifact Trial Rejection %%%%%%%%%%%%%%%%%
    
    % Calculate max-min amplitudes
    amplitudes = zeros(nChannels, nTrials);
    for tr = 1:nTrials
        for ch = 1:nChannels
            % [-1 +2] window needed for artifact rejection
            timeIdx = ft_data3_filt.time{tr} >= -1 & ft_data3_filt.time{tr} <= 2;
            signal = ft_data3_filt.trial{tr} (ch, timeIdx); % {trial} [channels x time]
            amplitudes(ch, tr) = max(signal) - min(signal);
        end
    end
    
    % Compute thresholds: mean + 5*std per channel
    mean_amp = mean(amplitudes, 2);           % [1 x channels]
    std_amp = std(amplitudes, 0, 2);          % [1 x channels]
    thresholds = mean_amp + 5 * std_amp;      % [1 x channels] - 5 stdevs
    
    % Identify artifact trials per channel
    artifact_mask = amplitudes > thresholds;      % [trials x channels]
    % bad_trials = any(artifact_mask, 1) | isnan(TrialDet(:,12))';
    bad_trials = any(artifact_mask, 1) | TrialDet(:,27)' ~= 1; % must = 1 for correct trials
    bad_trials = bad_trials';                     % convert to column [trials x 1]
    clean_trials_idx = find(~bad_trials);         % indices of good trials
        
    % convert rest of data
    ft_data_clean = ft_data3_filt;
    ft_data_clean.trial = ft_data3_filt.trial(clean_trials_idx);
    ft_data_clean.time = ft_data3_filt.time(clean_trials_idx);
    % ft_data_clean.sampleinfo = ft_data3_filt.sampleinfo(clean_trials_idx,:);
    responseTimes = TrialDet(clean_trials_idx,12);
    % meanResponseTime = mean(responseTimes);
    
    % New congruent/incongruent indices
    [~, loc_C] = ismember(Trials_C, clean_trials_idx);
    [~, loc_I] = ismember(Trials_I, clean_trials_idx);
    Trials_C_clean = loc_C(loc_C > 0); % remove rejected indices
    Trials_I_clean = loc_I(loc_I > 0);

    fprintf('\nRejected %d of %d trials (%.2f%%)\n', ...
        sum(bad_trials), nTrials, 100 * sum(bad_trials) / nTrials);

    %%%%%%%%%%%%%%%%% Feature Extraction %%%%%%%%%%%%%%%%%
    
    timeData = ft_data_clean.time;
    selectedChannels = ft_data3_filt.label;
    if exist('ch_ictal', 'var')
        mask1 = ~ismember(selectedChannels, ch_ictal);
        if exist('ParcellationValues_AllRegs', 'var')
            mask2 = ParcellationValues_AllRegs(:,1) ~= 1;
            selectedChannels = find(mask1 | mask2);
        else
            selectedChannels = find(mask1);
        end
    end
    [PowerFeatures, PowerData, PowerTimeData] = extractPowerFeatures3_1(ft_data_clean.trial, timeData, responseTimes);

    conPowerFeatures = PowerFeatures(Trials_C_clean);
    inPowerFeatures  = PowerFeatures(Trials_I_clean);

    features = struct();
    features.(['powerData_' patient]) = PowerData;
    features.(['powerTimeData_' patient]) = PowerTimeData;
    features.(['conPowerFeatures_' patient]) = conPowerFeatures;
    features.(['inPowerFeatures_' patient]) = inPowerFeatures;
    features.(['selectedChan_' patient]) = selectedChannels;
    features.(['responseTimes_' patient]) = responseTimes;
    features.(['trialsC_' patient]) = Trials_C_clean;
    features.(['trialsI_' patient]) = Trials_I_clean;

end

function [conflictModChans] = conflictModAnalysis(patient,PowerTimeData, PowerData, Trials_C_clean, Trials_I_clean, responseTimes)

% no more cell array, each trial of time x channels fed into spectrogram.
% output per cell: [freq x time × Nchannels]
% Changed all nChannels, dimensions accordingly

    nConTrials = length(Trials_C_clean);
    nInTrials = length(Trials_I_clean);
    nChannels = size(PowerData{1}, 3);
    nTrials = length(PowerData);
    meanRT = mean(responseTimes);

    minDuration = 0.15; % duration in s to see if congruent z score > 1 for
    % dt = mean(diff(PowerTimeData{1} (1,:)));   % assume uniform sampling
    dt = mean(diff(PowerTimeData{1} (1,:)));   % assume uniform sampling
    minSamples = round(minDuration / dt); % Convert time duration to number of samples
    responsiveChannels = false(nChannels, 1);  % initialize logical array per channel

    band_power_zscore = cell(1, nTrials);
    band_power_mean_max = cell(1, nTrials); % {trial} [column 1 = mean, column 2 = max]

    t = PowerTimeData{1}(1,:); % assume all times uniform


    %%%%%%%%%%%%%%%%% Electrode Responsiveness Analysis %%%%%%%%%%%%%%%%%

    for ch = 1:nChannels

        for tr=1:nTrials                        
            % Calculate mean band power during 500 ms baseline


            % ~~~~~~~~~~~~~~~~ Extract Baseline Power ~~~~~~~~~~~~~~~~~~~~
            
            % Cut time window to baseline: [-0.5s, 0]
            tBaselineWindow = t >= -0.5 & t <= 0; 
            % tBaselineWindow = t >= 1 & t <= 1.5; 
            S_baseline = PowerData{tr}(:,tBaselineWindow,ch);
            m = mean(S_baseline,2); % mean across time (1 value per frequency)
            expanded_m = repmat(m,1,size(S_baseline,2));

            % Power during baseline over time, per channel/trial
            S_baseline_norm = mean(S_baseline./expanded_m,1); 

            mu = mean(S_baseline_norm(:));
            sigma = std(S_baseline_norm(:));
    
            % ~~~~~~~~~~~~~~~~ Whole Signal Power ~~~~~~~~~~~~~~~~~~~~

            % Power during whole signal over time, per channel/trial
            
            % Try during [-0.5s onward] 
            % t_orig = t;
            % t = t(t>= -0.5);
            % t_idx = t_orig >= -0.5;
            

            % trying whole signal
            S_power = PowerData{tr}(:,:,ch);

            m = mean(S_power,2); 
            expanded_m = repmat(m,1,size(S_power,2));
            S_power_norm = mean(S_power./expanded_m,1); 

            % ~~~~~~~~~~~~~~~~~~~~ Calc Z Score ~~~~~~~~~~~~~~~~~~~~~~

            band_power_zscore{tr} (ch, :) = (S_power_norm - mu) / sigma;  % Store in 3D matrix

        end

        % ~~~~~~~~~~~~~~~~~~ Responsiveness Analysis ~~~~~~~~~~~~~~~~~~~~

        for tr_c = 1:nConTrials
            trialIdx = Trials_C_clean(tr_c);
            
            % Check where z-score > 1 during window [onset time, avg RT]
            meanRT = mean(responseTimes);
            above = band_power_zscore{trialIdx}(ch, t >= 0 & t <= meanRT) > 1;
    
            % Find all runs of true values, check if long enough
            d = diff([0, above, 0]); 
            startIdx = find(d == 1);
            endIdx   = find(d == -1) - 1;
            runLengths = endIdx - startIdx + 1;
            if any(runLengths >= minSamples)
                responsiveChannels(ch) = true;
                % break;  % if yes, stop checking other trials for this channel
            end
        end
    end


    %%%%%%%%%%%%%%%%%%%%% Feature Extraction From Z scores %%%%%%%%%%%%%%%% 

    for tr=1:nTrials

        band_power_mean_max{tr} = zeros(nChannels, 2);
        
        % for ch=1:nChannels
            
            % ~~~~~~~~~~~~~~~~~~~~~~ Save features ~~~~~~~~~~~~~~~~~~~~~~
            
            % Cut time window: only look at stimulus onset to behavioral response
            tWindow = t >= 0 & t <= responseTimes(tr); 

            S_epoched = band_power_zscore{tr}(:,tWindow);

            % [1st column - mean, 2nd column - max]
            band_power_mean_max{tr}(:,1) = mean(S_epoched, 2);  % mean over time
            band_power_mean_max{tr}(:,2) = max(S_epoched, [],2);   % max over time

        % end
    end

    conPowerFeatures = band_power_mean_max(Trials_C_clean);
    inPowerFeatures  = band_power_mean_max(Trials_I_clean);


    %%%%%%%%%%%%%%%%%%%%%% Conflict Modulation Analysis %%%%%%%%%%%%%%%%%

    % Do permutation test for:
    % [0s , 0s + time window] - stimulus align
    % [RT - time window , RT] - response align

    nChannels = size(PowerData{1}, 3);
    nTime = length(t);
    resT = t - meanRT;

    alpha = 0.1;
    nPermutations = 1000;
    conflictModChan = false(nChannels, 1); % final output per channel
    
    % Convert cell array to trial x time matrix per channel
    for ch = 1:nChannels
        
        stimConMatrix = zeros(nConTrials, nTime);
        stimInMatrix = zeros(nInTrials, nTime);
        resConMatrix = zeros(nConTrials, nTime);
        resInMatrix = zeros(nInTrials, nTime);
        
        % Congruent matrices
        for i = 1:nConTrials % Build stimulus aligned
            stimConMatrix(i, :) = band_power_zscore{Trials_C_clean(i)}(ch, :);

            rt = responseTimes(Trials_C_clean(i)); % Build response aligned
            t_resp = t - rt; % re-align t to response at t = 0
            resConMatrix(i, :) = interp1(t_resp, band_power_zscore{Trials_C_clean(i)}(ch,:), resT, 'linear', NaN);
        end
        
        % Repeat for incongruent
        for i = 1:nInTrials
            stimInMatrix(i, :) = band_power_zscore{Trials_I_clean(i)}(ch, :);

            rt = responseTimes(Trials_I_clean(i));
            t_resp = t - rt;
            resInMatrix(i, :) = interp1(t_resp, band_power_zscore{Trials_I_clean(i)}(ch,:), resT, 'linear', NaN);
        end

        timeWindowStim = find(t >= 0 & t <= meanRT);
        timeWindowRes = find(resT >= -meanRT & resT <= 0); % assume all NaNs will get cut
            
        % Run permutation tests at each time bin
        p1 = NaN(1, nTime); % initialize p-values vector at each time bin
        for ti = timeWindowStim
            p1(ti) = permutationTest(stimConMatrix(:, ti), stimInMatrix(:, ti), nPermutations);
        end

        p2 = NaN(1, nTime);
        for ti = timeWindowRes
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
        isConflictMod = any(h1sum15 >= windowSize) &&  any(h2sum15 >= windowSize);

        %%%%% can try changing to 100 ms --> 10 10ms bins %%%%
    
        conflictModChan(ch) = isConflictMod;
        conflictModChanIndices = find(conflictModChan==1);

    end
        
    finalChannelList = conflictModChan & responsiveChannels;
    fprintf('\n%d/%d channels (%.2f%%) are conflict modulated.\nWith alpha = %.2f, # permutations = %d\n', ...
    sum(finalChannelList), nChannels, 100 * sum(finalChannelList) / nChannels, alpha, nPermutations);
    
    conflictModChans = struct();
    conflictModChans.(['selectedChan_' patient '_confMod_a10']) = conflictModChanIndices;
    conflictModChans.(['conPowerFeatures_' patient]) = conPowerFeatures;
    conflictModChans.(['inPowerFeatures_' patient]) = inPowerFeatures;

end

function [band_power_mean_max, normalized_band_power, power_time_data] = extractPowerFeatures3_1(data, timeData, responseTimes)
    
% no more cell array, each trial of time x channels fed into spectrogram.
% output per cell: [freq x time × Nchannels]

    addpath(genpath('functions'))

    nTrials = length(data);
    nChannels = size(data{1}, 1);
            
    band_power_mean_max = cell(1, nTrials); % {trial} [column 1 = mean, column 2 = max]
    power_time_data = cell(1, nTrials);
    normalized_band_power = cell(1, nTrials); % {trial} [channels]

    for tr=1:nTrials
        band_power_mean_max{tr} = zeros(nChannels, 2);
        % normalized_band_power{tr} = cell(1, nChannels);
        % normalized_band_power{tr} = zeros(nChannels, 2);
            
        %%%%%%%%%%%%%%%%% Band Power Calculation %%%%%%%%%%%%%%%%%

        params = struct();
        % params.fpass = [70 110];  % frequency band you want to analyze

        [~,t,S,~]=computeNormalizedFreqMag(data{tr}',1000,params);
        
        % ~~~~~~~~~~~~ Save time-frequency power data ~~~~~~~~~~~~~~
        
        % Transpose each page: time x freq --> freq x time
        % Save all band power in frequency x time
        normalized_band_power {tr} = pagetranspose(S);
        t = t + timeData{1}(1);  % shift t
        power_time_data {tr} = t;

        % ~~~~~~~~~~~~~~~~~~~~~~ Save features ~~~~~~~~~~~~~~~~~~~~~~
        
        % Cut time window: only look at stimulus onset to behavioral response
        tWindow = t >= 0 & t <= responseTimes(tr); 
        S_epoched = S(tWindow,:,:);
        m = mean(S_epoched,1); % mean across time (1 value per frequency)
        expanded_m = repmat(m,size(S_epoched,1),1,1);

        % Power over time, per channel/trial
        % [time, freq, chan] --> (time, chan)
        S_norm_epoched = squeeze(mean(S_epoched./expanded_m,2)); 

        % [1st column - mean, 2nd column - max]
        band_power_mean_max{tr}(:,1,:) = mean(S_norm_epoched,1);  % mean over time
        band_power_mean_max{tr}(:,2,:) = max(S_norm_epoched,[],1);   % max over time

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
    m=mean(S,1); % collapses across time - 1 mean for each frequency bin
end

y=mean(S./repmat(m,size(S,1),1),2); % mean power for each trial and channel

% repmat is a vector of mean of S over time the size of frequency
% Divide by mean across time (scaling)  --> % of average power across time
% y=mean(S,2); % mean power for each trial and channel

end
