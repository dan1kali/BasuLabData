

%% Preprocess and extract features

% use preProcessTheirData if their sub

tic
files = {'MG116', 'MG117', 'MG118', 'MG120', ...
         }; 
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
        % patientFeatures = preProcessTheirData(files{i});

        fNames = fieldnames(patientFeatures);
        for k = 1:numel(fNames)
            features.(fNames{k}) = patientFeatures.(fNames{k});
        end
    catch ME
        fprintf('Error processing patient %s: %s\n', files{i}, ME.message);
        continue;
    end
end
toc

save('features_squeeze5.mat','features','-v7.3');


%% Testing conflict mod on one patient (must preprocess and extract features first)

% use this if sub16:
% conflictModAnalysisTheirData

patientList = {'MG116', 'MG117', 'MG120'};  % BW42, MG51b, sub16

tic
for i = 1:length(patientList)
    try
        patient = patientList{i};
        fprintf('\nRunning conflictModAnalysis for patient: %s\n', patient);
       %         patient, ...
    
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

save('conflictModChans_squeeze5.mat', 'conflictModChans');
% save('conflictModChans_squeeze2.mat', 'conflictModChans','-v7.3');


%% functions

function [features] = preProcess(patient)
    load(fullfile('patientData',patient));
    % [~, patient, ~] = fileparts(filename);

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
    
    % Calculate max-min amplitudes
    amplitudes = zeros(nChannels, nTrials);
    for tr = 1:nTrials
        for ch = 1:nChannels
            % [-1 +2] window needed for artifact rejection
            timeIdx = ft_data3.time{tr} >= -1 & ft_data3.time{tr} <= 2;
            signal = ft_data3.trial{tr} (ch, timeIdx); % {trial} [channels x time]
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
        
    % convert rest of ft_data3
    ft_data_clean = ft_data3;
    ft_data_clean.trial = ft_data3_filt.trial(clean_trials_idx);
    ft_data_clean.time = ft_data3.time(clean_trials_idx);
    ft_data_clean.sampleinfo = ft_data3.sampleinfo(clean_trials_idx,:);
    responseTimes = TrialDet(clean_trials_idx,12);
    meanResponseTime = mean(responseTimes);
    
    % New congruent/incongruent indices
    [~, loc_C] = ismember(Trials_C, clean_trials_idx);
    [~, loc_I] = ismember(Trials_I, clean_trials_idx);
    Trials_C_clean = loc_C(loc_C > 0); % remove rejected indices
    Trials_I_clean = loc_I(loc_I > 0);

    fprintf('\nRejected %d of %d trials (%.2f%%)\n', ...
        sum(bad_trials), nTrials, 100 * sum(bad_trials) / nTrials);

    %%%%%%%%%%%%%%%%% Feature Extraction %%%%%%%%%%%%%%%%%
    
    timeData = ft_data_clean.time;
    selectedChannels = ft_data3.label;
    if exist('ch_ictal', 'var')
        mask1 = ~ismember(selectedChannels, ch_ictal);
        if exist('Parcellation_Sided_v3', 'var')
            mask2 = ~strcmp(Parcellation_Sided_v3, 'RNan') & ~strcmp(Parcellation_Sided_v3, 'LNan');
            selectedChannels = find(mask1 & mask2);
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

function [features] = preProcessTheirData(patient)
    % Load all .mat files for patient and merge into one struct
    files = dir(fullfile('patientData', patient, '*.mat'));
    allData = struct();
    cleanedData = struct();
    
    for k = 1:length(files)
        data = load(fullfile(files(k).folder, files(k).name));
        fn = fieldnames(data);
        for f = 1:length(fn)
            allData.(fn{f}) = data.(fn{f});
        end
    end

    patientResName = ['A01_gamma_number_' patient];
    patientStimName = ['A02_gamma_number_' patient];
    artName = [patient '_art_number'];
    
    nChannels = length(allData.(patientResName));
    nTrials = length(allData.congruency_number);

    cleanedData.vrt_number = (allData.vrt_number/1000)';
    allData.vrt_number = repmat(allData.vrt_number', 1, nChannels);
    
    cleanedData.mean_vrt = zeros(1, nChannels);

    cleanedData.trial = cell(1, nTrials);
    cleanedData.time = cell(1, nTrials);
    
    % use their time windowing --> (93:605)
    t_now = allData.t_trial(93:605);
    t_now = t_now-(t_now(1)-(-2.5)); % make the t_now(1) = -2.5, now response_time = 0; 

    cleanedData.time = repmat({t_now}, 1,nTrials);

    % Get rid of bad channels in all variables
    % selectedChannels = find(~cellfun('isempty', allData.(patientResName)));
    selectedChannels = allData.good_channel; % Use good channels, not just ones not empty
    allData.(patientResName) = allData.(patientResName)(selectedChannels);
    allData.(patientStimName)  = allData.(patientStimName)(selectedChannels);
    allData.(artName)         = allData.(artName)(selectedChannels);
    % allData.correct_number         = allData.correct_number(selectedChannels);
    nChannels = length(selectedChannels);  % Update channel count

    nTime = size(allData.(patientResName){1}, 2);
    incorrect_indices = find(allData.correct_number==0);
    nanRow = NaN(1, nTime);

    % Populate with NaNs for artifact trial rows that were eliminated
    for ch = 1:nChannels
        channelData = allData.(patientResName){ch};
        for i_cor = 1:length(incorrect_indices)
             cor_loc = incorrect_indices(i_cor);
             channelData = [channelData(1:cor_loc-1, :); nanRow; channelData(cor_loc:end, :)];
             allData.(patientResName){ch} = channelData;
        end
    end

    for ch = 1:nChannels
        art_indices = allData.(artName){ch};
        art_indices = setdiff(art_indices, incorrect_indices);
        channelData = allData.(patientResName){ch};
        for i_art = 1:length(art_indices)
            art_loc = art_indices(i_art);
            channelData = [channelData(1:art_loc-1, :); nanRow; channelData(art_loc:end, :)];
            allData.(patientResName){ch} = channelData;
        end
    end

    cong_vals = allData.congruency_number;
    cleanedData.Trials_C = find(cong_vals == 1);
    cleanedData.Trials_I = find(cong_vals == -1);

    [cleanedData.trial,~,~] = channel2trial(allData.(patientResName));
    cleanedData.channel = allData.(patientResName);

    %%%%%%%%%%%%%%%%%%%%%%%%%%% Feature Extraction %%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    responseTimes = cleanedData.vrt_number;
    timeData = cleanedData.time;
    Trials_C_clean = cleanedData.Trials_C;
    Trials_I_clean = cleanedData.Trials_I;

    % no need to turn into power bc their data is already power z scores
    % [PowerFeatures, PowerData, PowerTimeData] = extractPowerFeatures(cleanedData.trial, timeData, responseTimes);

    % A01_gamma_number_sub16{i_channel} = trials_matrix; 
    % trials_matrix = A01_gamma_number_sub16{i_channel};
 
    %%%%% features calculation moved here from extractPowerFeatures %%%%%
    band_power_mean_max = cell(1, nTrials); % {trial} [column 1 = mean, column 2 = max]

    for tr=1:nTrials
        band_power_mean_max{tr} = zeros(nChannels, 2);
  
        for ch=1:nChannels
         
            % ~~~~~~~~~~~~~~~~~~~~~~ Save features ~~~~~~~~~~~~~~~~~~~~~~
            
            % Cut time window: only look at stimulus onset to behavioral response
            t=timeData{1} (1,:);
            meanRT = nanmean(responseTimes);

            tWindow = t >= -(meanRT) & t <= 0; 
            S_epoched = cleanedData.channel{ch}(tr,tWindow);

            % got rid of dividing by mean across time (scaling):  
                % --> % of average power across time
                % y=mean(S,2); % mean power for each trial and channel

            % [1st column - mean, 2nd column - max]
            band_power_mean_max{tr}(ch,1) = mean(S_epoched(:));  % mean over time
            band_power_mean_max{tr}(ch,2) = max(S_epoched(:));   % max over time

        end
    end

    conPowerFeatures = band_power_mean_max(Trials_C_clean);
    inPowerFeatures  = band_power_mean_max(Trials_I_clean);

    features = struct();
    features.(['powerData_' patient]) = cleanedData.trial;
    features.(['powerTimeData_' patient]) = timeData;
    features.(['conPowerFeatures_' patient]) = conPowerFeatures;
    features.(['inPowerFeatures_' patient]) = inPowerFeatures;
    features.(['selectedChan_' patient]) = selectedChannels;
    features.(['responseTimes_' patient]) = responseTimes;
    features.(['trialsC_' patient]) = cleanedData.Trials_C;
    features.(['trialsI_' patient]) = cleanedData.Trials_I;

end

function [conflictModChans] = conflictModAnalysisTheirData(patient,TimeData, PowerData, Trials_C_clean, Trials_I_clean, responseTimes)

    %%%%%%%%%%%%%%%%% Electrode Responsiveness Analysis %%%%%%%%%%%%%%%%%

    nConTrials = length(Trials_C_clean);
    nInTrials = length(Trials_I_clean);
    nTrials = length(PowerData);

    %%%%%%% diferent form, channels is first dimension here
    % because PowerData isnt a cell array of cells anymore
    nChannels = size(PowerData{1}, 1); 
    meanRT = mean(responseTimes);

    minDuration = 0.15; % duration in s to see if congruent z score > 1 for
    dt = mean(diff(TimeData{1} (1,:)));   % assume uniform sampling
    minSamples = round(minDuration / dt); % Convert time duration to number of samples
    responsiveChannels = false(nChannels, 1);  % initialize logical array per channel

    % They already converted to Z scores
    band_power_zscore = PowerData;

    % t = PowerTimeData{1}(1,:); % assume all times uniform
    stimT = TimeData{1} (1,:); % just normal time array
    

    for ch = 1:nChannels

        for tr_c = 1:nConTrials
            trialIdx = Trials_C_clean(tr_c);
            
            % Check where z-score > 1 during window [onset time, avg RT]
            meanRT = nanmean(responseTimes);
            above = band_power_zscore{trialIdx}(ch, stimT >= -meanRT & stimT <= 0) > 1;
    
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
    
    band_power_mean_max = cell(1, nTrials); % {trial} [column 1 = mean, column 2 = max]

    for tr=1:nTrials

        band_power_mean_max{tr} = zeros(nChannels, 2);
        
        % for ch=1:nChannels
            
            % ~~~~~~~~~~~~~~~~~~~~~~ Save features ~~~~~~~~~~~~~~~~~~~~~~
            
            t = stimT;

            % Cut time window: only look at stimulus onset to behavioral response
            % tWindow = t >= 0 & t <= responseTimes(tr); 
            tWindow = t >= 0 & t <= meanRT; 

            S_epoched = band_power_zscore{tr}(:,tWindow);

            % [1st column - mean, 2nd column - max]
            band_power_mean_max{tr}(:,1) = mean(S_epoched,2);  % mean over time
            band_power_mean_max{tr}(:,2) = max(S_epoched,[],2);   % max over time

        % end
    end

    conPowerFeatures = band_power_mean_max(Trials_C_clean);
    inPowerFeatures  = band_power_mean_max(Trials_I_clean);



    %%%%%%%%%%%%%%%%% Conflict Modulation Analysis %%%%%%%%%%%%%%%%%

    % Do permutation test for:
    % [0s , 0s + time window] - stimulus align
    % [RT - time window , RT] - response align

    % [-RT , -RT + time window] - stimulus align
    % [0s - time window , 0s] - response align

    nTime = size(TimeData{1}(1,:),2);
    resT = TimeData{1} (1,:);
    stimT = resT + meanRT; % now 0 at 0-meanRT

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
        for i = 1:nConTrials % Build response aligned
            resConMatrix(i, :) = band_power_zscore{Trials_C_clean(i)}(ch, :);

            rt = responseTimes(Trials_C_clean(i)); % Build stimulus aligned
            t_stim = resT + rt; % add rt to resT t to get stimulus at t = 0
            stimConMatrix(i, :) = interp1(t_stim, band_power_zscore{Trials_C_clean(i)}(ch,:), stimT, 'linear', NaN);
        end
        
        % Repeat for incongruent
        for i = 1:nInTrials
            resInMatrix(i, :) = band_power_zscore{Trials_I_clean(i)}(ch, :);

            rt = responseTimes(Trials_I_clean(i));
            t_stim = resT + rt;
            stimInMatrix(i, :) = interp1(t_stim, band_power_zscore{Trials_I_clean(i)}(ch,:), stimT, 'linear', NaN);
        end

        timeWindowStim = find(stimT >= 0 & stimT <= (meanRT+0));
        plot_tws = find(stimT >= 0 & stimT <= (meanRT+0.5));

        %%%%% They use -(meanRT+0.5) to 0 for a "plotting" time
        timeWindowRes = find(resT >= -(meanRT + 0) & resT <= 0);
        plot_twr = find(resT >= -(meanRT + 0.5) & resT <= 0);
            
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

function [newdata, nTrials, nChannels] = channel2trial(data)
    % channel2trial
    % start with each data{channel} is [nTrials_c x timePoints]
    
    nChannels = length(data);
    
    % Find max number of trials across all channels
    maxTrials = 0;
    maxTimePoints = 0;
    for c = 1:nChannels
        sz = size(data{c});
        if ~isempty(sz)
            maxTrials = max(maxTrials, sz(1));
            if length(sz) > 1
                maxTimePoints = max(maxTimePoints, sz(2));
            end
        end
    end
    
    % Initialize newdata: each cell is {trial} [nChannels x timePoints]
    newdata = cell(1, maxTrials);
    
    for t = 1:maxTrials
        % Pre-allocate trial matrix with NaNs for all channels and time points
        trialMatrix = nan(nChannels, maxTimePoints);
        
        for c = 1:nChannels
            if size(data{c},1) >= t && ~isempty(data{c}(t,:))
                rowData = data{c}(t, :);
                % If shorter than maxTimePoints, fill only that part
                lenData = length(rowData);
                trialMatrix(c,1:lenData) = rowData;
            else
                % No data for this trial in this channel: leave NaNs
            end
        end
        
        newdata{t} = trialMatrix;
    end
    
    nTrials = maxTrials;
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

% 5) fix the preprocessing window, fix trial rejection for incorrect trials


%% Questions

% 1) filter order?

% 2) mismatched time vectors?

% 3) for responsive testing, ok to see if zscore>1 etc in ANY trial, to greenlight electrode?

% 4) is the z score calculation correct.?

%% plotting with p test on bottom

% subplot(2,1,1)
% plot(t(timeWindowStim), mean(stimConMatrix(:, timeWindowStim), 1))
% hold on
% plot(t(timeWindowStim), mean(stimInMatrix(:, timeWindowStim), 1))
% legend({'Congruent', 'Incongruent'})
% title(sprintf('Stimulus Aligned Window - Mean Across Trials - Channel %d Example', ch))
% 
% subplot(2,1,2)
% plot(t,p1)
% yline(0.1,'m--','alpha=0.10','LabelHorizontalAlignment', 'left')
% yline(0.05,'r--')
% title('p value')
% 
% 
% xline(t(200),'b--','dt=0.08s','LabelOrientation', 'horizontal')
% xline(t(208),'b--')
% 
% 
% 
% subplot(2,1,1)
% plot(resT(timeWindowRes), mean(resConMatrix(:, timeWindowRes), 1))
% hold on
% plot(resT(timeWindowRes), mean(resInMatrix(:, timeWindowRes), 1))
% legend({'Congruent', 'Incongruent'})
% title(sprintf('Response Aligned Window - Mean Across Trials - Channel %d Example', ch))
% 
% subplot(2,1,2)
% plot(resT,p2)
% yline(0.1,'m--','alpha=0.10','LabelHorizontalAlignment', 'left')
% yline(0.05,'r--')
% title('p value')
% 
% 
% xline(resT(228),'b--','dt=0.13s','LabelOrientation', 'horizontal')
% xline(resT(241),'b--')

%% with stdevs shaded

% % Stimulus-aligned plot (subplot 1)
% subplot(2,1,1)
% 
% % Calculate mean and SEM for Congruent/Incongruent
% nTrialsCon = size(stimConMatrix, 1);
% nTrialsIn  = size(stimInMatrix, 1);
% 
% meanCon = mean(stimConMatrix(:, timeWindowStim), 1);
% semCon  = std(stimConMatrix(:, timeWindowStim), [], 1) / sqrt(nTrialsCon);
% 
% meanIn = mean(stimInMatrix(:, timeWindowStim), 1);
% semIn  = std(stimInMatrix(:, timeWindowStim), [], 1) / sqrt(nTrialsIn);
% 
% % Plot shaded SEM areas
% fill([t(timeWindowStim), fliplr(t(timeWindowStim))], ...
%      [meanCon + semCon, fliplr(meanCon - semCon)], ...
%      [0.8 0.8 1], 'EdgeColor', 'none', 'FaceAlpha', 0.3); hold on
% fill([t(timeWindowStim), fliplr(t(timeWindowStim))], ...
%      [meanIn + semIn, fliplr(meanIn - semIn)], ...
%      [1 0.8 0.8], 'EdgeColor', 'none', 'FaceAlpha', 0.3);
% 
% % Plot mean lines
% plot(t(timeWindowStim), meanCon, 'b-', 'LineWidth', 1.5)
% plot(t(timeWindowStim), meanIn, 'r-', 'LineWidth', 1.5)
% legend({'Congruent ± SEM', 'Incongruent ± SEM', 'Congruent Mean', 'Incongruent Mean'})
% title(sprintf('Stimulus Aligned Window - Mean Across Trials - Channel %d Example', ch))
% 
% 
% 
% % Response-aligned plot (subplot 2)
% subplot(2,1,2)
% 
% % Calculate mean and SEM for Congruent/Incongruent
% nTrialsCon = size(resConMatrix, 1);
% nTrialsIn  = size(resInMatrix, 1);
% 
% meanCon = mean(resConMatrix(:, timeWindowRes), 1);
% semCon  = std(resConMatrix(:, timeWindowRes), [], 1) / sqrt(nTrialsCon);
% 
% meanIn = mean(resInMatrix(:, timeWindowRes), 1);
% semIn  = std(resInMatrix(:, timeWindowRes), [], 1) / sqrt(nTrialsIn);
% 
% % Plot shaded SEM areas
% fill([resT(timeWindowRes), fliplr(resT(timeWindowRes))], ...
%      [meanCon + semCon, fliplr(meanCon - semCon)], ...
%      [0.8 0.8 1], 'EdgeColor', 'none', 'FaceAlpha', 0.3); hold on
% fill([resT(timeWindowRes), fliplr(resT(timeWindowRes))], ...
%      [meanIn + semIn, fliplr(meanIn - semIn)], ...
%      [1 0.8 0.8], 'EdgeColor', 'none', 'FaceAlpha', 0.3);
% 
% % Plot mean lines
% plot(resT(timeWindowRes), meanCon, 'b-', 'LineWidth', 1.5)
% plot(resT(timeWindowRes), meanIn, 'r-', 'LineWidth', 1.5)
% legend({'Congruent ± SEM', 'Incongruent ± SEM', 'Congruent Mean', 'Incongruent Mean'})
% title(sprintf('Response Aligned Window - Mean Across Trials - Channel %d Example', ch))


%% all together

% tiledlayout(2,1, 'TileSpacing', 'compact', 'Padding', 'compact');
% 
% % ----------------------------
% nexttile  % Top-left (Stimulus Aligned)
% nConStim = size(stimConMatrix, 1);
% nInStim  = size(stimInMatrix, 1);
% 
% meanCon = mean(stimConMatrix(:, timeWindowStim), 1);
% semCon  = std(stimConMatrix(:, timeWindowStim), [], 1) / sqrt(nConStim);
% 
% meanIn = mean(stimInMatrix(:, timeWindowStim), 1);
% semIn  = std(stimInMatrix(:, timeWindowStim), [], 1) / sqrt(nInStim);
% 
% fill([t(timeWindowStim), fliplr(t(timeWindowStim))], [meanCon + semCon, fliplr(meanCon - semCon)], ...
%     [0.8 0.8 1], 'EdgeColor', 'none', 'FaceAlpha', 0.3); hold on
% fill([t(timeWindowStim), fliplr(t(timeWindowStim))], [meanIn + semIn, fliplr(meanIn - semIn)], ...
%     [1 0.8 0.8], 'EdgeColor', 'none', 'FaceAlpha', 0.3);
% plot(t(timeWindowStim), meanCon, 'b-', 'LineWidth', 1.5)
% plot(t(timeWindowStim), meanIn, 'r-', 'LineWidth', 1.5)
% legend({'Con ± SEM', 'In ± SEM', 'Con Mean', 'In Mean'})
% title(sprintf('Stimulus Aligned Window - Channel %d', ch), 'FontSize', 16)
% 
% % ----------------------------
% nexttile  % Top-right (Response Aligned)
% nConRes = size(resConMatrix, 1);
% nInRes  = size(resInMatrix, 1);
% 
% meanCon = mean(resConMatrix(:, timeWindowRes), 1);
% semCon  = std(resConMatrix(:, timeWindowRes), [], 1) / sqrt(nConRes);
% 
% meanIn = mean(resInMatrix(:, timeWindowRes), 1);
% semIn  = std(resInMatrix(:, timeWindowRes), [], 1) / sqrt(nInRes);
% 
% fill([resT(timeWindowRes), fliplr(resT(timeWindowRes))], [meanCon + semCon, fliplr(meanCon - semCon)], ...
%     [0.8 0.8 1], 'EdgeColor', 'none', 'FaceAlpha', 0.3); hold on
% fill([resT(timeWindowRes), fliplr(resT(timeWindowRes))], [meanIn + semIn, fliplr(meanIn - semIn)], ...
%     [1 0.8 0.8], 'EdgeColor', 'none', 'FaceAlpha', 0.3);
% plot(resT(timeWindowRes), meanCon, 'b-', 'LineWidth', 1.5)
% plot(resT(timeWindowRes), meanIn, 'r-', 'LineWidth', 1.5)
% title(sprintf('Response Aligned Window - Channel %d', ch), 'FontSize', 16)
% 
% % ----------------------------
% nexttile  % Bottom-left (Stimulus p-values)
% plot(t,p1)
% yline(0.1,'m--','alpha=0.10','LabelHorizontalAlignment', 'left')
% yline(0.05,'r--')
% title('p value')
% 
% % ----------------------------
% nexttile  % Bottom-right (Response p-values)
% plot(resT,p2)
% yline(0.1,'m--','alpha=0.10','LabelHorizontalAlignment', 'left')
% yline(0.05,'r--')
% title('p value')
% 
