%% Instructions

% Assume your data is all in a subdirectory called "patientData"
% Must have following added to path: chronux_2_12

% Files:
% 'BW42', 'MG51b', 'MG79', 'MG86', ...
% 'MG89', 'MG90', 'MG91', 'MG95', ...
% 'MG96', 'MG99', 'MG102', 'MG104', ...
% 'MG105', 'MG106', 'MG111', 'MG112', ...
% 'MG116', 'MG117', 'MG118', 'MG120'

% 'UCMC01', 'UCMC02', 'UCMC03', 'UCMC05', ...
% 'UCMC06', 'UCMC07', 'UCMC08', 'UCMC09', ...
% 'UCMC11', 'UCMC13', 'UCMC14', 'UCMC15', 'UCMC17'


% 1) Preprocess for a given patient(s)
%   - Extracts power data, and generates features from that power data (not
%   z scores)
%   - saved in a directory called outputData

% 2) Then apply conflict modulation analysis
%   - this creates indices for conflict modulation
%   - and features from z-scored band power


%% 1) Preprocess

tic
files = {'UCMC15', 'UCMC17'};

% pool = gcp('nocreate');
% if isempty(pool)
%     parpool(4);
% else
%     disp('Pool already running!');
% end


for i = 1:length(files)
    try
        fprintf('\nProcessing patient %s\n', files{i});
        preProcess(files{i});
    catch ME
        fprintf('**** ERROR processing patient %s: %s *****\n', files{i}, ME.message);
        continue;
    end
end
toc

%% 2) Conflict Mod Analysis

tic
files = {'UCMC02', 'UCMC03', 'UCMC05', ...
'UCMC06', 'UCMC07', 'UCMC08', 'UCMC09', ...
'UCMC11', 'UCMC13', 'UCMC14'};

% pool = gcp('nocreate');
% if isempty(pool)
%     parpool(2);
% else
%     disp('Pool already running!');
% end

for i = 1:length(files)
    try
        fprintf('\nRunning conflictModAnalysis for patient: %s\n', files{i});
        conflictModAnalysis(files{i});
    catch ME
        fprintf('**** ERROR processing patient %s: %s *****\n', files{i}, ME.message);
        continue;
    end
end
toc

%% 3) Assign labels randomly

tic
files = {'MG91'};


pool = gcp('nocreate');
if isempty(pool)
    parpool(4);
else
    disp('Pool already running!');
end

parfor i = 1:length(files)
    try
        fprintf('\nRunning saveAllFeatures for patient: %s\n', files{i});
        saveAllFeatures(files{i});
    catch ME
        fprintf('**** ERROR processing patient %s: %s *****\n', files{i}, ME.message);
        continue;
    end
end
toc

%% functions

function preProcess(patient)
    load(fullfile('patientData', patient));
    outputName = patient;
    sr = 512;

    %%%%%%%%%%%%%%%%% Filtering %%%%%%%%%%%%%%%%%

    % Use pre pre-processed data

    ft_data3_filt = ft_data3_filt_rs;

    nTrials = numel(ft_data3_filt.trial);
    nChannels = numel(ft_data3_filt.label);

    %%%%%%%%%%%%%%%%% Artifact Trial Rejection %%%%%%%%%%%%%%%%%

    % Use pre-prepared column vector of bad trials
    if exist('TrialDet', 'var')
        bad_trials = unique([Trials_ZArt; find(isnan(TrialDet(:,12)))]); % convert to column [trials x 1]
    else
        bad_trials = Trials_ZArt;
    end
    all_trials = 1:nTrials;
    clean_trials_idx = setdiff(all_trials, bad_trials);  % indices of good trials

    % convert rest of data
    ft_data_clean = ft_data3_filt;
    ft_data_clean.trial = ft_data3_filt.trial(clean_trials_idx);
    ft_data_clean.time = ft_data3_filt.time(clean_trials_idx);
    responseTimes = TrialDet(clean_trials_idx,12);

    % New congruent/incongruent indices
    [~, loc_C] = ismember(Trials_C, clean_trials_idx);
    [~, loc_I] = ismember(Trials_I, clean_trials_idx);
    Trials_C_clean = loc_C(loc_C > 0); % remove rejected indices
    Trials_I_clean = loc_I(loc_I > 0);

    fprintf('\nRejected %d of %d trials (%.2f%%)\n', ...
        numel(bad_trials), nTrials, 100 * numel(bad_trials) / nTrials);


    %%%%%%%%%%%%%%%%% Channels for ROI %%%%%%%%%%%%%%%%%

    RegionsofInterest = {'dlPFC','dmPFC','OFC','vlPFC','STG','MTG','ITG','dACC','AMY','HIP'};
    ROIbyChannel = cell(length(RegionsofInterest),1);
    for ir = 1:length(RegionsofInterest)
        ROIbyChannel{ir} = find(strcmp(channelROI, RegionsofInterest{ir}))';
    end



    %%%%%%%%%%%%%%%%% Feature Extraction %%%%%%%%%%%%%%%%%

    timeData = ft_data_clean.time;
    selectedChans = ft_data3_filt.label;
    if exist('ch_ictal', 'var')
        mask1 = ~ismember(selectedChans, ch_ictal);
        if exist('ParcellationValues_AllRegs', 'var')
            mask2 = isnan(ParcellationValues_AllRegs(:,8))';
            selectedChans = find(mask1 | mask2);
        else
            selectedChans = find(mask1);
        end
    end
    [PowerFeatures, PowerData, PowerTimeData] = extractPowerFeatures(ft_data_clean, responseTimes);

    conBandPowerFeatures = PowerFeatures(Trials_C_clean);
    inBandPowerFeatures  = PowerFeatures(Trials_I_clean);

    outputFolder = fullfile('outputDataWaveletUC', outputName);
    if ~exist(outputFolder, 'dir')
        mkdir(outputFolder);
    end

    save(fullfile(outputFolder, 'powerData.mat'), 'PowerData','-v7.3');
    save(fullfile(outputFolder, 'powerTimeData.mat'), 'PowerTimeData');
    save(fullfile(outputFolder, 'conBandPowerFeatures.mat'), 'conBandPowerFeatures');
    save(fullfile(outputFolder, 'inBandPowerFeatures.mat'), 'inBandPowerFeatures');
    save(fullfile(outputFolder, 'selectedChans.mat'), 'selectedChans');
    save(fullfile(outputFolder, 'responseTimes.mat'), 'responseTimes');
    save(fullfile(outputFolder, 'trialsC.mat'), 'Trials_C_clean');
    save(fullfile(outputFolder, 'trialsI.mat'), 'Trials_I_clean');
    save(fullfile(outputFolder, 'ROIbyChannel.mat'), 'ROIbyChannel');

end


function saveAllFeatures(patient)

% no more cell array, each trial of time x channels fed into spectrogram.
% output per cell: [freq x time × Nchannels]
% Changed all nChannels, dimensions accordingly

    inputPath = fullfile('outputData', patient);
    outputName = patient;
    
    filesToLoad = {'powerData.mat', 'powerTimeData.mat', ...
        'selectedChans.mat','responseTimes.mat'};
    
    for i = 1:length(filesToLoad)
        load(fullfile(inputPath, filesToLoad{i}));
    end

    nChannels = size(PowerData{1}, 3);
    nTrials = length(PowerData);

    band_power_zscore = cell(1, nTrials);
    band_power_mean_max = cell(1, nTrials); % {trial} [column 1 = mean, column 2 = max]

    t = PowerTimeData{1}(1,:); % assume all times uniform

   %%%%%%%%%%%%%%%%% Electrode Responsiveness Analysis %%%%%%%%%%%%%%%%%

    for ch = 1:nChannels

        for tr=1:nTrials                        
            % Calculate mean band power during 500 ms baseline

            % ~~~~~~~~~~~~~~~~ Extract Baseline Power ~~~~~~~~~~~~~~~~~~~~
            
            % baseline extract: [-2s, +1s]
            % baseline cut: [-0.5s, 0]

            % t >= -0.5 & t <= 0;  % toggle to change window that is extracted
            tBaseline_extract = t >= -2 & t <= 1;
            tBaselinenew = t(tBaseline_extract);
            tBaseline_cut = tBaselinenew >= -0.5 & tBaselinenew <= 0;
            
            S_baseline = PowerData{tr}(:,tBaseline_extract,ch);

            m = mean(S_baseline,2); % mean across time (1 value per frequency)
            expanded_m = repmat(m,1,size(S_baseline,2));

            % Power during baseline over time, per channel/trial
            S_baseline_norm = mean(S_baseline./expanded_m,1); 
            S_baseline_cut = S_baseline_norm(tBaseline_cut);  %%%%%% final baseline %%%%%

            mu = mean(S_baseline_cut(:));
            sigma = std(S_baseline_cut(:));
    
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
            S_power_norm = mean(S_power./expanded_m,1);  %%%%%% final signal %%%%%

            % ~~~~~~~~~~~~~~~~~~~~ Calc Z Score ~~~~~~~~~~~~~~~~~~~~~~

            band_power_zscore{tr} (ch, :) = (S_power_norm - mu) / sigma;  % Store in 3D matrix

        end
    end




    %%%%%%%%%%%%%%%%%%%%% Feature Extraction From Z scores %%%%%%%%%%%%%%%% 

    for tr=1:nTrials

        band_power_mean_max{tr} = zeros(nChannels, 3);
        
        % for ch=1:nChannels
            
            % ~~~~~~~~~~~~~~~~~~~~~~ Save features ~~~~~~~~~~~~~~~~~~~~~~
            
            % Cut time window: only look at stimulus onset to behavioral response
            tWindow = t >= 0 & t <= responseTimes(tr); 
            % tWindow = t >= 0 & t <= nanmean(responseTimes); % MG91 only
            tWindownew = t(tWindow);
            S_epoched = band_power_zscore{tr}(:,tWindow);

            % [1st column - mean, 2nd column - max]
            band_power_mean_max{tr}(:,1) = mean(S_epoched, 2);  % mean over time
            band_power_mean_max{tr}(:,2) = max(S_epoched, [],2);   % max over time
            band_power_mean_max{tr}(:,3) = trapz(tWindownew, S_epoched,2);

        % end
    end

    allPowerFeatures = band_power_mean_max;

        
    outputFolder = fullfile('outputData', outputName);
    if ~exist(outputFolder, 'dir')
        mkdir(outputFolder);
    end

    save(fullfile(outputFolder, 'allPowerFeatures.mat'), 'allPowerFeatures');
end


function conflictModAnalysis(patient)

% no more cell array, each trial of time x channels fed into spectrogram.
% output per cell: [freq x time × Nchannels]
% Changed all nChannels, dimensions accordingly

    inputPath = fullfile('outputDataWaveletUC', patient);
    outputName = patient;
    
    filesToLoad = {'powerData.mat', 'powerTimeData.mat', ...
        'selectedChans.mat','responseTimes.mat', 'trialsC.mat', 'trialsI.mat'};
    
    for i = 1:length(filesToLoad)
        load(fullfile(inputPath, filesToLoad{i}));
    end

    nConTrials = length(Trials_C_clean);
    nInTrials = length(Trials_I_clean);
    nChannels = size(PowerData{1}, 3);
    nTrials = length(PowerData);
    nTime =  size(PowerData{1}, 2);
    meanRT = mean(responseTimes);

    minDuration = 0.15; % duration in s to see if congruent z score > 1 for
    % dt = mean(diff(PowerTimeData(1,:)));   % assume uniform sampling
    dt = mean(diff(PowerTimeData(1,:)));   % assume uniform sampling
    minSamples = round(minDuration / dt); % Convert time duration to number of samples
    responsiveChannels = false(nChannels, 1);  % initialize logical array per channel

    band_power_norm = cell(1, nTrials);
    res_band_power_zscore = cell(1, nTrials);
    band_power_mean_max = cell(1, nTrials); % {trial} [column 1 = mean, column 2 = max]

    t = PowerTimeData(1,:); % assume all times uniform




    % %%% ----- ANI BASELINE ------
    % [nTrials , nChannels, nFreqs, nTime]
    % power_time_data = ft_freq.time;
    % BLtimeIdx = find(power_time_data <= -0.1 & power_time_data >= -0.5);
    % TrialtimeIdx = find(power_time_data <= 1.4 & power_time_data >= 0.1);
    % 
    % PSD_baseline=mean(ft_freq.powspctrm(:,:,:,BLtimeIdx),4);
    % PSD_trial=ft_freq.powspctrm(:,:,:,TrialtimeIdx);
    % 
    % PSD_norm = mean(PSD_trial./repmat(PSD_baseline,1,1,1,numel(TrialtimeIdx)),3); % divide by mean BL, avg freq bin
    % PSD_norm = squeeze(log(PSD_norm)); % trl x chan x time
    % %%% ------------------------



    %%%%%%%%%%%%%%%%% Electrode Responsiveness Analysis %%%%%%%%%%%%%%%%%

    for ch = 1:nChannels

        for tr=1:nTrials                        
            % Calculate mean band power during 500 ms baseline
            % ~~~~~~~~~~~~~~~~ Extract Baseline Power ~~~~~~~~~~~~~~~~~~~~
            
            tBaseline_extract = t >= -0.5 & t <= -0.1;
            
            S_baseline = PowerData{tr}(:,tBaseline_extract,ch);
            S_baseline_meaned = mean(S_baseline,2);
            S_baseline_meaned_rep = repmat(S_baseline_meaned,1,nTime,1);
            % S_baseline_final = mean(S_baseline,1); % Average freq (1st dim)

            % mu = mean(S_baseline_final(:));
            % sigma = std(S_baseline_final(:));

            % ~~~~~~~~~~~~~~~~ Whole Signal Power ~~~~~~~~~~~~~~~~~~~~

            % Power during whole signal over time, per channel/trial
            % tSignal_extract = t >= 0.1 & t <= meanRT;
            % newT = t(tSignal_extract);

            S_power = PowerData{tr}(:,:,ch);
            S_power_final = mean(S_power,1); 

             % ~~~~~~~~~~~~~~~~~~~~ Calc Normalized Log Power ~~~~~~~~~~~~~~~~~

            band_power_norm{tr} (ch, :) = log(mean(S_power_final / S_baseline_meaned_rep),1);  % Store in 3D matrix
            % band_power_norm{tr} (ch, :) = (S_power_final - mu) / sigma;  % Store in 3D matrix


            % old
            % % ~~~~~~~~~~~~~~~~ Extract Baseline Power ~~~~~~~~~~~~~~~~~~~~
            % 
            % % baseline extract: [-2s, +1s]
            % % baseline cut: [-0.5s, 0]
            % 
            % tBaseline_extract = t >= -2 & t <= 1; % toggle to change window that is extracted
            % tBaselinenew = t(tBaseline_extract);
            % tBaseline_cut = tBaselinenew >= -0.5 & tBaselinenew <= 0;
            % 
            % S_baseline = PowerData{tr}(:,tBaseline_extract,ch);
            % 
            % m = mean(S_baseline,2); % mean across time (1 value per frequency)
            % expanded_m = repmat(m,1,size(S_baseline,2));
            % 
            % % Power during baseline over time, per channel/trial
            % S_baseline_norm = mean(S_baseline./expanded_m,1); 
            % S_baseline_cut = S_baseline_norm(tBaseline_cut); 
            % 
            % mu = mean(S_baseline_cut(:));
            % sigma = std(S_baseline_cut(:));
            %
            % % ~~~~~~~~~~~~~~~~ Whole Signal Power ~~~~~~~~~~~~~~~~~~~~
            %
            % % Power during whole signal over time, per channel/trial
            % 
            % % Try during [-0.5s onward] 
            % % t_orig = t;
            % % t = t(t>= -0.5);
            % % t_idx = t_orig >= -0.5;
            % 
            % % trying whole signal
            % S_power = PowerData{tr}(:,:,ch);
            % 
            % m = mean(S_power,2); 
            % expanded_m = repmat(m,1,size(S_power,2));
            % S_power_norm = mean(S_power./expanded_m,1); 
            %
            % % ~~~~~~~~~~~~~~~~~~~~ Calc Z Score ~~~~~~~~~~~~~~~~~~~~~~
            %
            % band_power_zscore{tr} (ch, :) = (S_power_final - mu) / sigma;  % Store in 3D matrix

        end



        % ~~~~~~~~~~~~~~~~~~ Responsiveness Analysis ~~~~~~~~~~~~~~~~~~~~

        % for tr_c = 1:nConTrials
        %     trialIdx = Trials_C_clean(tr_c);
        % 
        %     % Check where z-score > 1 during window [onset time, avg RT]
        %     meanRT = mean(responseTimes);
        %     above = band_power_norm{trialIdx}(ch, t >= 0 & t <= meanRT) > 1;
        % 
        %     % Find all runs of true values, check if long enough
        %     d = diff([0, above, 0]); 
        %     startIdx = find(d == 1);
        %     endIdx   = find(d == -1) - 1;
        %     runLengths = endIdx - startIdx + 1;
        %     if any(runLengths >= minSamples)
        %         responsiveChannels(ch) = true;
        %         % break;  % if yes, stop checking other trials for this channel
        %     end
        % end
    end



    %%%%%%%%%%%%%%%%%%%%% Feature Extraction From baselined power %%%%%%%%%%%%%%%% 
    
    resT = t - meanRT;
    for tr=1:nTrials

        band_power_mean_max{tr} = zeros(nChannels, 3);
        res_band_power_mean_max{tr} = zeros(nChannels, 3);
        
        for ch=1:nChannels
            
            % ~~~~~~~~~~~~~~~~~~~~~~ Save features ~~~~~~~~~~~~~~~~~~~~~~
            
            % Cut time window: only look at stimulus onset to behavioral response
            tWindow = t >= 0.1 & t <= responseTimes(tr); 
            tWindownew = t(tWindow);
            % S_epoched = band_power_norm{tr}(:,tWindow);
            S_epoched = band_power_norm{tr}(:,:);

            % [1st column - mean, 2nd column - max]
            band_power_mean_max{tr}(:,1) = mean(S_epoched, 2);  % mean over time
            band_power_mean_max{tr}(:,2) = max(S_epoched, [],2);   % max over time
            % band_power_mean_max{tr}(:,3) = trapz(tWindownew, S_epoched,2);


            % ~~~~~~~~~~~~~ Save response aligned features ~~~~~~~~~~~~~~~~~~
            % Build response aligned z scores
            rt = responseTimes(tr); % Build response aligned
            t_resp = t - rt;  % re-align t to response at t = 0
            res_band_power_zscore{tr}(ch,:) = interp1(t_resp, band_power_norm{tr}(ch,:), resT, 'linear', NaN);   
            S_res_epoched = res_band_power_zscore{tr}(:,:);

            % [1st column - mean, 2nd column - max]
            res_band_power_mean_max{tr}(:,1) = mean(S_res_epoched, 2);  % mean over time
            res_band_power_mean_max{tr}(:,2) = max(S_res_epoched, [],2);   % max over time
            % res_band_power_mean_max{tr}(:,3) = trapz(tWindownew, S_epoched,2);

        end
    end
    % normal features
    conPowerFeatures = band_power_mean_max(Trials_C_clean);
    inPowerFeatures  = band_power_mean_max(Trials_I_clean);

    % response aligned
    conResPowerFeatures = band_power_mean_max(Trials_C_clean);
    inResPowerFeatures  = band_power_mean_max(Trials_I_clean);


    %%%%%%%%%%%%%%%%%%%%%% Conflict Modulation Analysis %%%%%%%%%%%%%%%%%

    % Do permutation test for:
    % [0s , 0s + time window] - stimulus align
    % [RT - time window , RT] - response align

    nChannels = size(PowerData{1}, 3);
    nTime = length(t);

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
            stimConMatrix(i, :) = band_power_norm{Trials_C_clean(i)}(ch, :);

            rt = responseTimes(Trials_C_clean(i)); % Build response aligned
            t_resp = t - rt; % re-align t to response at t = 0
            resConMatrix(i, :) = interp1(t_resp, band_power_norm{Trials_C_clean(i)}(ch,:), resT, 'linear', NaN);
        end

        % Repeat for incongruent
        for i = 1:nInTrials
            stimInMatrix(i, :) = band_power_norm{Trials_I_clean(i)}(ch, :);

            rt = responseTimes(Trials_I_clean(i));
            t_resp = t - rt;
            resInMatrix(i, :) = interp1(t_resp, band_power_norm{Trials_I_clean(i)}(ch,:), resT, 'linear', NaN);
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
        conflictModChans = find(conflictModChan==1);
        conflictModChans = conflictModChans(ismember(conflictModChans, selectedChans));

    end

    finalChannelList = conflictModChan & responsiveChannels;
    fprintf('\n%d/%d channels (%.2f%%) are conflict modulated.\nWith alpha = %.2f, # permutations = %d\n', ...
    sum(finalChannelList), nChannels, 100 * sum(finalChannelList) / nChannels, alpha, nPermutations);

    outputFolder = fullfile('outputDataWaveletUC', outputName);
    if ~exist(outputFolder, 'dir')
        mkdir(outputFolder);
    end

    save(fullfile(outputFolder, 'conflictModChans.mat'), 'conflictModChans');
    save(fullfile(outputFolder, 'powerTimeData.mat'), 'PowerTimeData');
    save(fullfile(outputFolder, 'conPowerFeatures.mat'), 'conPowerFeatures');
    save(fullfile(outputFolder, 'inPowerFeatures.mat'), 'inPowerFeatures');
    save(fullfile(outputFolder, 'conResPowerFeatures.mat'), 'conResPowerFeatures');
    save(fullfile(outputFolder, 'inResPowerFeatures.mat'), 'inResPowerFeatures');
end


function [band_power_mean_max, band_power, power_time_data] = extractPowerFeatures(data, responseTimes)
    
    % Use wavelet decomp

    % Calculate PSD 
    cfg = [];
    cfg.output = 'pow';
    cfg.keeptrials = 'yes';
    cfg.method = 'wavelet';
    cfg.pad='nextpow2'; %pads data to match desired freq resolution (https://www.fieldtriptoolbox.org/faq/spectral/freqanalysis_paddinginsufficient/)
    cfg.foi = [70:120]; %1-110Hz w/1Hz steps
    cfg.toi = [-1:(1/512*4):2]; % -1s to +2s w.r.t. image onset, w/7.8ms steps. If you leave sampling rate at 1000Hz, would just do [-1:0.01:2] instead. 
    
    % 4D Matrix: [trials × channels × frequencies × time]
    ft_freq = ft_freqanalysis(cfg,data);
    
    power_time_data = ft_freq.time;
    meanRT = mean(responseTimes);
    [nTrials , nChannels, nFreqs, nTime] = size(ft_freq.powspctrm);
    band_power_mean_max = cell(1, nTrials); % {trial} [column 1 = mean, column 2 = max]
    band_power = cell(1, nTrials); % {trial} [channels]

    for tr=1:nTrials
        band_power_mean_max{tr} = zeros(nChannels, 2);
        
        % Extract and squeeze to get [channel x freq x time]
        trial_data = squeeze(ft_freq.powspctrm(tr, :, :, :));  
        
        % Rearrange to [freq x time x channel]
        band_power{tr} = permute(trial_data, [2, 3, 1]); 

        % ~~~~~~~~~~~~~~~~~~~ Save non z-scored features ~~~~~~~~~~~~~~~~~~~        
        tWindow_cut = power_time_data >= 0.1 & power_time_data <= meanRT; % stim onset to behavioral response

        % Average across frequencies
        band_power_cut = band_power{tr}(:,tWindow_cut,:); 

        % Mean across time (dim 2): result is [freq x 1 x channels]
        final_band_power_cut = squeeze(mean(band_power_cut, 1)); % can change to normalized_band_power_cut
        
        % [1st column - mean, 2nd column - max] over time
        band_power_mean_max{tr}(:,1) = mean(final_band_power_cut,1);  % mean over time
        band_power_mean_max{tr}(:,2) = max(final_band_power_cut,[],1);   % max over time

    end
    
   
    
    
    % old

    % power_time_data = ft_freq.time;
    % 
    % [nTrials , nChannels, nFreqs, nTime] = size(ft_freq.powspctrm);
    % band_power_mean_max = cell(1, nTrials); % {trial} [column 1 = mean, column 2 = max]
    % normalized_band_power = cell(1, nTrials); % {trial} [channels]
    % 
    % for tr=1:nTrials
    %     band_power_mean_max{tr} = zeros(nChannels, 2);
    % 
    %     % Extract and squeeze to get [channel x freq x time]
    %     trial_data = squeeze(ft_freq.powspctrm(tr, :, :, :));  
    % 
    %     % Rearrange to [freq x time x channel]
    %     normalized_band_power{tr} = permute(trial_data, [2, 3, 1]); 
    % 
    %     % ~~~~~~~~~~~~~~~~~~~ Save non z-scored features ~~~~~~~~~~~~~~~~~~~        
    %     tWindow_cut = power_time_data >= 0 & power_time_data <= responseTimes(tr); % stim onset to behavioral response
    % 
    %     band_power_whole = normalized_band_power{tr}(:,:,:); % [freq x time x channels]
    %     band_power_cut = normalized_band_power{tr}(:,tWindow_cut,:); 
    % 
    %     % Mean across time (dim 2): result is [freq x 1 x channels]
    %     m = mean(band_power_whole, 2); % can change to normalized_band_power_cut
    % 
    %     % Expand to match original size: [freq x time x channels]
    %     expanded_m = repmat(m, 1, length(find(tWindow_cut)), 1);
    % 
    %     % [time x channels]
    %     final_band_power_cut = squeeze(mean(band_power_cut./expanded_m,1)); 
    % 
    %     % normalized_band_power_cut = normalized_band_power{tr}(:,:,tWindow_cut);
    % 
    %     % [1st column - mean, 2nd column - max] over time
    %     band_power_mean_max{tr}(:,1) = mean(final_band_power_cut,1);  % mean over time
    %     band_power_mean_max{tr}(:,2) = max(final_band_power_cut,[],1);   % max over time
    % 
    % end
end