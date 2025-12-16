%% Instructions

% Files:
% 'BW42', 'MG51b', 'MG79', 'MG86', ...
% 'MG89', 'MG90', 'MG91', 'MG95', ...
% 'MG96', 'MG99', 'MG102', 'MG104', ...
% 'MG105', 'MG106', 'MG111', 'MG112', ...
% 'MG116', 'MG117', 'MG118', 'MG120',...
% 'UCMC01', 'UCMC02', 'UCMC03', 'UCMC05', ...
% 'UCMC06', 'UCMC07', 'UCMC08', 'UCMC09', ...
% 'UCMC11', 'UCMC13', 'UCMC14', 'UCMC15', 'UCMC17',...



%% Extract Features
tic
files = {...
'BW42', 'MG51b', 'MG79', 'MG86', ...
'MG89', 'MG90', 'MG91', 'MG95', ...
'MG96', 'MG99', 'MG102', 'MG104', ...
'MG105', 'MG106', 'MG111', 'MG112',...
'MG116', 'MG117', 'MG118', 'MG120',...
'UCMC01', 'UCMC02', 'UCMC03', 'UCMC05', ...
'UCMC06', 'UCMC07', 'UCMC08', 'UCMC09', ...
'UCMC11', 'UCMC13', 'UCMC14', 'UCMC15', 'UCMC17',...
};

% pool = gcp('nocreate');
% if isempty(pool)
%     parpool(2);
% else
%     disp('Pool already running!');
% end

for i = 1:length(files)
    % try
        fprintf('\nRunning extractFeatures for patient: %s\n', files{i});
        extractFeatures(files{i});
    % catch ME
    %     fprintf('**** ERROR processing patient %s: %s *****\n', files{i}, ME.message);
    %     continue;
    % end
end
toc

%% functions

function extractFeatures(patient)

% patientPowerData form: nTrials x nChannels x nTime
% Changed all nChannels, dimensions accordingly

    inputPath = fullfile('c_inputPowerData', patient);
    outputName = patient;
    
    filesToLoad = {'normPow_Alpha.mat','channelROI.mat'...
        'normPow_Beta.mat', 'normPow_Gamma.mat', 'normPow_HighGamma.mat',...
        'normPow_Theta.mat', 'responseTime.mat', 'time_trial.mat',...
        'Trials_C_clean.mat', 'Trials_I_clean.mat'};
    
    for i = 1:length(filesToLoad)
        load(fullfile(inputPath, filesToLoad{i}));
    end


    PowerData = normPow_HighGamma;
    PowerTimeData = time_trial;

    [nTrials, nChannels, nTime] = size(PowerData);

    nConTrials = length(Trials_C_clean);
    nInTrials = length(Trials_I_clean);
    meanRT = mean(responseTime);

    res_band_power = cell(1, nTrials);
    band_power_mean_max = cell(1, nTrials); % {trial} [column 1 = mean, column 2 = max]

    t = PowerTimeData; % assume all times uniform


    %%%%%%%%%%%%%%%%% Channels for ROI %%%%%%%%%%%%%%%%%

    RegionsofInterest = {'dlPFC','dmPFC','OFC','vlPFC','STG','MTG','ITG','dACC','AMY','HIP'};
    ROIbyChannel = cell(length(RegionsofInterest),1);
    for ir = 1:length(RegionsofInterest)
        ROIbyChannel{ir} = find(strcmp(channelROI, RegionsofInterest{ir}))';
    end


    %%%%%%%%%%%%%%%%%%%%% Feature Extraction From baselined power %%%%%%%%%%%%%%%% 
    
    resT = t - meanRT;
    for tr=1:nTrials

        band_power_mean_max{tr} = zeros(nChannels, 3);
        res_band_power_mean_max{tr} = zeros(nChannels, 3);
        
        for ch=1:nChannels
            
            % ~~~~~~~~~~~~~~~~~~~~~~ Save features ~~~~~~~~~~~~~~~~~~~~~~
            
            % Cut time window: only look at stimulus onset to behavioral response
            % tWindow = t >= 0.1 & t <= responseTime(tr); 
            tWindow = t >= 0.1 & t <= meanRT; 
            S_epoched = (PowerData(tr,:,tWindow));

            % [1st column - mean, 2nd column - max]
            band_power_mean_max{tr}(:,1) = mean(S_epoched,3);  % mean over time
            band_power_mean_max{tr}(:,2) = max(S_epoched,[],3);   % max over time
            % band_power_mean_max{tr}(:,3) = trapz(tWindownew, S_epoched,2);


            % ~~~~~~~~~~~~~ Save response aligned features ~~~~~~~~~~~~~~~~~~
            % Build response aligned z scores
            rt = responseTime(tr); % Build response aligned
            t_resp = t - rt;  % re-align t to response at t = 0
            PowerData_reshaped = reshape(PowerData(tr,ch,:), 1, []); 
            res_band_power{tr}(ch,:) = interp1(t_resp, PowerData_reshaped, resT, 'linear', NaN);   
            S_res_epoched = res_band_power{tr} (ch,tWindow);


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

    outputFolder = fullfile('c_outputPowerData_log','highGamma',outputName);
    if ~exist(outputFolder, 'dir')
        mkdir(outputFolder);
    end

    save(fullfile(outputFolder, 'conPowerFeatures.mat'), 'conPowerFeatures');
    save(fullfile(outputFolder, 'inPowerFeatures.mat'), 'inPowerFeatures');
    save(fullfile(outputFolder, 'conResPowerFeatures.mat'), 'conResPowerFeatures');
    save(fullfile(outputFolder, 'inResPowerFeatures.mat'), 'inResPowerFeatures');
    save(fullfile(outputFolder, 'ROIbyChannel.mat'), 'ROIbyChannel');
end

%% resave data

% Folder containing the original .mat files
% inputFolder = '/Users/macbook/Desktop/research/BasuLabData/patientPowerDataOG'; 
% 
% % Folder where each patient's variable files will be saved
% outputFolder = '/Users/macbook/Desktop/research/BasuLabData/patientPowerData';
% 
% % Get list of all .mat files in the input folder
% matFiles = dir(fullfile(inputFolder, '*.mat'));
% 
% for i = 1:length(matFiles)
%     % Load the current .mat file
%     matFileName = fullfile(inputFolder, matFiles(i).name);
%     data = load(matFileName);
% 
%     % Get the patient code from the file name ("UC01.mat" -> "UC01")
%     [~, patientCode, ~] = fileparts(matFiles(i).name);
% 
%     % Create a folder for the patient 
%     patientFolder = fullfile(outputFolder, patientCode);
%     if ~exist(patientFolder, 'dir')
%         mkdir(patientFolder);
%     end
% 
%     % Save each variable in its own .mat file within the patient's folder
%     variableNames = fieldnames(data);
%     for j = 1:length(variableNames)
%         varName = variableNames{j};
%         varValue = data.(varName);
% 
%         % Save each variable to a separate .mat file named after the variable
%         saveFileName = fullfile(patientFolder, [varName, '.mat']);
%         eval([varName ' = data.' varName ';']);
%         save(saveFileName, varName);
%     end
% end
% 
% disp('All variables saved successfully!');