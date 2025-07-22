tic
files = {'BW42.mat', 'MG51b.mat'};
features = struct();

for i = 1:length(files)
    patientFeatures = processPatient(files{i});
    % patientFeatures = processPatientold(files{i});
    
    fNames = fieldnames(patientFeatures);
    for k = 1:numel(fNames)
        features.(fNames{k}) = patientFeatures.(fNames{k});
    end
end

save('features.mat','features');
clear
toc


%% functions

function [features] = processPatient(filename)
    load(fullfile(filename));
    [~, patient, ~] = fileparts(filename);
    
    congruentData = ft_data3.trial(Trials_C);
    incongruentData = ft_data3.trial(Trials_I);
    selectedChannels = ft_data3.label;
    mask1 = ~ismember(selectedChannels, ch_ictal);
    mask2 = ~strcmp(Parcellation_Sided_v3, 'RNan') & ~strcmp(Parcellation_Sided_v3, 'LNan');
    selectedChannels = find(mask1 & mask2);

    [conMeanPower, conMaxPower] = extractPowerFeatures(congruentData);
    [inMeanPower, inMaxPower] = extractPowerFeatures(incongruentData);
    
    features = struct();
    features.(['conMeanPower_' patient]) = conMeanPower;
    features.(['conMaxPower_' patient]) = conMaxPower;
    features.(['inMeanPower_' patient]) = inMeanPower;
    features.(['inMaxPower_' patient]) = inMaxPower;
    features.(['selectedChannels_' patient]) = selectedChannels;

end

function [meanPower, maxPower] = extractPowerFeatures(data)
    
    % Keep only if no switchChannels needed
    % nTrials = length(data);
    % nChannels = size(data{1}, 1);

    [data,nTrials,nChannels] = switchChannelsTrials(data);

    meanPower = cell(1, nChannels);
    maxPower = cell(1, nChannels);
    
    for c = 1:nChannels
    meanPower{c}=zeros(nTrials,1);
    maxPower{c}=zeros(nTrials,1);
    

    %%%%%%%%%%%%%%%%%% change to chronux soon %%%%%%%%%%%%%%%%


    meanPower{c} = mean(data{c}, 2);
    maxPower{c} = max(data{c},[],2);
    

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


    end
end

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


%%

% Notes
% 
% Potentially bad
% % BW42 channel 7 LA_03 LA_04
% % BW42 channel 15 LA_06 LA_07

% MG51b
% 
% plot(ft_data3.time{1},ft_data3.trial{1}(1:20,:))
%  plot(ft_data3.time{1},ft_data3.trial{1}(1:length(ft_data3.trial{1}(:,1)),:))
% % plot(ft_data3.time{1},ft_data3.trial{1}(features.selectedChannels_BW42,:))
% hold on
% legend(ft_data3.label);



%% unused old functions

function features = processPatientold(filename)
    load(fullfile(filename));
    [~, patient, ~] = fileparts(filename);
    
    congruentData = ft_data3.trial(Trials_C);
    incongruentData = ft_data3.trial(Trials_I);
    selectedChannels = ft_data3.label;
    mask1 = ~ismember(selectedChannels, ch_ictal);
    mask2 = ~strcmp(Parcellation_Sided_v3, 'RNan') & ~strcmp(Parcellation_Sided_v3, 'LNan');
    selectedChannels = find(mask1 & mask2);

    [conMeanPower, conMaxPower] = extractPowerFeaturesold(congruentData);
    [inMeanPower, inMaxPower] = extractPowerFeaturesold(incongruentData);
    
    features = struct();
    features.(['conMeanPower_' patient]) = conMeanPower;
    features.(['conMaxPower_' patient]) = conMaxPower;
    features.(['inMeanPower_' patient]) = inMeanPower;
    features.(['inMaxPower_' patient]) = inMaxPower;
    features.(['selectedChannels_' patient]) = selectedChannels;

end

function [meanPower, maxPower] = extractPowerFeaturesold(data)
    nTrials = length(data);
    nChannels = size(data{1}, 1);
    meanPower = cell(1, nChannels);
    maxPower = cell(1, nChannels);
    
    for c = 1:nChannels
        meanPowerVals = zeros(nTrials, 1);
        maxPowerVals = zeros(nTrials, 1);
        for t = 1:nTrials
            meanPowerVals(t) = mean(data{t}(c, :));
            maxPowerVals(t) = max(data{t}(c, :));
        end
        meanPower{c} = meanPowerVals;
        maxPower{c} = maxPowerVals;
    end
end
