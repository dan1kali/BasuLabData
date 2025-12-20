
% Assume your data is all in a subdirectory called "patientData"
% Must have functions added to path: generateCrossValInd.m , permutationTest.m

% Files:
% 'BW42', 'MG51b', 'MG79', 'MG86', ...
% 'MG89', 'MG90', 'MG91', 'MG95', ...
% 'MG96', 'MG99', 'MG102', 'MG104', ...
% 'MG105', 'MG106', 'MG111', 'MG112',...
% 'MG116', 'MG117', 'MG118', 'MG120',...
% 'UCMC01', 'UCMC02', 'UCMC03', 'UCMC05', ...
% 'UCMC06', 'UCMC07', 'UCMC08', 'UCMC09', ...
% 'UCMC11', 'UCMC13', 'UCMC14', 'UCMC15', 'UCMC17',...


% 'selChans - z-scored power','confChans - z-scored power', ...
% 'selChans - normalized power',  'confChans - normalized power'
% 'selChans - non z-scored power','confChans - non z-scored power' ...
% 'selChans - z-scored power - rand','confChans - z-scored power - rand', ...
% 'selChans - z-scored power - res','confChans - z-scored power - res', ...
% 'all chans - normalized power',...
% 'region chans - normalized power',..


%% Plot Patient Data

subjects = {... 
'BW42', 'MG51b', 'MG79', 'MG86', ...
'MG89', 'MG90', 'MG91', 'MG95', ...
'MG96', 'MG99', 'MG102', 'MG104', ...
'MG105', 'MG106', 'MG111', 'MG112',...
'MG116', 'MG117', 'MG118', 'MG120',...
'UCMC01', 'UCMC02', 'UCMC03', 'UCMC05', ...
'UCMC06', 'UCMC07', 'UCMC08', 'UCMC09', ...
'UCMC11', 'UCMC13', 'UCMC14', 'UCMC15', 'UCMC17',...
};

config = {'allChans - normalized power',...
};


nBars = length(subjects);
% nBars = 1;
nGroups = length(config);
groupedBars = zeros(nBars,nGroups);
groupedErr = zeros(nBars,nGroups);
shapVals = cell(nBars,nGroups);
meanBetas = cell(nBars,nGroups);
maxBetas = cell(nBars,nGroups);

for igroup = 1:nGroups
    % try
    if nBars ==1
        [y, err,shap,mean_weights,max_weights,sel_chan_number]  = SVM(subjects(:),config(igroup));
        groupedBars(:,igroup) = y;
        groupedErr(:,igroup) = err;
        meanBetas = mean_weights;
        maxBetas = max_weights;
    else
        for ibar=1:nBars
            [y, err,shap,mean_weights,max_weights,sel_chan_number]  = SVM(subjects(ibar),config(igroup));
            groupedBars(ibar,igroup) = y; % 33 patients x 1, each 1 value is mean of 500 values
            groupedErr(ibar,igroup) = err;
            shapVals{ibar,igroup} = shap;
            meanBetas{ibar,igroup} = mean_weights;
            maxBetas{ibar,igroup} = max_weights;
        end
    end
    % catch ME
    % disp(getReport(ME, 'extended'));
    % end
end

%% plot

barplot(groupedBars, subject, groupedErr, config)
% weightsplot(meanBetas,maxBetas,sel_chan_number,subjects) % only most recent sub and condition

%%

% % initializing
% for i = 1:33
%     groupedRegionAccuracies{i} = zeros(7,2);
% end

% % making values the same in certain columns
% for i = 1:33
%     groupedRegionAccuracies{i}(1,:) = regionAccuracies{i}(1,:);
% end


for i = 1:33
    groupedRegionAccuracies{i}(2,1) = groupedBars(i);
    groupedRegionAccuracies{i}(2,2) = groupedErr(i);
end

save('groupedRegionAccuracies_wavelet_nolog_theta.mat','groupedRegionAccuracies')
% save('regionAccuracies_wavelet_nolog_theta.mat','regionAccuracies')


%% Plot number of channels/min number of Trials

subjects = {'BW42', 'MG51b', 'MG79', 'MG86', ...
'MG89', 'MG90',  'MG91', 'MG95', ...
'MG96', 'MG99', 'MG102', 'MG104', ...
'MG105', 'MG106', 'MG111', 'MG112',...
'MG116', 'MG117', 'MG118', 'MG120'};

nChannels = nTrialsMin(subjects);

barplot(nChannels, subjects);
ylim([0 175]); xlim([0 21]); title('Min Number of Trials per Patient')


%% functions

function [fea_number_con, fea_number_in, m_number_out,sel_chan_number,n] = concatenateFeatures(m_number, config,inputPath)
    
    % inputPath = fullfile('outputDataChronux_zscore',subject);
    % inputPath = fullfile('outputPowerData_nolog','highGamma',subject);
    
    filesToLoad = { ...
        'conPowerFeatures.mat','inPowerFeatures.mat',...
                'ROIbyChannel.mat',...
        };        

        % 'allPowerFeatures.mat', ...
        % 'taskChans.mat', ...
        % 'ROIbyChannel.mat',...
        % 'conLogPowerFeatures.mat','inLogPowerFeatures.mat', ...
        % 'conBandPowerFeatures.mat', 'inBandPowerFeatures.mat', ...
        % 'selectedChans.mat','conflictModChans.mat', ...
        % 'conResPowerFeatures.mat','inResPowerFeatures.mat',...

    for i = 1:length(filesToLoad)
        load(fullfile(inputPath, filesToLoad{i}));
    end

    % nTrials = numel(allPowerFeatures); 
    % randomOrder = randperm(nTrials);
    % half = floor(nTrials/2);
    % conPowerFeaturesRand = allPowerFeatures(randomOrder(1:half));
    % inPowerFeaturesRand = allPowerFeatures(randomOrder(half+1:end));


     switch config{1}
         case 'confChans - normalized power'
            sel_chan_number = conflictModChans;
            conPower = conPowerFeatures;
            inPower = inPowerFeatures;
         case 'selChans - normalized power'
            sel_chan_number = selectedChans;
            conPower = conPowerFeatures;
            inPower = inPowerFeatures;

         % % Log power only (non response aligned)
         % case 'confChans - normalized log power' % only log for UC data
         %    sel_chan_number = conflictModChans;
         %    conPower = conPowerFeatures;
         %    inPower = inPowerFeatures;
         % case 'selChans - normalized log power'
         %    sel_chan_number = selectedChans;
         %    conPower = conPowerFeatures;
         %    inPower = inPowerFeatures;

         % Random labels
         case 'confChans - z-scored power - rand'
            sel_chan_number = conflictModChans;
            conPower = conPowerFeaturesRand;
            inPower = inPowerFeaturesRand;
         case 'selChans - z-scored power - rand'
            sel_chan_number = selectedChans;
            conPower = conPowerFeaturesRand;
            inPower = inPowerFeaturesRand;

         % Response aligned
         case 'confChans - z-scored power - res'
            sel_chan_number = conflictModChans;
            conPower = conResPowerFeatures;
            inPower = inResPowerFeatures;
         case 'selChans - z-scored power - res'
            sel_chan_number = selectedChans;
            conPower = conResPowerFeatures;
            inPower = inResPowerFeatures;

         % Band power only 
         case 'confChans - non z-scored power'
            sel_chan_number = conflictModChans;
            conPower = conBandPowerFeatures;
            inPower = inBandPowerFeatures;
         case 'selChans - non z-scored power'
            sel_chan_number = selectedChans;
            conPower = conBandPowerFeatures;
            inPower = inBandPowerFeatures;
         
         % all chans 
         case 'allChans - normalized power'
            sel_chan_number = 'all';
            conPower = conPowerFeatures;
            inPower = inPowerFeatures;

         % all chans 
         case 'region chans - normalized power'
            RegionLabels = {'dlPFC', 'dmPFC', 'OFC', 'vlPFC', 'STG', 'MTG', 'ITG', 'dACC', 'AMY', 'HIP'};
            ROI = {'dlPFC', 'dmPFC', 'OFC', 'vlPFC', 'STG', 'MTG', 'ITG', 'dACC', 'AMY', 'HIP'}; %%%%% change to region you want to graph
            ROIdx = ismember(RegionLabels, ROI);
            ROI_sel_chan_number = [ROIbyChannel{ROIdx}];
            
            sel_chan_number = ROI_sel_chan_number;
            conPower = conPowerFeatures;
            inPower = inPowerFeatures;
        
          % task chans 
         case 'taskChans - normalized power'
            sel_chan_number = taskChans;
            conPower = conPowerFeatures;
            inPower = inPowerFeatures;

         otherwise
            error('Unknown config: %s\n', config{1});
    end

    fea_number_con = [];
    fea_number_in = [];

    if strcmp(sel_chan_number, 'all')
        sel_chan_number = 1:length(conPowerFeatures{1});
    end

    n = min(length(conPower), length(inPower));

    % n=40;

    if ~isempty(sel_chan_number)


        for i = 1:length(sel_chan_number)
            %%%%%%%%%%%%%% move randsample outside the loop %%%%%%%%%%%%%%%
             
            ch = sel_chan_number(i); 
            % --- Max Power: Pull from Column 2 ---
            conMaxVals = cellfun(@(x) x(ch,2), conPower);  % con max across trials
            inMaxVals = cellfun(@(x) x(ch,2), inPower);   % in max across trials
            % --- Mean Power: Pull from Column 1 ---
            conMeanVals = cellfun(@(x) x(ch,1), conPower);  % con mean across trials
            inMeanVals = cellfun(@(x) x(ch,1), inPower);   % in mean across trials

            fea_number_con(1:n, m_number) = randsample(conMaxVals, n);
            fea_number_in(1:n, m_number) = randsample(inMaxVals, n);
            m_number = m_number + 1;
            fea_number_con(1:n, m_number) = randsample(conMeanVals, n);
            fea_number_in(1:n, m_number) = randsample(inMeanVals, n);
            m_number = m_number + 1;

        end

        m_number_out = m_number;
    else
        m_number_out = 0;
    end
end

function nChannels = numberChannels(subjects)
nChannels = zeros(length(subjects));

    for i_sub = 1:length(subjects)
        inputPath = fullfile('outputData', subjects);
        
        filesToLoad = {'conPowerFeatures.mat'};
        
        for i = 1:length(filesToLoad)
            load(fullfile(inputPath{i_sub}, filesToLoad{i}));
        end

    nChannels(i_sub) = size(conPowerFeatures{1},1);
    end

end

function nTrialsMin = nTrialsMin(subjects)
nTrialsMin = zeros(length(subjects));

    for i_sub = 1:length(subjects)
        inputPath = fullfile('outputData', subjects);
        
        filesToLoad = {'conPowerFeatures.mat','inPowerFeatures.mat'};
        
        for i = 1:length(filesToLoad)
            load(fullfile(inputPath{i_sub}, filesToLoad{i}));
        end

    nConTrials = numel(conPowerFeatures);
    nInTrials = numel(inPowerFeatures);

    nTrialsMin(i_sub) = min(nConTrials, nInTrials);
    end

end

function [y,err,shap,mean_weights,max_weights,sel_chan_number] = SVM(subjects,config)

    for i_randsamp = 1:50
    % m_number = 1;
    fea_number_con = [];
    fea_number_in = [];
    
        for i_sub = 1:length(subjects)
            m_number = 1;

            % inputPath = fullfile('outputDataChronux_zscore',subject);
            inputPath1 = fullfile('c_outputPowerData_log','highGamma',subjects{i_sub});
            % inputPath2 = fullfile('outputPowerData_nolog','theta',subjects{i_sub});

            [fea_con_tmp, fea_in_tmp, m_number_out,sel_chan_number, n] = concatenateFeatures(m_number,config,inputPath1);
            fea_number_con = [fea_number_con, fea_con_tmp]; % Concatenate horizontally
            fea_number_in = [fea_number_in, fea_in_tmp];
            
            if exist('inputPath2', 'var') % do this to do theta and gamma together
                [fea_con_tmp, fea_in_tmp, m_number_out,sel_chan_number, n] = concatenateFeatures(m_number,config,inputPath2);
                fea_number_con = [fea_number_con, fea_con_tmp]; % Concatenate horizontally
                fea_number_in = [fea_number_in, fea_in_tmp];
            end
        end
    
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% SVM %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
        n_sample = n;
   
        if m_number_out ~= 0

            X = real([fea_number_con; fea_number_in]);   % samples x features
            Y = [zeros(n_sample,1); ones(n_sample,1)];  % labels
            cv = cvpartition(Y,'KFold',10);

            for i = 1:cv.NumTestSets
                train_idx = training(cv,i);
                test_idx  = test(cv,i);
                X_train = X(train_idx,:);
                Y_train = Y(train_idx);
                X_test  = X(test_idx,:);
                Y_test  = Y(test_idx);
                % Train SVM
                Mdl = fitcsvm(X_train, Y_train, 'Standardize', true, 'KernelFunction', 'linear');
            
                % Beta Coefficients
                beta = Mdl.Beta;
                % abs_beta = abs(beta); % abs value taken later
            
                nChannels = length(beta);
                means_idx = 1:2:nChannels;
                max_idx   = 2:2:nChannels;
            
                if i_randsamp == 1 && i == 1
                    mean_beta = zeros(length(means_idx), 10, 50);
                    max_beta  = zeros(length(max_idx), 10, 50);
                    allShapVals = cell(10,50);
                    meanShapVals = zeros(nChannels, 10, 50);
                end
            
                mean_beta(:,i,i_randsamp) = beta(means_idx);
                max_beta(:,i,i_randsamp)  = beta(max_idx);
            
                % Prediction accuracy
                labels = predict(Mdl, X_test);
                correct_number(i_randsamp,i) = mean(labels == Y_test)*100;
            
                % ~~~~ Shapley Values ~~~~~
                S = shapley(Mdl, X_train);  % X_train is  background
                S = fit(S, X_test);         % X_test is query points
            
                fixedClass = 1;
                nsamples = size(X_test,1);
                npredictors = size(X_test,2);
            
                shapVals = zeros(npredictors, nsamples);
                for il = 1:nsamples
                    for f = 1:npredictors
                        shapVals(f,il) = S.Shapley{f, num2str(fixedClass)}(il);
                    end
                end
            
                allShapVals{i,i_randsamp} = abs(shapVals);
                meanShapVals(:,i,i_randsamp) = mean(abs(shapVals),2);
            
                clear Mdl
            end

        else
            correct_number = 0;
            mean_beta = 0;
            max_beta = 0;

        end
       
    end 

    % meanShapVals = squeeze(mean(meanShapVals, [2 3]));

    y = [mean(correct_number(:))];
    % calculate SE with n=50 random samplings --> divide by /sqrt(50)
    err = std(correct_number(:))/sqrt(length(correct_number));
    shap = meanShapVals;
    % shap = mean(meanShapVals(:));
    mean_weights = mean_beta;
    max_weights = max_beta;
end




function barplot(y, xlabels, err,config)

    [nBars, nGroups] = size(y);

    figure;
    if size(y,1)==1 % if b only has one grouping
        xPos = 1:1.5:nGroups*1.5;  % space bars 1.5 units apart
        for i = 1:nGroups
            % b(i) = bar(i, y(i), 'BarWidth', 0.5);
            b(i) = bar(xPos(i), y(i), 'BarWidth', 0.5);
            hold on
        end
        % x = 1:nGroups;  % positions match bar indices
        % xlim([0 (nBars+2)]);
        xlim([0 xPos(end) + 1]);
        xticks(xPos);
        for i = 1:nGroups
            b(i).Labels = b(i).YData;
        end
        line([0 (nGroups+1)],[50 50],'color','k','linestyle','--','linewidth',1);
        xticks(1:nGroups);xlabel('Feature Extraction Method');
        ylim([00 125]); 
    else
        b = bar(y, 'BarWidth', 1);
        xlim([0 (nBars+1)]);
        xticks(1:nBars); xtickangle(45); xticklabels(xlabels);xlabel('Patient');
        x = nan(nBars, nGroups);
        for i = 1:nGroups
          x(:,i) = b(i).XEndPoints;
          % b(i).Labels = b(i).YData;
        end
        line([0 (nBars+1)],[50 50],'color','k','linestyle','--','linewidth',1);
        ylim([00 100]);
    end

    for iGroup = 1:nGroups
        color = [0.25 + 0.75 * (iGroup / nGroups), 0, 0];  % Dark red to light red gradient
        b(iGroup).FaceColor = color;
    end
    
    hold on;

    err=err.*10; % fix error from /500 to /50
    if exist('err', 'var')
        er = errorbar(x,y,err); 
        for ier = 1:numel(er)
        er(ier).LineStyle = 'none'; er(ier).CapSize = 5; er(ier).Color = [0 0 0];
        end
    end

    ylabel('Accuracy (%)');
    title('10 fold C-V, 50 sessions','FontSize',16);   
    set(gca,'box','off','tickDir','out')
    grid on;
    
    if exist('config', 'var')
        legend(b,config,'Location', 'northwest')
    end
end

function weightsplot(mean_weights,max_weights,sel_chan_number,subject)

    RegionLabels = {'dlPFC', 'dmPFC', 'OFC', 'vlPFC', 'STG', 'MTG', 'ITG', 'dACC', 'AMY', 'HIP'};
    % inputPath = fullfile('outputDataChronux_ratio', subject);
    inputPath = fullfile('outputPowerData_nolog','theta',subject);
    filesToLoad = {'ROIbyChannel.mat'};
    
    for i = 1:length(filesToLoad)
        filePath = char(fullfile(inputPath, filesToLoad{i}));  % Convert to char
        ROIchannels = cell(1,length(filesToLoad));
        for i_sub = 1:length(subject)
            load((filePath(i_sub, :)));
            ROIchannels{i_sub} = ROIbyChannel;
        end
    end

    nPatients = length(ROIchannels); % Number of patients
    nROIs = numel(ROIchannels{1});  % Number of ROIs, inferred from one patient's data (e.g., first cell)
    channels_per_ROI = zeros(nROIs,nPatients);  % Matrix to store the number of channels for each ROI
    for i = 1:nPatients
        for j = 1:nROIs
            channels_per_ROI(j,i) = numel(ROIchannels{i}{j});  % Store the number of channels for each ROI
        end
    end
    figure;
    b0 = bar(channels_per_ROI, 'grouped');  % Transpose so that each ROI is a group, and each bar in a group is a patient
    xlabel('ROI'); ylabel('Number of Channels');
    xticks(1:nROIs); xticklabels(RegionLabels);
    title('Number of Channels per ROI for Each Patient');
    legend(b0, subject, 'Location', 'Best');
    for ig = 1:size(channels_per_ROI,2)
        color = [0.25 + 0.75 * (ig / size(channels_per_ROI,2)), 0, 0];  % Dark red to light red gradient
        b0(ig).FaceColor = color;
    end

    mean_abs_beta_flat = reshape(mean_weights, size(mean_weights,1), []);  % now size = [60 x 500]
    max_abs_beta_flat = reshape(max_weights, size(max_weights,1), []);

    beta_mean = squeeze(mean(mean_weights,[2 3]));
    beta_max = squeeze(mean(max_weights,[2 3]));

    beta_mean_err = std(mean_abs_beta_flat,0,2)/sqrt(size(mean_abs_beta_flat, 2)); 
    beta_max_err = std(max_abs_beta_flat,0,2)/sqrt(size(max_abs_beta_flat, 2)); 


    if exist("ROIbyChannel",'var')
        beta_mean_plot = zeros(10,length(subject));
        beta_max_plot = zeros(10,length(subject));
        beta_mean_err_plot = zeros(10,length(subject));
        beta_max_err_plot = zeros(10,length(subject));
            
        for i_sub = 1:length(subject)        
            nROI = numel(ROIchannels{i_sub});
        
            for iROI = 1:nROI
                
                subjectROI = ROIchannels{i_sub};
                channelsinROI = find(ismember(sel_chan_number, subjectROI{iROI}));

                if isempty(channelsinROI)
                    beta_mean_plot (iROI,i_sub) = 0;
                    beta_max_plot (iROI,i_sub) = 0;
                    beta_mean_err_plot (iROI,i_sub) = 0;
                    beta_max_err_plot (iROI,i_sub) = 0;
                else
                    beta_mean_plot (iROI,i_sub) = mean(beta_mean(channelsinROI));
                    beta_max_plot (iROI,i_sub) = mean(beta_max(channelsinROI));
                    beta_mean_err_plot (iROI,i_sub) = mean(beta_mean_err(channelsinROI));
                    beta_max_err_plot (iROI,i_sub) = mean(beta_max_err(channelsinROI));
                end
            end
        end
    end
    
    beta_mean_plot = mean(beta_mean_plot,2);
    beta_mean_err_plot = mean(beta_mean_err_plot,2);
    beta_max_plot = mean(beta_max_plot,2);
    beta_max_err_plot = mean(beta_max_err_plot,2);


    [nBars, nGroups] = size(beta_mean_plot);

    figure;
    tiledlayout(1, 2, 'TileSpacing', 'compact');
    
    meanplot = nexttile;
    b1 = bar(beta_mean_plot, 'grouped','BarWidth', 1);
    xlim([0 (nBars+1)]); xticks(1:nBars); xticklabels(RegionLabels);

    x = nan(nBars, nGroups);
    for i = 1:nGroups
        x(:,i) = b1(i).XEndPoints;
    end
        
    hold on;
    er = errorbar(x, beta_mean_plot, beta_mean_err_plot); 
    for ier = 1:numel(er)
        er(ier).LineStyle = 'none'; er(ier).CapSize = 5; er(ier).Color = [0 0 0];
    end
   
    for ig = 1:nGroups
        color = [0.25 + 0.75 * (ig / nGroups), 0, 0];  % Dark red to light red gradient
        b1(ig).FaceColor = color;
    end

    set(gca, 'FontSize', 10);
    xlabel('Region of Interest', 'FontSize', 12);
    ylabel('Absolute Beta Weight', 'FontSize', 12);
    title('Mean Beta Weights', 'FontSize', 14);
    grid on;

    maxplot = nexttile;
    b2 = bar(beta_max_plot, 'grouped','BarWidth', 1);
    xlim([0 (nBars+1)]); xticks(1:nBars); xticklabels(RegionLabels);
    hold on;
    er = errorbar(x, beta_max_plot, beta_max_err_plot); 
    for ier = 1:numel(er)
        er(ier).LineStyle = 'none'; er(ier).CapSize = 5; er(ier).Color = [0 0 0];
    end

    for ig = 1:nGroups
        color = [0.25 + 0.75 * (ig / nGroups), 0, 0];  % Dark red to light red gradient
        b2(ig).FaceColor = color;
    end

    xlabel('Region of Interest', 'FontSize', 12);
    title('Max Beta Weights', 'FontSize', 14);
    grid on;
    
    if nGroups ~= 1
        legend(b1, subject, 'Location', 'eastoutside');
    else
    end

    linkaxes([meanplot, maxplot], 'y');
    sgtitle('Feature Contributions Across Regions', 'FontSize', 16);
end
