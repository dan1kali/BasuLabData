
% Assume your data is all in a subdirectory called "patientData"
% Must have functions added to path: generateCrossValInd.m , permutationTest.m

% Files:
% 'BW42', 'MG51b', 'MG79', 'MG86', ...
% 'MG89', 'MG90', 'MG91', 'MG95', ...
% 'MG96', 'MG99', 'MG102', 'MG104', ...
% 'MG105', 'MG106', 'MG111', 'MG112',...
% 'MG116', 'MG117', 'MG118', 'MG120'


% 'UCMC02', 'UCMC03', 'UCMC05', ...
% 'UCMC06', 'UCMC07', 'UCMC08', 'UCMC09', ...
% 'UCMC11', 'UCMC13', 'UCMC14'



% 'UCMC01', 'UCMC02', 'UCMC03', 'UCMC05', ...
% 'UCMC06', 'UCMC07', 'UCMC08', 'UCMC09', ...
% 'UCMC11', 'UCMC13', 'UCMC14', 'UCMC15', 'UCMC17'


% 'selChans - z-scored power','confChans - z-scored power', ...
%   'selChans - normalized log power',  'confChans - normalized log power'
%     'selChans - non z-scored power','confChans - non z-scored power' ...
%      'selChans - z-scored power - rand','confChans - z-scored power - rand', ...
%      'selChans - z-scored power - res','confChans - z-scored power - res', ...

%% Plot Patient Data
tic

n = 40; % minumum # of trials
subjects = {'UCMC01', 'UCMC02', 'UCMC03', 'UCMC05', ...
'UCMC06', 'UCMC07', 'UCMC08', 'UCMC09', ...
'UCMC11', 'UCMC13', 'UCMC14', 'UCMC15', 'UCMC17'};

config = {'confChans - z-scored power'};

nBars = length(subjects);
% nBars = 1;
nGroups = length(config);

groupedBars = zeros(nBars,nGroups);
groupedErr = zeros(nBars,nGroups);

for igroup = 1:nGroups
    % try
    if nBars ==1
        [y, err,mean_weights,max_weights,sel_chan_number]  = SVM(subjects(:),n,config(igroup));
        groupedBars(:,igroup) = y;
        groupedErr(:,igroup) = err;
    else
        for ibar=1:nBars
            [y, err,mean_weights,max_weights,sel_chan_number]  = SVM(subjects(ibar),n,config(igroup));
            groupedBars(ibar,igroup) = y;
            groupedErr(ibar,igroup) = err;
        end
    end
    % catch ME
    % disp(getReport(ME, 'extended'));
    % end
end

% barplot(groupedBars, subjects, groupedErr, config)
weightsplot(mean_weights,max_weights,sel_chan_number,subjects) % only most recent sub and condition
toc

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

function [fea_number_con, fea_number_in, m_number_out,sel_chan_number] = concatenateFeatures(subject, m_number, n, config)
    
    inputPath = fullfile('outputDataChronuxUC', subject);
    
    filesToLoad = {'selectedChans.mat','conflictModChans.mat', ...
        'conBandPowerFeatures.mat', 'inBandPowerFeatures.mat', ...
        'conPowerFeatures.mat','inPowerFeatures.mat', ...
        'conResPowerFeatures.mat','inResPowerFeatures.mat'};        

        % 'allPowerFeatures.mat', ...
        % 'ROIbyChannel.mat',...
        % 'conLogPowerFeatures.mat','inLogPowerFeatures.mat', ...
       
    for i = 1:length(filesToLoad)
        load(fullfile(inputPath, filesToLoad{i}));
    end


    % nTrials = numel(allPowerFeatures); 
    % randomOrder = randperm(nTrials);
    % half = floor(nTrials/2);
    % conPowerFeaturesRand = allPowerFeatures(randomOrder(1:half));
    % inPowerFeaturesRand = allPowerFeatures(randomOrder(half+1:end));


     switch config{1}
         case 'confChans - z-scored power'
            sel_chan_number = conflictModChans;
            conPower = conPowerFeatures;
            inPower = inPowerFeatures;
         case 'selChans - z-scored power'
            sel_chan_number = selectedChans;
            conPower = conPowerFeatures;
            inPower = inPowerFeatures;

         % Log power only (non response aligned)
         case 'confChans - normalized log power'
            sel_chan_number = conflictModChans;
            conPower = conPowerFeatures;
            inPower = inPowerFeatures;
         case 'selChans - normalized log power'
            sel_chan_number = selectedChans;
            conPower = conPowerFeatures;
            inPower = inPowerFeatures;

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
         

         otherwise
            error('Unknown config: %s\n', config{1});
    end

    % sel_chan_number = selectedChans;
    % conPower = conBandPowerFeatures;
    % inPower = inBandPowerFeatures;

    fea_number_con = [];
    fea_number_in = [];

    if ~isempty(sel_chan_number)
        for i = 1:length(sel_chan_number)
            ch = sel_chan_number(i);

            % --- Max Power: Pull from Column 2 ---
            tent1 = cellfun(@(x) x(ch,2), conPower);  % con max across trials
            tent2 = cellfun(@(x) x(ch,2), inPower);   % in max across trials
            fea_number_con(1:n, m_number) = randsample(tent1, n);
            fea_number_in(1:n, m_number) = randsample(tent2, n);
            m_number = m_number + 1;

            % --- Mean Power: Pull from Column 1 ---
            tent1 = cellfun(@(x) x(ch,1), conPower);  % con mean across trials
            tent2 = cellfun(@(x) x(ch,1), inPower);   % in mean across trials
            fea_number_con(1:n, m_number) = randsample(tent1, n);
            fea_number_in(1:n, m_number) = randsample(tent2, n);
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

function [y, err,mean_weights,max_weights,sel_chan_number] = SVM(subjects,n,config)

    for i_randsamp = 1:50
    % m_number = 1;
    fea_number_con = [];
    fea_number_in = [];
    
        for i_sub = 1:length(subjects)
            m_number = 1;
            [fea_con_tmp, fea_in_tmp, m_number_out,sel_chan_number] = concatenateFeatures(subjects{i_sub}, m_number, n,config);
            fea_number_con = [fea_number_con, fea_con_tmp]; % Concatenate horizontally
            fea_number_in = [fea_number_in, fea_in_tmp];
        end
    
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% SVM %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
        n_sample = n;
    


        if m_number_out ~= 0

            for i_iter = 1
                % Number
                [train_ind, test_ind,n_test] = generateCrossValInd(n_sample); % n_sample = 52;
                for i = 1:10 % 10-fold 
                    X_train = [fea_number_con(train_ind(i,:),:);fea_number_in(train_ind(i,:),:)]; % made real
                    Y_train = [zeros(n_sample-n_test,1);ones(n_sample-n_test,1)];
                    Mdl = fitcsvm(real(X_train),Y_train,'Standardize',true,'KernelFunction','linear');
        
                    beta = Mdl.Beta;
                    abs_beta = abs(beta);
    
                    nChannels = floor(length(abs_beta));
                    means_idx = 1:2:nChannels;
                    max_idx = 2:2:nChannels;
    
                    if i_randsamp == 1 && i == 1
                        mean_abs_beta = zeros(length(means_idx), 10, 50);
                        max_abs_beta  = zeros(length(max_idx), 10, 50);
                    end
    
                    % mean_abs_beta (:,i,i_randsamp) = abs_beta(means_idx);
                    % max_abs_beta (:,i,i_randsamp) = abs_beta(max_idx);
                    mean_abs_beta (:,i,i_randsamp) = beta(means_idx);
                    max_abs_beta (:,i,i_randsamp) = beta(max_idx);
    
                    X_test = [fea_number_con(test_ind(i,:),:);fea_number_in(test_ind(i,:),:)];
                    labels = predict(Mdl,real(X_test)); % made real
                    Y_test = [zeros(n_test,1);ones(n_test,1)]; % ground truth
                    n_correct = 0;
                    for j = 1:length(labels)
                        if labels(j)==Y_test(j)
                            n_correct = n_correct+1;
                        end
                    end
                    correct_number(i_randsamp,i) = n_correct/length(Y_test)*100;
                    clear Mdl
                end
            end
        else
            correct_number = 0;
            mean_abs_beta = 0;
            max_abs_beta = 0;

        end
       
    end 

    y = [mean(correct_number(:))];
    err = std(correct_number(:))/sqrt(numel(correct_number)); 
    mean_weights = mean_abs_beta;
    max_weights = max_abs_beta;
end

function barplot(y, xlabels, err,config)

    [nBars, nGroups] = size(y);

    if size(y,1)==1 % if b only has one grouping
        for i = 1:nGroups
            b(i) = bar(i, y(i), 'BarWidth', 1);
            hold on
        end
        x = 1:nGroups;  % positions match bar indices
        xlim([0 (nBars+1)]);
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


    % Bar Shading
    baseColors = [0.1   0.4   0.7; 0.8   0.3   0.25]; % set shades, [blue;red]
    targetColors = [0.5 0.7 1;1   0.5 0.3   ];
    nShades = 2;
    % (:,:,1) for blues, (:,:,2) for reds - (shades x RGB x groups)
    shades = zeros(nShades, 3, size(baseColors,1));
    for g = 1:size(baseColors,1)
        for c = 1:3  % R, G, B channels
            shades(:, c, g) = linspace(baseColors(g, c), targetColors(g, c), nShades);
        end
    end

    if nGroups ==4
        for i = 1:nGroups
            if i <= 2
                b(i).FaceColor = shades(i,:,1); % first group blue shades
            else
                b(i).FaceColor = shades(i-2,:,2); % second group red shades
            end
        end
    end
    
    hold on;

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

    inputPath = fullfile('outputDataChronuxUC', subject);
    filesToLoad = {'ROIbyChannel.mat'};
    
    for i = 1:length(filesToLoad)
        filePath = char(fullfile(inputPath, filesToLoad{i}));  % Convert to char
        ROIchannels = cell(1,length(filesToLoad));
        for i_sub = 1:length(subject)
            load((filePath(i_sub, :)));
            ROIchannels{i_sub} = ROIbyChannel;
        end
    end

    % for i_sub = 1:length(subjects)
    %     inputPath = fullfile('outputDataWaveletUC', subject);
    %     filesToLoad = {'ROIbyChannel.mat'};
    %     for i = 1:length(filesToLoad)
    %         load(fullfile(inputPath{i_sub}, filesToLoad{i}));
    %     end
    % end

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

    RegionLabels = {'dlPFC', 'dmPFC', 'OFC', 'vlPFC', 'STG', 'MTG', 'ITG', 'dACC', 'AMY', 'HIP'};

    
    beta_mean_plot = mean(beta_mean_plot,2);
    beta_mean_err_plot = mean(beta_mean_err_plot,2);
    beta_max_plot = mean(beta_max_plot,2);
    beta_max_err_plot = mean(beta_max_err_plot,2);


    [nBars, nGroups] = size(beta_mean_plot);

    figure;
    tiledlayout(1, 2);
    
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
    set(gca, 'FontSize', 10);
    xlabel('Region of Interest', 'FontSize', 12);
    ylabel('Absolute Beta Weight', 'FontSize', 12);
    title('Feature Importance - Mean', 'FontSize', 14);
    grid on;
    
    maxplot = nexttile;
    b2 = bar(beta_max_plot, 'grouped','BarWidth', 1);
    xlim([0 (nBars+1)]); xticks(1:nBars); xticklabels(RegionLabels);
    hold on;
    er = errorbar(x, beta_max_plot, beta_max_err_plot); 
    for ier = 1:numel(er)
        er(ier).LineStyle = 'none'; er(ier).CapSize = 5; er(ier).Color = [0 0 0];
    end
    xlabel('Region of Interest', 'FontSize', 12);
    title('Feature Importance - Max', 'FontSize', 14);
    legend(b2, subject, 'Location', 'northwest');
    grid on;
    
    if nGroups ~= 1
        legend(b1, subject, 'Location', 'northwest');
        legend(b2, subject, 'Location', 'northwest');
    else
    end

    linkaxes([meanplot, maxplot], 'y');
end
