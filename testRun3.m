
% Assume your data is all in a subdirectory called "patientData"
% Must have functions added to path: generateCrossValInd.m , permutationTest.m

% Files:
% 'BW42', 'MG51b', 'MG79', 'MG86', ...
% 'MG89', 'MG90', 'MG91', 'MG95', ...
% 'MG96', 'MG99', 'MG102', 'MG104', ...
% 'MG105', 'MG106', 'MG111', 'MG112',...
% 'MG116', 'MG117', 'MG118', 'MG120'

tic

n = 40; % # correct trials
subjects = {'BW42', 'MG51b', 'MG79', 'MG86', ...
'MG89', 'MG90', 'MG95', ...
'MG96', 'MG99', 'MG102', 'MG104', ...
'MG105', 'MG106', 'MG111', 'MG112',...
'MG116', 'MG117', 'MG118', 'MG120'};

config = {'selChans - z-scored power','confChans - z-scored power', ...
    'selChans - non z-scored power','confChans - non z-scored power'};

% nBars = length(subjects);
nBars = 1;
nGroups = length(config);

groupedBars = zeros(nBars,nGroups);
groupedErr = zeros(nBars,nGroups);

for igroup = 1:nGroups
    if nBars ==1
        [y, err] = SVM(subjects(:),n,config(igroup));
        groupedBars(:,igroup) = y;
        groupedErr(:,igroup) = err;
    else
        for ibar=1:nBars
            [y, err] = SVM(subjects(ibar),n,config(igroup));
            groupedBars(ibar,igroup) = y;
            groupedErr(ibar,igroup) = err;
        end
    end
end

barplot(groupedBars, groupedErr, subjects, config)
toc
%% functions

function [fea_number_con, fea_number_in, m_number_out] = concatenateFeatures(subject, m_number, n, config)
    
    inputPath = fullfile('outputData', subject);
    
    filesToLoad = {'selectedChans.mat','conflictModChans.mat', ...
        'conPowerFeatures.mat','inPowerFeatures.mat', ...
        'conBandPowerFeatures.mat', 'inBandPowerFeatures.mat'};
    
    for i = 1:length(filesToLoad)
        load(fullfile(inputPath, filesToLoad{i}));
    end

     switch config{1}
         case 'confChans - z-scored power'
            sel_chan_number = conflictModChans;
            conPower = conPowerFeatures;
            inPower = inPowerFeatures;
         case 'selChans - z-scored power'
            sel_chan_number = selectedChans;
            conPower = conPowerFeatures;
            inPower = inPowerFeatures;
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
    end
end

function [y, err] = SVM(subjects,n,config)

    for i_randsamp = 1:50
    m_number = 1;
    fea_number_con = [];
    fea_number_in = [];
    
        for i_sub = 1:length(subjects)
            [fea_con_tmp, fea_in_tmp, m_number] = concatenateFeatures(subjects{i_sub}, m_number, n,config);
            fea_number_con = [fea_number_con, fea_con_tmp]; % Concatenate horizontally
            fea_number_in = [fea_number_in, fea_in_tmp];
        end
    
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% SVM %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
        n_sample = n;
    
        for i_iter = 1
            % Number
            [train_ind, test_ind,n_test] = generateCrossValInd(n_sample); % n_sample = 52;
            for i = 1:10 % 10-fold 
                X_train = [fea_number_con(train_ind(i,:),:);fea_number_in(train_ind(i,:),:)]; % made real
                Y_train = [zeros(n_sample-n_test,1);ones(n_sample-n_test,1)];
                Mdl = fitcsvm(real(X_train),Y_train,'Standardize',true,'KernelFunction','linear');
    
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
        
    end 

    y = [mean(correct_number(:))];
    err = std(correct_number(:))/sqrt(numel(correct_number)); 
end

function barplot(y, err, xlabels,config)

    [nBars, nGroups] = size(y);

    if size(y,1)==1 % if b only has one grouping
        for i = 1:nGroups
            b(i) = bar(i, y(i), 'BarWidth', 1);
            hold on
        end
        x = 1:nGroups;  % positions match bar indices
        for i = 1:nGroups
            b(i).Labels = b(i).YData;
        end
        line([0 (nGroups+1)],[50 50],'color','k','linestyle','--','linewidth',1);
        xticks(1:nGroups);xlabel('Feature Extraction Method');
        ylim([00 125]); 
    else
        b = bar(y, 'BarWidth', 1);
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

    for i = 1:4
        if i <= 2
            b(i).FaceColor = shades(i,:,1); % first group blue shades
        else
            b(i).FaceColor = shades(i-2,:,2); % second group red shades
        end
    end
    
    hold on;

    er = errorbar(x,y,err); 

    for ier = 1:numel(er)
        er(ier).LineStyle = 'none'; er(ier).CapSize = 5; er(ier).Color = [0 0 0];
    end

    ylabel('Accuracy (%)');
    title('Group 2 - 10 fold C-V, 50 sessions','FontSize',16);   
    set(gca,'box','off','tickDir','out')
    legend(b,config,'Location', 'northwest')
end

