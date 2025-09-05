
% Assume your data is all in a subdirectory called "patientData"
% Must have functions added to path: (1) generateCrossValInd.m , (2)
% permutationTest.m

tic

n = 40;   % At least 34 correct trials
subjects = {'BW42'}; % BW42, MG51b

% Files:
% 'BW42', 'MG51b', 'MG79', 'MG86', 
% 'MG89', 'MG90', 'MG91', 'MG95', 
% 'MG96', 'MG99', 'MG102', 'MG104', 
% 'MG105', 'MG106', 'MG111', 'MG112',
% 'MG116', 'MG117', 'MG118', 'MG120'

for i_randsamp = 1:50
m_number = 1;
fea_number_con = [];
fea_number_in = [];

    for i_sub = 1:length(subjects)
        [fea_con_tmp, fea_in_tmp, m_number] = concatenateFeatures(subjects{i_sub}, m_number, n);
        fea_number_con = [fea_number_con, fea_con_tmp]; % Concatenate horizontally
        fea_number_in = [fea_number_in, fea_in_tmp];
    end

    %% SVM

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

%% Plot

y = [mean(correct_number(:))];
err = std(correct_number(:))/sqrt(numel(correct_number)); 

x = 1;

b = bar(x,y,'FaceColor',[0.511 0.515 1],'BarWidth', 0.4);

hold on;
er = errorbar(x,y,err,err); er.Color = [0 0 0]; er.LineStyle = 'none'; er.CapSize = 5;
for i = 1:length(x)
    text(x(i), y(i) + err(i) + 1, sprintf('%.1f', y(i)),'HorizontalAlignment', ...
        'center','VerticalAlignment', 'bottom','FontSize', 10);
end

xticks([1 2 3 4 5]); xlim([0 3]);
xticklabels({'Band Power','Z Scores'});xlabel('Features Extraction Location');

ylim([00 100]); line([0 7],[50 50],'color','k','linestyle','--','linewidth',1.5)
ylabel('Accuracy (%)');

title('Group 2 - 10 fold C-V, 50 sessions','FontSize',16);

set(gca,'fontsize', 10,'box','off','FontName','Arial','tickDir','out')
toc

%% functions

function [fea_number_con, fea_number_in, m_number_out] = concatenateFeatures(subject, m_number, n)
    
    inputPath = fullfile('outputData', subject);
    
    filesToLoad = {'selectedChans.mat','conflictModChans.mat', ...
        'conPowerFeatures.mat','inPowerFeatures.mat', ...
        'conBandPowerFeatures.mat', 'inBandPowerFeatures.mat'};
    
    for i = 1:length(filesToLoad)
        load(fullfile(inputPath, filesToLoad{i}));
    end

    sel_chan_number = selectedChans;
    conPower = conPowerFeatures;
    inPower = inPowerFeatures;

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