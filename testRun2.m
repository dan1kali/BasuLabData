%addpath_recurse('functions_general');
load(fullfile('features.mat'));
tic

n = 59;   % At least 34 correct trials
subjects = {'BW42','MG51b'}; % BW42, MG51b

for i_randsamp = 1:50
m_number = 1;
fea_number_con = [];
fea_number_in = [];

    for i_sub = 1:length(subjects)
        [fea_con_tmp, fea_in_tmp, m_number] = concatenateFeatures(features,subjects{i_sub}, m_number, n);
        fea_number_con = [fea_number_con, fea_con_tmp]; % Concatenate horizontally
        fea_number_in = [fea_number_in, fea_in_tmp];
    end

    %% SVM

    n_sample = n;

    for i_iter = 1
        % Number
        [train_ind, test_ind,n_test] = generateCrossValInd(n_sample); % n_sample = 52;
        for i = 1:10 % 10-fold 
            X_train = [fea_number_con(train_ind(i,:),:);fea_number_in(train_ind(i,:),:)];
            Y_train = [zeros(n_sample-n_test,1);ones(n_sample-n_test,1)];
            Mdl = fitcsvm(X_train,Y_train,'Standardize',true,'KernelFunction','linear');

            X_test = [fea_number_con(test_ind(i,:),:);fea_number_in(test_ind(i,:),:)];
            labels = predict(Mdl,X_test);
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

xticks([1 2 3]); xticklabels({'Original','Permutation Test','Visual'});xlabel('Channel Rejection Method');

ylim([00 100])
xlim([0 4])
line([0 6],[50 50],'color','k','linestyle','--','linewidth',1.5)
ylabel('Accuracy (%)');
title('10 fold cross validation SVM with 50 sessions');
set(gca,'fontsize', 10,'box','off','FontName','Arial','tickDir','out')
toc

%% functions

function [fea_number_con, fea_number_in, m_number_out] = concatenateFeatures(features, subject, m_number, n)
    fea_number_con = [];
    fea_number_in = [];

    sel_chan_number = features.(['selectedChan_' subject]);
    % sel_chan_number = features.(['selectedChan_' subject '_confMod_a10' ]);
    % sel_chan_number = features.selectedChan_BW42_confMod_vis;

    if ~isempty(sel_chan_number)
        conPower = features.(['conPowerFeatures_' subject]);
        inPower = features.(['inPowerFeatures_' subject]);

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
