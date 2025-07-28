%addpath_recurse('functions_general');
load(fullfile('features.mat'));
tic

n = 34; 
m_number = 1;
subjects = {'BW42', 'MG51b'};
fea_number_con = [];
fea_number_in = [];

for i_randsamp = 1:50

    for i_sub = 1:length(subjects)
        [fea_con_tmp, fea_in_tmp, m_number] = concatenateFeatures(features,subjects{i_sub}, m_number, n);
        fea_number_con = [fea_number_con, fea_con_tmp]; % Concatenate horizontally
        fea_number_in = [fea_number_in, fea_in_tmp];
    end

    %% SVM

        fea_number_con2 = fea_number_con;
        fea_number_in2 = fea_number_in;
        n_sample = n;

    for i_iter = 1
        % Number
        [train_ind, test_ind,n_test] = generateCrossValInd(n_sample); % n_sample = 52;
        for i = 1:10 % 10-fold 
            X_train = [fea_number_con2(train_ind(i,:),:);fea_number_in2(train_ind(i,:),:)];
            Y_train = [zeros(n_sample-n_test,1);ones(n_sample-n_test,1)];
            Mdl = fitcsvm(X_train,Y_train,'Standardize',true,'KernelFunction','linear');

            X_test = [fea_number_con2(test_ind(i,:),:);fea_number_in2(test_ind(i,:),:)];
            labels = predict(Mdl,X_test);
            Y_test = [zeros(n_test,1);ones(n_test,1)];
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

correct_number = correct_number(:);

%% Plot

y = [mean(correct_number(:))];

err = std(correct_number(:))/sqrt(numel(correct_number)); 

x = 1;

b = bar(x,y,'FaceColor',[0.511 0.515 1],'BarWidth', 0.4);

hold on;
er = errorbar(x,y,err,err);    
er.Color = [0 0 0];                            
er.LineStyle = 'none'; 
er.CapSize = 5;

for i = 1:length(x)
    text(x(i), y(i) + err(i) + 1, sprintf('%.1f', y(i)), ...
        'HorizontalAlignment', 'center','VerticalAlignment', 'bottom', ...
        'FontSize', 10);
end


ylim([00 100])
xlim([0 2])
line([0 6],[50 50],'color','k','linestyle','--','linewidth',1.5)
set(gca,'XTickLabel',[]);set(gca,'XTick',[]);
ylabel('Accuracy (%)');
xlabel('MGH 2 subjects');
title('10 fold cross validation SVM with 50 sessions');
set(gca,'fontsize', 10,'box','off','FontName','Arial','tickDir','out')
toc

%% functions

function [fea_number_con, fea_number_in, m_number_out] = concatenateFeatures(features, subject, m_number, n)
    fea_number_con = [];
    fea_number_in = [];

    sel_chan_number = features.(['selectedChannels_' subject]);

    if ~isempty(sel_chan_number)
        conPower = features.(['conPowerFeatures_' subject]);
        inPower = features.(['inPowerFeatures_' subject]);

        for i = 1:length(sel_chan_number)
            ch = sel_chan_number(i);

            % --- Max Power: Column 2 ---
            tent1 = cellfun(@(x) x(ch,2), conPower);  % con max across trials
            tent2 = cellfun(@(x) x(ch,2), inPower);   % in max across trials
            fea_number_con(1:n, m_number) = randsample(tent1, n);
            fea_number_in(1:n, m_number) = randsample(tent2, n);
            m_number = m_number + 1;

            % --- Mean Power: Column 1 ---
            tent1 = cellfun(@(x) x(ch,1), conPower);  % con mean across trials
            tent2 = cellfun(@(x) x(ch,1), inPower);   % in mean across trials
            fea_number_con(1:n, m_number) = randsample(tent1, n);
            fea_number_in(1:n, m_number) = randsample(tent2, n);
            m_number = m_number + 1;
        end

        m_number_out = m_number;
    end
end
