
% download data here: https://drive.google.com/drive/folders/1RSbYruG8tFLI_rQjXn1X09Lo4jbq-ZfT?usp=sharing

%% high gamma
% Load the region accuracies data
load('regionAccuracies_nolog_hg.mat');    % select for high gamma

% run the plot of interest
allRegionPlot
groupRegionPlot

% Load the subject accuracies data
load('subjectAccuracies_nolog_hg.mat');    % select for high gamma
barplot

%% theta
% Load the region accuracies data
load('regionAccuracies_nolog_theta.mat')  % select for theta

% run the plot of interest
allRegionPlot
groupRegionPlot

% Load the subject accuracies data
load('subjectAccuracies_nolog_theta.mat')  % select for theta
barplot



%% functions
function allRegionPlot

nPatients = numel(regionAccuracies);
nRegions = size(regionAccuracies{1}, 1);
RegionLabels = {'dlPFC', 'dmPFC', 'OFC', 'vlPFC', 'STG', 'MTG', 'ITG', 'dACC', 'AMY', 'HIP'};

allAcc = nan(nPatients, nRegions);
allErr = nan(nPatients, nRegions);
for i = 1:nPatients
    temp = regionAccuracies{i};
    temp(temp == 0) = NaN; % ignore zeros
    allAcc(i, :) = temp(:, 1)';
    allErr(i, :) = temp(:, 2)';
end

figure;

%%%%%% Choose Plot Type by uncommmenting: %%%%%%%%

% Bar
% regionAcc = mean(allAcc, 'omitnan');
% regionErr = std(allAcc, 'omitnan') ./ sqrt(sum(~isnan(allAcc)));
% bar(regionAcc, 'FaceColor', [0.2 0.6 0.8]); hold on;
% errorbar(1:nRegions, regionAcc, regionErr, 'k.', 'LineWidth', 1.2);

% % Boxplot
% boxplot(allAcc, 'Labels', RegionLabels); hold on;

% % Scatter Plot
rng(1); % for reproducibility
jitterAmount = 0.1;
x=repmat(1:length(RegionLabels),size(allAcc,1),1);
for j=1:size(x,2)
    x(:,j) = x(:,j) + jitterAmount * randn(size(x(:,j)));
    hold on;
    scatter(x(:,j), allAcc(:,j), 25, 'r', 'filled');
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

set(gca, 'XTick', 1:nRegions, 'XTickLabel', RegionLabels);
xlabel('Region');
ylabel('Accuracy (%)');
title('Patient Accuracies by Region (n=33)');
ylim([0 100])
xlim([0 nRegions+1])
grid on; box off;
end


function groupRegionPlot
nPatients = numel(groupedRegionAccuracies);
nRegions = size(groupedRegionAccuracies{1}, 1);
RegionLabels = {'dlPFC', 'dmPFC + dACC', 'OFC', 'vlPFC', 'STG + MTG + ITG', 'AMY', 'HIP'};

allAcc = nan(nPatients, nRegions);
allErr = nan(nPatients, nRegions);

for i = 1:nPatients
    temp = groupedRegionAccuracies{i};
    temp(temp == 0) = NaN; % ignore zeros
    allAcc(i, :) = temp(:, 1)';
    allErr(i, :) = temp(:, 2)';
end

figure;

%%%%%% Choose Plot Type by uncommmenting: %%%%%%%%

% % Bar
% regionAcc = mean(allAcc, 'omitnan');
% regionErr = std(allAcc, 'omitnan') ./ sqrt(sum(~isnan(allAcc)));
% bar(regionAcc, 'FaceColor', [0.2 0.6 0.8]); hold on;
% errorbar(1:nRegions, regionAcc, regionErr, 'k.', 'LineWidth', 1.2);

% % Boxplot
% boxplot(allAcc, 'Labels', RegionLabels); hold on;

% % Scatter Plot
rng(1); % for reproducibility
jitterAmount = 0.1;
x=repmat(1:length(RegionLabels),size(allAcc,1),1);
for j=1:size(x,2)
    x(:,j) = x(:,j) + jitterAmount * randn(size(x(:,j)));
    hold on;
    scatter(x(:,j), allAcc(:,j), 25, 'r', 'filled');
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
set(gca, 'XTick', 1:nRegions, 'XTickLabel', RegionLabels);
xlabel('Region');
ylabel('Accuracy (%)');
title('Patient Accuracies by Grouped Region (n=33)');
ylim([0 100])
xlim([0 nRegions+1])
grid on; box off;


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