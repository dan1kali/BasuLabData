%% Beta weights plot

subject = {... 
'BW42', 'MG51b', 'MG79', 'MG86', ...
'MG89', 'MG90', 'MG91', 'MG95', ...
'MG96', 'MG99', 'MG102', 'MG104', ...
'MG105', 'MG106', 'MG111', 'MG112',...
'MG116', 'MG117', 'MG118', 'MG120',...
'UCMC01', 'UCMC02', 'UCMC03', 'UCMC05', ...
'UCMC06', 'UCMC07', 'UCMC08', 'UCMC09', ...
'UCMC11', 'UCMC13', 'UCMC14', 'UCMC15', 'UCMC17',...
};


RegionLabels = {'dlPFC', 'dmPFC', 'OFC', 'vlPFC', 'STG', 'MTG', 'ITG', 'dACC', 'AMY', 'HIP'};

    beta_mean_plot = zeros(10,length(subject));
    beta_max_plot = zeros(10,length(subject));
    beta_mean_err_plot = zeros(10,length(subject));
    beta_max_err_plot = zeros(10,length(subject));
        
    for i_sub = 1:length(subject)        
        nROI = numel(ROIchannels{i_sub});
    
        mean_weights = abs(meanBetas{i_sub});
        max_weights = abs(maxBetas{i_sub});
    
        mean_abs_beta_flat = reshape(mean_weights, size(mean_weights,1), []);  % now size = [60 x 500]
        max_abs_beta_flat = reshape(max_weights, size(max_weights,1), []);
        
        beta_mean = squeeze(mean(mean_weights,[2 3]));
        beta_max = squeeze(mean(max_weights,[2 3]));
        
        beta_mean_err = std(mean_abs_beta_flat,0,2)/sqrt(size(mean_abs_beta_flat, 2)); 
        beta_max_err = std(max_abs_beta_flat,0,2)/sqrt(size(max_abs_beta_flat, 2)); 
        
        
        for iROI = 1:nROI
            
            subjectROI = ROIchannels{i_sub};
            % channelsinROI = find(ismember(allTaskChans{isub}, subjectROI{iROI}));
            channelsinROI = subjectROI{iROI};

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


%% generate shap_vals_plot

subject = {... 
'BW42', 'MG51b', 'MG79', 'MG86', ...
'MG89', 'MG90', 'MG91', 'MG95', ...
'MG96', 'MG99', 'MG102', 'MG104', ...
'MG105', 'MG106', 'MG111', 'MG112',...
'MG116', 'MG117', 'MG118', 'MG120',...
'UCMC01', 'UCMC02', 'UCMC03', 'UCMC05', ...
'UCMC06', 'UCMC07', 'UCMC08', 'UCMC09', ...
'UCMC11', 'UCMC13', 'UCMC14', 'UCMC15', 'UCMC17',...
};


RegionLabels = {'dlPFC', 'dmPFC', 'OFC', 'vlPFC', 'STG', 'MTG', 'ITG', 'dACC', 'AMY', 'HIP'};

    shap_vals_plot = zeros(10,length(subject));
   
    for i_sub = 1:length(subject)        
        nROI = numel(ROIchannels{i_sub});
    
        shap_vals_orig = abs(shapVals{i_sub});
            
        shap_vals = squeeze(mean(shap_vals_orig,[2 3]));
             
        for iROI = 1:nROI
            
            subjectROI = ROIchannels{i_sub};
            % channelsinROI = find(ismember(allTaskChans{isub}, subjectROI{iROI}));
            channelsinROI = subjectROI{iROI};
            % each channel in ROI corresponds to 2 rows in the shap matrix, which
            % has mean and max features for every channel
            shapchannelsinROI = reshape([2*channelsinROI-1; 2*channelsinROI], 1, []);

            if isempty(channelsinROI)
                shap_vals_plot (iROI,i_sub) = 0;
            else
                shap_vals_plot (iROI,i_sub) = mean(shap_vals(shapchannelsinROI));
            end
        end
    end

%% Shapley Bar Graph

load('c_plotData.mat');
shap_vals_plot_orig = shap_vals_plot;
shap_vals_plot = mean(shap_vals_plot,2); % comment out for grouped bars

[nBars, nGroups] = size(shap_vals_plot);

shap_err = std(shap_vals_plot_orig, 0, 2) ./ sqrt(nGroups);  % Standard error

figure;
b1 = bar(shap_vals_plot, 'grouped','BarWidth', 1);
xlim([0 (nBars+1)]); xticks(1:nBars); xticklabels(RegionLabels);

x = nan(nBars, nGroups);
for i = 1:nGroups
    x(:,i) = b1(i).XEndPoints;
end
hold on;
er = errorbar(x, shap_vals_plot, shap_err); 
for ier = 1:numel(er)
    er(ier).LineStyle = 'none'; er(ier).CapSize = 5; er(ier).Color = [0 0 0];
end

set(gca, 'FontSize', 10);
xlabel('Region of Interest', 'FontSize', 12);
ylabel('Shap Values', 'FontSize', 12);
title('Shapley Values by Region', 'FontSize', 14);
grid on;

%% Shapley Heatmap

RegionLabels = {'dlPFC', 'dmPFC', 'OFC', 'vlPFC', 'STG', 'MTG', 'ITG', 'dACC', 'AMY', 'HIP'};
load('c_plotData.mat');

shap_vals_plot_nan = shap_vals_plot;
shap_vals_plot_nan(shap_vals_plot_nan == 0) = NaN;

figure;
imagesc(shap_vals_plot_nan) %%% change to heatmap %%%
colorbar
xlabel('Patient')
ylabel('Regions')
yticks(1:10)
xticks(1:33)
yticklabels(RegionLabels)
xticklabels(subject)
title('SHAP values','FontSize', 14);


%% Shapley Boxplot


load('c_plotData.mat'); % shap_vals_plot_orig: 10 x 33
% Rows = regions, Columns = patients
[nRegions, nPatients] = size(shap_vals_plot);

figure;
boxplot(shap_vals_plot', 'Labels', RegionLabels); 

ylabel('SHAP Values', 'FontSize', 12);
xlabel('Region of Interest', 'FontSize', 12);
title('Shapley Values by Region', 'FontSize', 14);
grid on;
yline(0,'k--');
ylim([0 0.45])


%% Shapley Jitter scatter

load('c_plotData.mat'); % shap_vals_plot_orig: 10 x 33
[nRegions, nPatients] = size(shap_vals_plot);

figure; hold on;

rng(1); 
jitterAmount = 0.1;

for r = 1:nRegions
    x_jitter = r + jitterAmount*randn(1, nPatients); 
    scatter(x_jitter, shap_vals_plot(r,:), 50, 'filled');
end

xlabel('Region of Interest', 'FontSize', 12);
ylabel('SHAP Values', 'FontSize', 12);
title('SHAP Values per Region, n=33', 'FontSize', 14);
xticks(1:nRegions);
xticklabels(RegionLabels);
grid on;
yline(0,'k--');
set(gca, 'FontSize', 10);
hold off;
xlim([0 11])



%% combined 

subject = {... 
'BW42', 'MG51b', 'MG79', 'MG86', ...
'MG89', 'MG90', 'MG91', 'MG95', ...
'MG96', 'MG99', 'MG102', 'MG104', ...
'MG105', 'MG106', 'MG111', 'MG112',...
'MG116', 'MG117', 'MG118', 'MG120',...
'UCMC01', 'UCMC02', 'UCMC03', 'UCMC05', ...
'UCMC06', 'UCMC07', 'UCMC08', 'UCMC09', ...
'UCMC11', 'UCMC13', 'UCMC14', 'UCMC15', 'UCMC17',...
};


RegionLabels = {'dlPFC', 'dmPFC', 'OFC', 'vlPFC', 'STG', 'MTG', 'ITG', 'dACC', 'AMY', 'HIP'};

    beta_mean_plot = zeros(10,length(subject));
    beta_max_plot = zeros(10,length(subject));
    beta_mean_err_plot = zeros(10,length(subject));
    beta_max_err_plot = zeros(10,length(subject));
        
    for i_sub = 1:length(subject)        
        nROI = numel(ROIchannels{i_sub});
    
        mean_weights = abs(meanBetas{i_sub});
        max_weights = abs(maxBetas{i_sub});
    
        mean_abs_beta_flat = reshape(mean_weights, size(mean_weights,1), []);  % now size = [60 x 500]
        max_abs_beta_flat = reshape(max_weights, size(max_weights,1), []);
        
        beta_mean = squeeze(mean(mean_weights,[2 3]));
        beta_max = squeeze(mean(max_weights,[2 3]));
        
        beta_mean_err = std(mean_abs_beta_flat,0,2)/sqrt(size(mean_abs_beta_flat, 2)); 
        beta_max_err = std(max_abs_beta_flat,0,2)/sqrt(size(max_abs_beta_flat, 2)); 
        
        
        for iROI = 1:nROI
            
            subjectROI = ROIchannels{i_sub};
            % channelsinROI = find(ismember(allTaskChans{isub}, subjectROI{iROI}));
            channelsinROI = subjectROI{iROI};

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

beta_mean_plot = mean(beta_mean_plot,2);
beta_mean_err_plot = mean(beta_mean_err_plot,2);
beta_max_plot = mean(beta_max_plot,2);
beta_max_err_plot = mean(beta_max_err_plot,2);


[nBars, nGroups] = size(beta_mean_plot);

figure;
tiledlayout(1, 3, 'TileSpacing', 'compact');

meanplot = nexttile;
b1 = bar(beta_mean_plot, 'grouped','BarWidth', 1);
xlim([0 (nBars+1)]); xticks(1:nBars); xticklabels(RegionLabels);


set(gca, 'FontSize', 10);
xlabel('Region of Interest', 'FontSize', 12);
ylabel('Absolute Beta Weight', 'FontSize', 12);
title('Mean Beta Weights', 'FontSize', 14);
grid on;

maxplot = nexttile;
b2 = bar(beta_max_plot, 'grouped','BarWidth', 1);
xlim([0 (nBars+1)]); xticks(1:nBars); xticklabels(RegionLabels);

xlabel('Region of Interest', 'FontSize', 12);
title('Max Beta Weights', 'FontSize', 14);
grid on;

if nGroups ~= 1
    legend(b1, subject, 'Location', 'eastoutside');
else
end


shapplot = nexttile;
load('c_plotData.mat');

shap_vals_plot_orig = shap_vals_plot;
shap_vals_plot = mean(shap_vals_plot,2);

[nBars, nGroups] = size(shap_vals_plot);

shap_err = std(shap_vals_plot_orig, 0, 2) ./ sqrt(nGroups);  % Standard error

b3 = bar(shap_vals_plot, 'grouped','BarWidth', 1);
xlim([0 (nBars+1)]); xticks(1:nBars); xticklabels(RegionLabels);


set(gca, 'FontSize', 10);
xlabel('Region of Interest', 'FontSize', 12);
ylabel('Shap Values', 'FontSize', 12);
title('Shapley Values by Region', 'FontSize', 14);
grid on;

linkaxes([meanplot, maxplot,shapplot], 'y');
sgtitle('Feature Contributions Across Regions', 'FontSize', 16);
