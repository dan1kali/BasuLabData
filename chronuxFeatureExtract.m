
% time_indices = find(ft_data3.time >= 0 & t <= response_time);


%% 

% Need to preprocess first

% Then need to cut into congruent and incongruent

% then run for every channel/trial? 
% In all the decoding analyses, a support vector machine (SVM) classifier 
% with a linear kernel was trained after concatenating all the conflict 
% modulated electrodes in each task.


%% Extract max/mean power

nTrials = length(ft_data3.trial);
nChannels = length(ft_data3.trial{1}(:,1));

band_power = zeros(( nTrials .* nChannels ) ,2); % column 1 mean, column 2 max

for c=1:nChannels
    for tr=1:nTrials
        [y,t,S,f,Snorm]=computeNormalizedFreqMag(ft_data3.trial{tr}(c,:),1000);
    
        t = t + ft_data3.time{1}(1);  % shift t so its zero corresponds to 0 s in original time
        
        % time_indices = find(t >= 0 & t <= response_time);
        time_indices = find(t >= 0);
        
        S_epoched = Snorm(time_indices, :);  
    
    concat_index = ((c-1)*nTrials)+ tr; % concatenate all electrodes
    band_power(concat_index,1) = mean(S_epoched(:));  % mean over time
    band_power(concat_index,2) = max(S_epoched(:));    % max over time
    end
    % band_power = [mean power, max power] for all trials per channel
end

% per one CHANNEL
save('BW42_band_power.mat', 'band_power');

avg_band_power = mean(band_power, 1);  % [average mean power, average max power]
sem_band_power = std(band_power, 0, 1) ./ sqrt(size(band_power, 1));  % SEM for each

%% Plot bar graph

bar([1 2], avg_band_power); hold on
errorbar([1 2], avg_band_power, sem_band_power, 'k', 'LineStyle', 'none', 'LineWidth', 1.5);

xticks([1 2]);
xticklabels({'Mean Power', 'Max Power'});
ylabel('Normalized Power');
title('Average Band Power with SEM');
legend({'Average Power', 'SEM Error Bars'}, 'Location', 'best');

hold off

fprintf('mean power = %.2f\n max power = %.2f\n',avg_band_power(1), avg_band_power(2));


%% functions

function [y,t,S,f,Snorm]=computeNormalizedFreqMag(x,Fs,params)

if(nargin < 4)
    params.Fs = Fs;
end
if(~isfield(params,'Fs'))
    params.Fs = Fs;
end
if(~isfield(params,'fpass'))
    params.fpass = [70 120];
end
if(~isfield(params,'tapers'))
    params.tapers = [5 7];
end
if(~isfield(params,'movingwin'))
    params.movingwin = [0.2 0.01]; % [window size, step size]
end
if ~isfield(params,'norm_option')
    params.norm_option='mean';
end
if ~isfield(params,'err')
    params.err=[1 0.05];
end

[S,t,f,~] = mtspecgramc(x,params.movingwin,params); % Chronux function
% units of spectrum? 

if strcmp(params.norm_option,'median')
    m=median(S,1);
elseif strcmp(params.norm_option,'mean')
    m=mean(S,1);
end


Snorm = S./repmat(m, size(S,1), 1);  % normalized power (this is implicit in your code)

y=mean(Snorm,2); % mean power for each trial and channel

% get max and min of S

end