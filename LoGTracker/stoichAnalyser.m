function [meanStoich, medianStoich, modeStoich, SD, SEM, stoich_all, iSingle] = stoichAnalyser(Spots, varargin)
% STOICHANALYSER calculates stoichiometry statistics and displays KDE distribution.

    %% INPUT PARSING
    mask = [];
    params = struct();
    if nargin >= 2, mask = varargin{1}; end
    if nargin >= 3, params = varargin{2}; end

    % Set default parameters
    params = setDefault(params, 'iSingle', []);
    params = setDefault(params, 'lateFrameFraction', 0.8);
    params = setDefault(params, 'nFrames', 5);
    params = setDefault(params, 'minFrames', 5);
    params = setDefault(params, 'threshold', 150);
    params = setDefault(params, 'showPlot', 1);
    params = setDefault(params, 'kdeBandwidth', 3);

    %% FILTERING
    if ~isempty(mask)
        Spots = filterSpotsByMask(Spots, mask);
    end

    if isempty(Spots)
        warning('No spots to analyse');
        [meanStoich, medianStoich, modeStoich, SD, SEM, iSingle] = deal(NaN);
        stoich_all = [];
        return;
    end

    %% CALCULATION
    % Calls external function for linear fit analysis
    [~, iSingle, stoich_all] = generalLinearStoich(Spots, ...
        'iSingle', params.iSingle, ...
        'lateFrameFraction', params.lateFrameFraction, ...
        'nFrames', params.nFrames, ...
        'minFrames', params.minFrames, ...
        'threshold', params.threshold);

    if isempty(stoich_all)
        warning('No valid trajectories found');
        [meanStoich, medianStoich, modeStoich, SD, SEM] = deal(NaN);
        return;
    end

    %% STATISTICS
    meanStoich = mean(stoich_all);
    medianStoich = median(stoich_all);
    SD = std(stoich_all);
    SEM = SD / sqrt(length(stoich_all));

    % Determine mode from Kernel Density Estimate peak
    [kdeY, kdeX] = ksdensity(stoich_all, 'Bandwidth', params.kdeBandwidth);
    [~, maxIdx] = max(kdeY);
    modeStoich = kdeX(maxIdx);

    % Console output
    fprintf('Stoichiometry Analysis (n=%d)\n', length(stoich_all));
    fprintf('  Mean: %.2f | Median: %.2f | Mode: %.2f\n', meanStoich, medianStoich, modeStoich);
    fprintf('  SD: %.2f | SEM: %.2f | iSingle: %.1f\n', SD, SEM, iSingle);

    %% VISUALISATION
    if params.showPlot
        renderKdePlot(kdeX, kdeY, meanStoich, medianStoich, modeStoich, length(stoich_all));
    end
end

%% HELPER FUNCTIONS

function p = setDefault(p, field, value)
    % Assigns default value if field does not exist
    if ~isfield(p, field), p.(field) = value; end
end

function renderKdePlot(x, y, avg, med, mod, n)
    % Generates probability density plot with statistical markers
    figure('Name', 'Stoichiometry Distribution', 'NumberTitle', 'off');
    plot(x, y, 'b-', 'LineWidth', 2);
    hold on;
    
    yLim = max(y) * 1.1;
    plot([avg, avg], [0, yLim], 'r--', 'LineWidth', 1.5);
    plot([med, med], [0, yLim], 'g--', 'LineWidth', 1.5);
    plot([mod, mod], [0, yLim], 'm--', 'LineWidth', 1.5);
    
    xlabel('Stoichiometry');
    ylabel('Probability Density');
    title(sprintf('Stoichiometry Distribution (n=%d)', n));
    legend('KDE', sprintf('Mean=%.2f', avg), sprintf('Median=%.2f', med), sprintf('Mode=%.2f', mod));
    xlim([0, max(x)]);
    ylim([0, yLim]);
    hold off;
end

function SpotsFiltered = filterSpotsByMask(Spots, mask)
    % Keeps trajectories whose first detected spot is within the mask
    Spots(Spots(:,10) == 0, :) = [];
    if isempty(Spots), SpotsFiltered = []; return; end
    
    trajIDs = unique(Spots(:,10));
    keepFlags = false(size(trajIDs));
    
    for i = 1:length(trajIDs)
        tSpots = Spots(Spots(:,10) == trajIDs(i), :);
        [~, firstIdx] = min(tSpots(:, 9));
        
        x = round(tSpots(firstIdx, 1));
        y = round(tSpots(firstIdx, 2));
        
        % Check image boundaries and mask value
        if x >= 1 && x <= size(mask, 2) && y >= 1 && y <= size(mask, 1)
            if mask(y, x) > 0, keepFlags(i) = true; end
        end
    end
    
    SpotsFiltered = Spots(ismember(Spots(:,10), trajIDs(keepFlags)), :);
end