function [stoichiometry, iSingle, stoich_all] = generalLinearStoich(Spots, varargin)
% GENERALLINEARSTOICH calculates stoichiometry by extrapolating intensity to t=0.

    %% PARSE PARAMETERS
    p = inputParser;
    p.addParameter('iSingle', []);
    p.addParameter('lateFrameFraction', 0.8);
    p.addParameter('nFrames', 5);
    p.addParameter('minFrames', 5);
    p.addParameter('threshold', 150);
    p.parse(varargin{:});
    res = p.Results;

    %% CALCULATE SINGLE MOLECULE INTENSITY
    % iSingle is derived from late frames where trajectories are likely single molecules
    iSingle = res.iSingle;
    if isempty(iSingle)
        maxFrame = max(Spots(:, 9));
        lateCutoff = maxFrame * res.lateFrameFraction;
        lateData = Spots(Spots(:, 9) > lateCutoff & Spots(:, 10) > 0, 5);
        iSingle = mean(lateData);
    end

    %% ANALYZE TRAJECTORIES
    trajIDs = unique(Spots(Spots(:, 10) > 0, 10));
    raw_results = zeros(length(trajIDs), 1);
    validCount = 0;

    for i = 1:length(trajIDs)
        % Extract data for current trajectory
        trajData = Spots(Spots(:, 10) == trajIDs(i), :);
        nAvailable = size(trajData, 1);
        nFit = min(res.nFrames, nAvailable);

        if nFit >= res.minFrames
            % Calculate time relative to the start of imaging
            relTime = trajData(1:nFit, 9) - trajData(1, 12);
            intensity = trajData(1:nFit, 5);

            % Perform linear fit to find initial intensity (intercept)
            coeffs = polyfit(relTime, intensity, 1);
            slope = coeffs(1);
            intercept = coeffs(2);
            val = intercept / iSingle;

            % Apply quality filters: non-increasing intensity and within stoichiometry bounds
            if val > 0 && val < res.threshold && slope <= 0
                validCount = validCount + 1;
                raw_results(validCount) = val;
            end
        end
    end

    %% FINAL OUTPUTS
    stoich_all = raw_results(1:validCount);
    if validCount > 0
        stoichiometry = mean(stoich_all);
    else
        stoichiometry = NaN;
    end
end