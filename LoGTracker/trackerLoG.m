function [Spots, frame_average, p, meta_data, image_data] = trackerLoG(fileName, p)
% TRACKERLOG performs spot tracking using LoG detection and LAP linking.

    %% SET DEFAULT PARAMETERS
    if ~exist('p','var') || isempty(p), p = struct(); end

    % Detection and Localisation
    p = setDefault(p, 'DetectionMethod', 'LoG');
    p = setDefault(p, 'disk_radius', 5);
    p = setDefault(p, 'gaussian', 0);
    p = setDefault(p, 'Candidate_d_min', 0);
    p = setDefault(p, 'LocalisationMethod', 'iterative');
    p = setDefault(p, 'subarray_halfwidth', 8);
    p = setDefault(p, 'inner_circle_radius', 5);
    p = setDefault(p, 'gauss_mask_sigma', 2);
    p = setDefault(p, 'error_set', 0.05);
    p = setDefault(p, 'SNR_min', 0.4);
    p = setDefault(p, 'sigmaFit_min', 0);
    p = setDefault(p, 'sigmaFit_max', p.inner_circle_radius);

    % Linking
    p = setDefault(p, 'LinkingMethod', 'LAP');
    p = setDefault(p, 'd_01_max', 5);
    p = setDefault(p, 'Iratio_01_min', 0.5);
    p = setDefault(p, 'Iratio_01_max', 3);
    p = setDefault(p, 'gap_frames', 3);
    p = setDefault(p, 'd_min', 1);

    % Frame selection and Output
    p = setDefault(p, 'all', 1);
    p = setDefault(p, 'startFrame', 1);
    p = setDefault(p, 'endFrame', 999);
    p = setDefault(p, 'noFrames', 5);
    p = setDefault(p, 'DetermineFirstFrames', 0);
    p = setDefault(p, 'firstFrame', 1);
    p = setDefault(p, 'FramesToTrack', 0);
    p = setDefault(p, 'satPixelVal', 10^10);
    p = setDefault(p, 'show_output', 0);
    p = setDefault(p, 'show_text_output', 1);

    %% LOAD DATA
    if isa(fileName, 'char')
        if p.all == 1
            [image_data, meta_data] = imEx1(fileName);
        else
            [image_data, meta_data] = imEx1(fileName, p.startFrame, p.endFrame);
        end
    else
        image_data = fileName;
        meta_data = 'Image data passed directly';
    end
    numFrames = size(image_data, 3);

    %% DEFINE TRACKING RANGE
    if p.DetermineFirstFrames == 1
        [firstFrame, ~, ~, ~] = LaserOn3(image_data, 0, 0);
    else
        firstFrame = p.firstFrame;
    end

    if p.FramesToTrack == 0
        p.FramesToTrack = numFrames - firstFrame;
    end

    startFrame = firstFrame;
    endFrame = min(firstFrame + p.FramesToTrack, numFrames);

    % Calculate background/average frame
    avgIdx = firstFrame:min(firstFrame + p.noFrames - 1, numFrames);
    frame_average = mean(image_data(:,:,avgIdx), 3);

    %% MAIN TRACKING LOOP
    if p.show_text_output
        fprintf('\nTracking: Detection=%s, Localisation=%s, Linking=%s\n', ...
            p.DetectionMethod, p.LocalisationMethod, p.LinkingMethod);
    end

    all_spots = cell(endFrame - startFrame + 1, 1);
    tic;

    for i = startFrame:endFrame
        if p.show_text_output, fprintf('Frame %d/%d\r', i, endFrame); end
        
        frame = image_data(:,:,i);
        if max(frame(:)) > p.satPixelVal, continue; end
        
        % Detection
        if strcmpi(p.DetectionMethod, 'log')
            result = findSpotsLoG_internal(frame, p.disk_radius, p.show_output);
        else
            result = findSpotsOriginal_internal(frame, p.disk_radius, p.gaussian, p.show_output);
        end
        
        [y_est, x_est] = find(result);
        if isempty(x_est), continue; end
        if p.Candidate_d_min > 0
            [y_est, x_est] = mergeCandidates(y_est, x_est, p.Candidate_d_min);
        end
        
        % Localisation
        frame_spots = zeros(length(x_est), 12);
        valid_count = 0;
        
        for j = 1:length(x_est)
            if strcmpi(p.LocalisationMethod, 'quadratic')
                [xc, yc, clip, Ibg, Isp, ~, bg_std, m_pix, noConv] = ...
                    gaussianMaskQuadratic_internal(frame, x_est(j), y_est(j), p);
            else
                [xc, yc, clip, Ibg, Isp, ~, bg_std, m_pix, noConv] = ...
                    gaussianMaskIterative_internal(frame, x_est(j), y_est(j), p);
            end
            
            snr = Isp / (bg_std * m_pix + eps);
            
            if ~noConv && ~clip && snr > p.SNR_min
                [sdx, sdy, Icent] = estimatePSFWidth(frame, xc, yc, Ibg, p);
                valid_count = valid_count + 1;
                frame_spots(valid_count,:) = [xc, yc, clip, Ibg, Isp, sdx, sdy, Icent, i, 0, snr, startFrame];
            end
        end
        
        frame_spots = frame_spots(1:valid_count, :);
        if ~isempty(frame_spots) && p.d_min > 0
            frame_spots = mergeCoincidentSpots(frame_spots, p.d_min);
        end
        all_spots{i - startFrame + 1} = frame_spots;
    end

    Spots = cell2mat(all_spots);

    %% LINKING
    if strcmpi(p.LinkingMethod, 'greedy') && ~isempty(Spots)
        % Greedy linking is performed sequentially on the resulting array
        frames = unique(Spots(:,9));
        for i = 2:length(frames)
            Spots = linkSpotsGreedy_internal(Spots, frames(i), frames(i-1), p);
        end
    elseif strcmpi(p.LinkingMethod, 'LAP') && ~isempty(Spots)
        Spots = linkSpotsLAP_internal(Spots, p);
    end

    if p.show_text_output
        fprintf('\nTracking complete in %.2f seconds. Tracks: %d\n', toc, max(Spots(:,10)));
    end
end

function p = setDefault(p, field, value)
    if ~isfield(p, field), p.(field) = value; end
end

%% INTERNAL HELPER FUNCTIONS

function result = findSpotsLoG_internal(frame, disk_radius, show_output)
    % Laplacian of Gaussian spot detection
    frame = double(frame);
    sigma = max(1, disk_radius / 2.5);
    kSize = max(5, ceil(6 * sigma));
    if mod(kSize, 2) == 0, kSize = kSize + 1; end
    
    filtered = -imfilter(frame, fspecial('log', kSize, sigma), 'replicate');
    
    % Adaptive threshold based on background statistics
    medV = median(filtered(:));
    stdV = std(filtered(:));
    bg = filtered(filtered < medV + 2*stdV);
    thresh = median(bg) + 5 * std(bg);
    
    result = imregionalmax(filtered) & (filtered > thresh);
    
    % Clear borders
    b = ceil(disk_radius);
    result([1:b, end-b+1:end], :) = 0;
    result(:, [1:b, end-b+1:end]) = 0;
end

function result = findSpotsOriginal_internal(frame, disk_radius, gaussian, show_output)
    % COMBINED DETECTION: Merges Background Masking (findSpots2LoG) and Histogram (findSpots3LoG).
    result1 = findSpots2LoG(frame, 2, disk_radius, gaussian, 0);
    result2 = findSpots3LoG(frame, 2, disk_radius, gaussian, 0);
    result_combined = result1 + result2;
    % Ensure binary output (where 1 means found in either or both methods)
    result_combined(result_combined > 1) = 1;
    % Reduces any multi-pixel clusters to a single center point
    result = bwmorph(result_combined, 'thin', Inf);
end

function [xc, yc, clip, Ibg, Isp, Idata, bg_std, m_pix, noConv] = ...
    gaussianMaskQuadratic_internal(img, x_est, y_est, p)
    % Localisation using quadratic fit to local peak
    d = p.subarray_halfwidth;
    r = p.inner_circle_radius;
    
    % Subarray extraction
    x_est = max(d+1, min(size(img,2)-d, round(x_est)));
    y_est = max(d+1, min(size(img,1)-d, round(y_est)));
    I = double(img(y_est-d:y_est+d, x_est-d:x_est+d));
    Idata = I;
    
    [Xs, Ys] = meshgrid(-d:d, -d:d);
    inner = (Xs.^2 + Ys.^2) <= r^2;
    Ibg = mean(I(~inner));
    bg_std = std(I(~inner));
    I2 = I - Ibg;
    Isp = sum(I2(inner));
    m_pix = sum(inner(:));
    
    [~, maxIdx] = max(I2(:));
    [pr, pc] = ind2sub(size(I2), maxIdx);
    
    clip = 0; noConv = 0;
    if pc > 1 && pc < size(I2,2) && pr > 1 && pr < size(I2,1)
        % Quadratic interpolation for sub-pixel peak
        dx = 0.5 * (I2(pr, pc-1) - I2(pr, pc+1)) / (I2(pr, pc-1) - 2*I2(pr, pc) + I2(pr, pc+1) + eps);
        dy = 0.5 * (I2(pr-1, pc) - I2(pr+1, pc)) / (I2(pr-1, pc) - 2*I2(pr, pc) + I2(pr+1, pc) + eps);
        xc = x_est + (pc - d - 1) + max(-0.5, min(0.5, dx));
        yc = y_est + (pr - d - 1) + max(-0.5, min(0.5, dy));
    else
        xc = x_est; yc = y_est; clip = 1; noConv = 1;
    end
end

function [xc, yc, clip, Ibg, Isp, Idata, bg_std, m_pix, noConv] = ...
    gaussianMaskIterative_internal(img, x_est, y_est, p)
    % Localisation using iterative Gaussian-weighted centroid
    d = p.subarray_halfwidth;
    r = p.inner_circle_radius;
    sig = p.gauss_mask_sigma;
    
    x_est = max(d+1, min(size(img,2)-d, round(x_est)));
    y_est = max(d+1, min(size(img,1)-d, round(y_est)));
    I = double(img(y_est-d:y_est+d, x_est-d:x_est+d));
    Idata = I;
    [Xs, Ys] = meshgrid(x_est-d:x_est+d, y_est-d:y_est+d);
    
    xc = x_est; yc = y_est; clip = 0; noConv = 0;
    
    for k = 1:300
        distSq = (Xs - xc).^2 + (Ys - yc).^2;
        inner = distSq <= r^2;
        Ibg = mean(I(~inner));
        bg_std = std(I(~inner));
        I2 = I - Ibg;
        Isp = sum(I2(inner));
        m_pix = sum(inner(:));
        
        mask = exp(-distSq / (2*sig^2));
        I3 = I2 .* (mask / sum(mask(:)));
        
        xn = sum(I3(:) .* Xs(:)) / sum(I3(:));
        yn = sum(I3(:) .* Ys(:)) / sum(I3(:));
        
        if sqrt((xn-xc)^2 + (yn-yc)^2) <= p.error_set && k >= 5, break; end
        xc = xn; yc = yn;
        
        if abs(xc-x_est) > (d-r+1) || abs(yc-y_est) > (d-r+1), clip = 1; break; end
        if k == 300, noConv = 1; end
    end
end

function [sdx, sdy, Icent] = estimatePSFWidth(frame, xc, yc, Ibg, p)
    % Estimates PSF width from cross-section profiles
    d = min(p.subarray_halfwidth, 5);
    x = max(d+1, min(size(frame,2)-d, round(xc)));
    y = max(d+1, min(size(frame,1)-d, round(yc)));
    
    xp = double(frame(y, x-d:x+d)) - Ibg;
    yp = double(frame(y-d:y+d, x)) - Ibg;
    
    sdx = estimateSigma(xp);
    sdy = estimateSigma(yp);
    Icent = double(frame(y, x)) - Ibg;
end

function s = estimateSigma(profile)
    p = max(0, profile(:)');
    [mv, mi] = max(p);
    l = find(p(1:mi) < mv/2, 1, 'last');
    r = find(p(mi:end) < mv/2, 1, 'first') + mi - 1;
    if isempty(l), l = 1; end
    if isempty(r), r = length(p); end
    s = max(0.5, min(5, (r-l)/2.355));
end

function spots = linkSpotsGreedy_internal(spots, f1, f0, p)
    % Basic frame-to-frame proximity linking
    s1 = spots(spots(:,9) == f1, :);
    s0 = spots(spots(:,9) == f0, :);
    if isempty(s1) || isempty(s0), return; end
    
    dist = pdist2(s1(:,1:2), s0(:,1:2));
    [mv, mi] = min(dist, [], 2);
    
    for j = 1:size(s1,1)
        if mv(j) <= p.d_01_max
            ratio = s0(mi(j), 5) / (s1(j, 5) + eps);
            if ratio >= p.Iratio_01_min && ratio <= p.Iratio_01_max
                if s0(mi(j), 10) == 0, s0(mi(j), 10) = max(spots(:,10)) + 1; end
                s1(j, 10) = s0(mi(j), 10);
            end
        end
    end
    spots(spots(:,9) == f1, :) = s1;
    spots(spots(:,9) == f0, :) = s0;
end

function spots = linkSpotsLAP_internal(spots, p)
    % Linear Assignment Problem linking with gap closing
    frames = sort(unique(spots(:,9)));
    spots(:,10) = 0;
    nextID = 1;
    
    for i = 2:length(frames)
        idx1 = find(spots(:,9) == frames(i-1));
        idx2 = find(spots(:,9) == frames(i));
        if isempty(idx1) || isempty(idx2), continue; end
        
        D = pdist2(spots(idx2,1:2), spots(idx1,1:2));
        D(D > p.d_01_max) = Inf;
        
        % Filter by intensity ratio
        I1 = spots(idx1, 5); I2 = spots(idx2, 5);
        for a = 1:size(D,1)
            for b = 1:size(D,2)
                r = I1(b) / (I2(a) + eps);
                if r < p.Iratio_01_min*0.5 || r > p.Iratio_01_max*2, D(a,b) = Inf; end
            end
        end
        
        M = greedyAssignment(D);
        for k = 1:size(M,1)
            if D(M(k,1), M(k,2)) < Inf
                if spots(idx1(M(k,2)), 10) == 0
                    spots(idx1(M(k,2)), 10) = nextID;
                    nextID = nextID + 1;
                end
                spots(idx2(M(k,1)), 10) = spots(idx1(M(k,2)), 10);
            end
        end
    end
    
    % Finalise unlinked and close gaps
    unlinked = spots(:,10) == 0;
    spots(unlinked, 10) = nextID:(nextID + sum(unlinked) - 1);
    if p.gap_frames > 0, spots = closeGaps_internal(spots, p); end
    [~, ~, spots(:,10)] = unique(spots(:,10));
end

function M = greedyAssignment(D)
    [n2, n1] = size(D);
    M = []; free = true(1, n1);
    for j = 1:n2
        [mv, mi] = min(D(j, free));
        if mv < Inf
            idx = find(free);
            real_idx = idx(mi);
            M = [M; j, real_idx];
            free(real_idx) = false;
        end
    end
end

function spots = closeGaps_internal(spots, p)
    % Connects trajectory ends to starts within frame and distance limits
    traj = unique(spots(:,10));
    n = length(traj);
    if n < 2, return; end
    
    ends = zeros(n, 4); starts = zeros(n, 4);
    for i = 1:n
        tS = spots(spots(:,10) == traj(i), :);
        [~, li] = max(tS(:,9)); [~, fi] = min(tS(:,9));
        ends(i,:) = [traj(i), tS(li, 9), tS(li, 1:2)];
        starts(i,:) = [traj(i), tS(fi, 9), tS(fi, 1:2)];
    end
    
    gapCost = Inf(n, n);
    for i = 1:n
        for j = 1:n
            gap = starts(j,2) - ends(i,2);
            if gap > 0 && gap <= p.gap_frames + 1
                d = sqrt(sum((ends(i,3:4) - starts(j,3:4)).^2));
                if d < p.d_01_max * sqrt(gap + 1) * 2
                    gapCost(i,j) = d + gap * 0.3;
                end
            end
        end
    end
    
    [~, sIdx] = sort(gapCost(:));
    [ei, si] = ind2sub(size(gapCost), sIdx);
    merged = false(n, 1);
    for k = 1:length(sIdx)
        if gapCost(sIdx(k)) == Inf, break; end
        if ~merged(ei(k)) && ~merged(si(k))
            spots(spots(:,10) == starts(si(k), 1), 10) = ends(ei(k), 1);
            merged(si(k)) = true; merged(ei(k)) = true;
        end
    end
end

function [yo, xo] = mergeCandidates(yi, xi, d_min)
    keep = true(length(xi), 1);
    for i = 1:length(xi)
        if ~keep(i), continue; end
        for j = i+1:length(xi)
            if keep(j) && sqrt((xi(i)-xi(j))^2 + (yi(i)-yi(j))^2) < d_min, keep(j) = false; end
        end
    end
    xo = xi(keep); yo = yi(keep);
end

function s_out = mergeCoincidentSpots(s, d_min)
    keep = true(size(s,1), 1);
    for i = 1:size(s,1)
        if ~keep(i), continue; end
        for j = i+1:size(s,1)
            if keep(j) && sqrt(sum((s(i,1:2)-s(j,1:2)).^2)) < d_min
                if s(j,5) > s(i,5), keep(i) = false; else, keep(j) = false; end
            end
        end
    end
    s_out = s(keep, :);
end