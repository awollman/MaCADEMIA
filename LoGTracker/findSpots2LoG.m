function [result] = findSpots2LoG(frame, method, disk_radius, gaussian, output)
% FINDSPOTS2LOG finds candidate bright spots using morphological background estimation.

    if output == 1, pause on; end

    % Convert to grayscale and optionally filter
    frame0 = mat2gray(frame);
    if gaussian == 1
        frame0 = imfilter(frame0, fspecial('gaussian', 3));
    end
    
    show_step(frame0, 'Initial processed frame', output);

    % Background estimation using morphological opening
    se = strel('disk', disk_radius);
    signal = imtophat(frame0, se);
    bgnd = imopen(frame0, se);
    
    show_step(signal, 'Top-hat enhanced signal', output);
    show_step(bgnd, 'Background estimate', output);

    % Create signal mask based on background threshold
    signal_mask = imbinarize(bgnd, graythresh(bgnd));
    signal_only = signal .* signal_mask;
    
    show_step(signal_mask, 'Background threshold mask', output);
    show_step(signal_only, 'Masked signal', output);

    % Binarise and clean up signal
    sig_bw = imbinarize(signal_only, graythresh(signal_only));
    b_clean = imopen(sig_bw, strel('disk', 1));
    b_clean = bwmorph(b_clean, 'fill', 1);
    b_clean = bwmorph(b_clean, 'open', 1);
    
    show_step(b_clean, 'Cleaned binary image', output);

    % Extract candidates based on method
    if method == 1
        % Shrink objects to points via erosion
        b_erode = imerode(b_clean, ones(2,2));
        b_erode = imerode(b_erode, ones(2,2));
        result = bwmorph(b_erode, 'shrink', 1);
    else
        % Ultimate erosion to find regional maxima
        if max(b_clean(:)) > 0
            result = bwulterode(b_clean);
        else
            result = b_clean;
        end
    end

    % Apply final mask if defined in scope
    if exist('CellMask', 'var')
        result = result .* signal_mask;
    end

    show_step(result, 'Final candidate spots', output);
end

function show_step(img, msg, output)
    % Helper to display images if output is enabled
    if output == 1
        imshow(img, []);
        title(msg);
        pause;
    end
end