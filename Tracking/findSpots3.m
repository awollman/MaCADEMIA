function [result] = findSpots3(frame, method, disk_radius, gaussian, output)
% FINDSPOTS3 finds candidate bright spots using histogram-based thresholding.

    if output == 1, pause on; end

    % Convert to grayscale and optionally filter
    frame0 = mat2gray(frame);
    if gaussian == 1
        frame0 = imfilter(frame0, fspecial('gaussian', 3));
    end
    
    show_step(frame0, 'Initial processed frame', output);

    % Top-hat transformation to flatten background
    se = strel('disk', disk_radius);
    signal = imtophat(frame0, se);
    
    show_step(signal, 'Top-hat enhanced signal', output);

    % Calculate threshold from histogram peak using FWHM
    [hist_y, hist_x] = imhist(signal(signal > 0));
    [width_x, max_val] = fwhm(hist_x, hist_y);
    
    if width_x > 0
        auto_thresh = max_val;
    else
        auto_thresh = hist_x(hist_y == max(hist_y));
    end

    % Filter and binarise
    gauss_signal = imfilter(signal, fspecial('gaussian'));
    sig_bw = imbinarize(gauss_signal, auto_thresh(1));
    
    show_step(gauss_signal, 'Gaussian filtered signal', output);
    show_step(sig_bw, 'Binarised signal', output);

    % Morphological cleanup
    b_clean = imopen(sig_bw, strel('disk', 1));
    b_clean = bwmorph(b_clean, 'fill', 1);
    
    show_step(b_clean, 'Cleaned binary image', output);

    % Extract candidates based on method
    if method == 1
        % Successive erosions and shrinking
        b_erode = imerode(b_clean, ones(2,2));
        b_erode = imerode(b_erode, ones(2,2));
        result = bwmorph(b_erode, 'shrink', 1);
    else
        % Ultimate erosion
        if max(b_clean(:)) > 0
            result = bwulterode(b_clean);
        else
            result = b_clean;
        end
    end

    % Check for CellMask in scope
    if exist('CellMask', 'var') && exist('SignalMask', 'var')
        result = result .* SignalMask;
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