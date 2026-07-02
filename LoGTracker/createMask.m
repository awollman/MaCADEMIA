function mask = createMask(filePath)
% CREATEMASK converts a TIFF image into a binary mask for stoichiometry analysis.

    % Load the image data from the specified TIFF file
    rawMask = imread(filePath);

    % Convert the image to a logical array where non-zero pixels are true
    % This matches the stoichAnalyser requirement: mask(y, x) > 0
    mask = rawMask > 0;

    % Ensure the mask is 2D by selecting the first plane if the file is multi-channel or RGB
    if size(mask, 3) > 1
        mask = mask(:, :, 1);
    end
end