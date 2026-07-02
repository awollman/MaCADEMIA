% CREATEPLOG initializes parameters for trackerLoG.

%% Image Loading and Frame Selection
p.useBioFormats = 1;         % Use BioFormats for non-TIFF compatibility
p.all = 1;                   % Load the entire image file
p.startFrame = 1;            % Starting frame if 'all' is 0
p.endFrame = 999;            % Ending frame if 'all' is 0
p.DetermineFirstFrames = 1;  % Auto-detect frame where laser turns on
p.firstFrame = 1;            % Manual start frame if auto-detect is 0
p.FramesToTrack = 0;         % Frames to track (0 = track until end of movie)
p.noFrames = 5;              % Number of frames to average for background
p.satPixelVal = 10^10;       % Pixel value threshold for saturation

%% Detection (Finding Candidates)
p.DetectionMethod = 'LoG';   % Method: 'LoG' (recommended) or 'original'
p.disk_radius = 5;           % Size of morphological element for spot finding
p.gaussian = 0;              % Apply Gaussian filter before detection
p.Candidate_d_min = 0;       % Minimum distance between candidate spots

%% Localisation (Refining Positions)
p.LocalisationMethod = 'iterative'; % Method: 'iterative' or 'quadratic'
p.subarray_halfwidth = 8;           % Half-width of the processing window
p.inner_circle_radius = 5;          % Radius for signal integration
p.gauss_mask_sigma = 2;             % Sigma for the Gaussian tracking mask
p.error_set = 0.05;                 % Convergence threshold for iterative masking
p.SNR_min = 0.4;                    % Minimum signal-to-noise ratio to accept a spot
p.sigmaFit_min = 0;                 % Minimum acceptable spot width
p.sigmaFit_max = 5;                 % Maximum acceptable spot width

%% Linking (Building Trajectories)
p.LinkingMethod = 'LAP';     % Method: 'LAP' (recommended) or 'greedy'
p.d_01_max = 5;              % Maximum distance for frame-to-frame linking
p.gap_frames = 3;            % Maximum frames a spot can disappear (LAP only)
p.Iratio_01_min = 0.5;       % Minimum intensity ratio between linked frames
p.Iratio_01_max = 3;         % Maximum intensity ratio between linked frames
p.d_min = 1;                 % Distance threshold to merge coincident spots

%% Output and Display
p.show_output = 0;           % Show diagnostic plots during tracking
p.show_text_output = 1;      % Print progress to the command window