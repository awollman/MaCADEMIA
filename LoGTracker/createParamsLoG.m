% STOICH_PARAMS initializes parameters for stoichAnalyser and generalLinearStoich.
%% Stoichiometry Calculation
params.iSingle = [];              % Single molecule intensity (empty = auto-calculate)
params.lateFrameFraction = 0.8;   % Movie fraction used to find single molecules
params.nFrames = 5;               % Frames used for initial intensity extrapolation
params.minFrames = 5;             % Minimum trajectory length required for analysis
params.threshold = 150;           % Maximum stoichiometry allowed

%% Visualization
params.showPlot = 1;              % Enable or disable KDE distribution plot
params.kdeBandwidth = 5;          % Smoothing factor for the KDE curve