% Parameters for filter
cfg.filter = 1; % Switch to filter/smooth profile
cfg.ftype = 'gaussian'; % Rotationally symmetric Gaussian lowpass filter
cfg.fsize = 7; % Filter size
cfg.fsigma = 1; % Standard deviation threshold
cfg.fwindow = 4; % Window size in convolution

% Parameters for peak identification
cfg.pkselect = 1; % Use selected peaks from 2D processing
cfg.pktolm = 5; % 2D tolerance for peaks (bins)
cfg.pkthresh = -50; % dB threshold level for peaks
cfg.pkprom = 0; % Filter by prominence threshold
cfg.phiLim = [225 315]; % Limits for phi (degrees)
cfg.rThresh = 5; % Search range for pkselect; 

% Parameters for plotting
cfg.doPlot = 0; % Turn on intermediate plotting

% Activate config
if cfg.pkselect == 1
    load('intSelect.mat')
end
if cfg.doPlot == 0
    set(0,'DefaultFigureVisible','off')
end