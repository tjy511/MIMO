function cfg = mimo_config_field
% Config for profile recreations

% Default config
cfg.notes = 'Setup for field recreation';

% Layer profiles
cfg.noiseFloor = -100; % Noise floor level [dBV]
cfg.intLevel = -40; % Internal reflector reflectivity [dBV]
cfg.bedLevel = -20; % Bed reflector reflectivity [dBV]
cfg.bedDepth = 600; % Depth of bed reflector [m]
cfg.bedRoughness = [pi 10]; % Bed roughness; [amplitude period]
cfg.bedMountain = []; % Subglacial bump; [centre-location width height]

% Backscatter
cfg.scatterType = 'gaussian'; % 'gaussian' 'cosine'
cfg.sigmai = 0.25; % Specularity factor for internal layer (standard bell curve = 1) [degrees]
cfg.sigmab = 1; % Specularity factor for bed layer (standard bell curve = 1) [degrees]

% Antenna location
cfg.antSetup = 'store1'; % 'store1' 'real' 'virtual'
cfg.txrx = [8 8]; % Row/Column of antenna array; Default: [8 8]
cfg.quadrant = 2; % Default: 2
cfg.off = [2.49 -2.49]; % Offset between Tx1 / Rx1; Default: [2.49 -2.49] [m]
cfg.beamwidth = 1; % Variable for beam width (default = 1)
cfg.dPhy = 0.83; % Antenna separation; Default: 0.83 [m]

% Beam forming
cfg.freq = 3e8; % MIMO cfg.frequency [Hz]
cfg.thetaStA = -30:30; % Steering angle range [degrees]
cfg.antSelect = 'bowtie'; % Antenna selection ('isotropic' 'pencil' 'bowtie' 'helix' 'dipole')
cfg.c = 3e8/sqrt(3.1); % Speed of wave propagation (3e8/sqrt(3.18) in ice);
cfg.weight = 1; % Weighting beam with array factor [binary]

% Saving
cfg.doSave = 1; % Save output
cfg.fileLoc = '~/Google Drive/Academic/papers/paper3/figs/';

% Switches
cfg.doConvolution = 1; % Do 2D convolution smoothing of profile
cfg.doPlot = 1; % Plot intermediate steps