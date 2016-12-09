

%% Calculate and load array factor / beampattern
%
% TJ Young (adapted from arrayFactor scripts)
% 05 December 2016

%% Load bowtie antenna radiation pattern
[ant.theta,ant.RA,ant.RE] = antennaBP('bowtie',1);

% Radiation parameter computation (this doesn't work, fix this)
ant.hpbwRA(1) = hpbw(dB(ant.RA),ant.theta); % Azimuth HPBW [deg]
ant.hpbwRE(2) = hpbw(dB(ant.RE),ant.theta); % Elevation HPBW [deg]

%% Load 2D array

% Calculate virtual antenna locations
pos = array2d_final_gland;
xPos = pos(:,1); yPos = pos(:,2);

% Calculate weights
w = ones(1,numel(xPos))/numel(xPos);

%% Calculate array geometry and array factor for different frequencies

% Wave-frequency and wave-speed
freqs = [2e8 3e8 4e8];
c = 3e8/sqrt(3.18); % Eq. 6.2 of Hubbard & Glasser 2005

% Scanning angles
thetaScA = ant.theta;
phiScA = 0;

% Steering angle
thetaStA = 0; 

% Calculate array pattern and compute radiation pattern 
%bowtieWeight = [bowtieInterp(3,271:360) bowtieInterp(3,1:91)];

for ii = 1:numel(freqs)
    freq = freqs(ii);
    W(ii,:) = arrayFactor(xPos,yPos,w,freq,c,thetaScA,phiScA,thetaStA); % Unweighted
    WWeight(ii,:) = W(ii,:) .* ant.RE; % Weight W with antenna beamform
end

% Calculate HPBW
W_hpbw = hpbw(dB(WWeight),thetaScA);
%%
% Plot the array pattern for various frequencies 
plotBeamPattern(db(W(:,901:2701)), freqs, thetaScA(:,901:2701)) % Unweighted
plotBeamPattern(db(WWeight(:,901:2701)), freqs, thetaScA(:,901:2701)) % Weighted
