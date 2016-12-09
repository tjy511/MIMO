%% Calculate and load array factor / beampattern
%
% TJ Young (adapted from arrayFactor scripts)
% 05 December 2016

%% Create 2D array

% Obtain parameters
%vdat = fmcw_load('DATA2014-08-03-1755.DAT', 1);

% Calculate virtual antenna locations
array2d_final_gland
xPos = ve(:,1); yPos = ve(:,2);
clearvars -except xPos yPos dphy

% Calculate weights
w = ones(1,numel(xPos))/numel(xPos);

%% Recreate bowtie antenna radiation pattern
% Temporarily eyeballed radiation pattern from graph, waiting on Lai Bun to
% send over data...
clear bowtie
bowtieRaw(1,:) = deg2rad(0:10:360); % theta
bowtieRaw(2,:) = [0 -0.5 -1.5 -3 -5 -7.5 -11 -15 -25 ...
    -50 -20 -16 -15 -16 -18 -22 -14 -10 ...
    -9 -10 -14 -22 -18 -16 -15 -16 -20 ...
    -50 -25 -15 -11 -7.5 -5 -3 -1.5 -0.5 0]; % Azimuth
bowtieRaw(3,:) = [0 -0.01 -0.025 -0.1 -0.25 -0.5 -1 -2 -3.25 -5 ...
    -6.75 -8.75 -11 -12.75 -14 -15 -12 -10 -9 ...
    -10 -12 -15 -14 -12.75 -11 -8.75 -6.75 ...
    -5 -3.25 -2 -1 -0.5 -0.25 -0.1 -0.025 -0.01 0]; % Elevation
bowtieRaw(2,:) = db2mag(bowtieRaw(2,:)); bowtieRaw(3,:) = db2mag(bowtieRaw(3,:)); 
bowtieInterp(1,:) = deg2rad(0:1:360); % theta
bowtieInterp(2,:) = interp1(bowtieRaw(1,:),bowtieRaw(2,:),bowtieInterp(1,:),'pchip');
bowtieInterp(3,:) = interp1(bowtieRaw(1,:),bowtieRaw(3,:),bowtieInterp(1,:),'pchip');

% Plotting fancies
figure
pax = polaraxes; subplot(2,1,1,pax);
polarplot(bowtieInterp(1,:),db(bowtieInterp(2,:)));
pax.RLim = [-50 0]; pax.ThetaZeroLocation = 'bottom';
title('Azimuth')
pax = polaraxes; subplot(2,1,2,pax);
polarplot(bowtieInterp(1,:),db(bowtieInterp(3,:)));
pax.RLim = [-20 0]; pax.ThetaZeroLocation = 'bottom';
title('Elevation')

% Radiation parameter computation (this doesn't work, fix this)
ant_hpbw(1) = hpbw(dB(bowtieInterp(2,:)),bowtieInterp(1,:)); % Azimuth HPBW [deg]
ant_hpbw(2) = hpbw(dB(bowtieInterp(3,:)),bowtieInterp(1,:)); % Elevation HPBW [deg]

%% Plot array geometry and array factor for different frequencies

% Wave-frequency and wave-speed
freqs = [2e8 3e8 4e8];
%f = vdat.f;
c = 3e8/sqrt(3.18); % Eq. 6.2 of Hubbard & Glasser 2005

% Scanning angles
thetaScA = -90:90;
phiScA = 0;

% Steering angle
thetaStA = 0; 

% Calculate array pattern and compute radiation pattern 
bowtieWeight = [bowtieInterp(3,271:360) bowtieInterp(3,1:91)];

for ii = 1:numel(freqs)
    freq = freqs(ii);
    W(ii,:) = arrayFactor(xPos,yPos,w,freq,c,thetaScA,phiScA,thetaStA); % Unweighted
    WWeight(ii,:) = W(ii,:) .* bowtieWeight; % Weight W with antenna beamform
end

% Calculate HPBW
W_hpbw = hpbw(dB(WWeight),thetaScA);

% Plot the array pattern for various frequencies 
%plotBeampattern(db(W), w, freqs, c, thetaScA,thetaStA) % Unweighted
plotBeampattern(db(WWeight), w, freqs, c, thetaScA,thetaStA) % Weighted

%% Plot the 3D polar beampattern of the array with plotBeampattern3D()

plotBeampattern3D(xPos, yPos, w)
