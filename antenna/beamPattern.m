%% Calculate and load array factor / beampattern
%
% TJ Young (adapted from arrayFactor scripts)
% 05 December 2016

ccc

%% Input parameters

% Antenna parameters
ant.type = 'isotropic'; % 'bowtie' 'helix' 'pencil' 'isotropic'
thetaStA = 30; % Steering angle
freqs = [2e8 3e8 4e8]; % Wave frequency [Hz]
c = 3e8/sqrt(3.18); % Wave speed [m/s] (Eq. 6.2 of Hubbard & Glasser 2005)

% Antenna locations
ant.loc = 'other'; % 'store1' 'other'
txrx = [32 1]; % Default: [8 8]
quadrant = 4; % Default: 4
off = [2.49 -2.49]; % Default: [2.49 -2.49]
dPhy = 0.83; % Default: 0.83

doPlot = 1; 
doHPBW = 0; 
doSave = 1;

%% Load bowtie antenna radiation pattern
[ant.theta,ant.RA,ant.RE,fig_bp] = antennaBP(ant.type,doPlot);

%% Load 2D array

% Calculate virtual antenna locations
if strcmp(ant.loc,'store1') == 1
    [ve,fig_ant] = antennaLoc(ant.loc,doPlot);
else
    [ve,fig_ant] = antennaLoc(txrx,quadrant,off,dPhy,doPlot);
end
xPos = ve(:,1); yPos = ve(:,2);

% Calculate weights
w = ones(1,numel(xPos))/numel(xPos);

%% Calculate array geometry and array factor for different frequencies

% Scanning angles
thetaScA = ant.theta;
phiScA = 0;

% Calculate array pattern and compute radiation pattern 
%bowtieWeight = [bowtieInterp(3,271:360) bowtieInterp(3,1:91)];

for ii = 1:numel(freqs)
    freq = freqs(ii);
    W(ii,:) = arrayFactor(xPos,yPos,w,freq,c,thetaScA,phiScA,thetaStA); % Unweighted
    WWeight(ii,:) = W(ii,:) .* ant.RE; % Weight W with antenna beamform
end

%% Calculate HPBW
if doHPBW
    ant.hpbwRA(1) = hpbw(dB(ant.RA),ant.theta); % Azimuth HPBW [deg]
    ant.hpbwRE(2) = hpbw(dB(ant.RE),ant.theta); % Elevation HPBW [deg]
    W_hpbw = hpbw(dB(WWeight),thetaScA); % Array factor HPBW [deg]
end

%% Plot the array pattern for various frequencies 
fig_uw = plotBeamPattern(db(W(:,901:2701)), freqs, thetaScA(:,901:2701)); % Unweighted
title(['Array factor for unweighted phased array of ',ant.type,' antennas'])
set(gcf,'units','normalized','position',[0 0 1/2 1])
fig_w = plotBeamPattern(db(WWeight(:,901:2701)), freqs, thetaScA(:,901:2701)); % Weighted
title(['Array factor for weighted phased array of ',ant.type,' antennas'])
set(gcf,'units','normalized','position',[1/2 0 1/2 1])

%% Saving fancies

if doSave
    % Create and cd to folder
    fileLoc = strcat('~/Google Drive/Academic/papers/paper3/figs/recreation/',...
        num2str(txrx(1)),'-',num2str(txrx(2)),'_',num2str(dPhy*1000));
    try
        cd(fileLoc);
    catch
        mkdir(fileLoc); cd(fileLoc);
    end
    set([fig_w fig_uw],'color','w')
    export_fig(fig_uw,strcat(ant.type,'_AF_uw_',num2str(thetaStA),'.png'),'-m2');
    export_fig(fig_w,strcat(ant.type,'_AF_w_',num2str(thetaStA),'.png'),'-m2');
    if doPlot
        set([fig_bp fig_ant],'color','w')
        export_fig(fig_bp,strcat(ant.type,'_RP','.png'),'-m2')
        export_fig(fig_ant,strcat('arraySetup_',num2str(txrx(1)),'-',...
            num2str(txrx(2)),'_',num2str(dPhy*1000),'.png'),'-m2');
    end
end

