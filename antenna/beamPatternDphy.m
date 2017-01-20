%function beamPatternDphy(antType,txrx,dPhy)

%% Calculate and load array factor / beampattern
%
% TJ Young (adapted from arrayFactor scripts)
% 05 December 2016

%% Input parameters

% Antenna parameters
%ant.type = antType;
ant.type = 'isotropic'; % 'bowtie' 'helix' 'pencil' 'isotropic' 'dipole'
thetaStA = 0; % Steering angle
freqs = 3e8; % Wave frequency [Hz]
c = 3e8; %/sqrt(3.18); % Wave speed [m/s] (Eq. 6.2 of Hubbard & Glasser 2005)

% Antenna locations
ant.loc = 'virtual'; % 'store1' 'other'
txrx = [16 1]; % Default: [8 8]
quadrant = 4; % Default: 4
off = [2.49 -2.49]; % Default: [2.49 -2.49]
dPhy = [2 1 1/2]; % Default: 0.83 % NOTE that dPhy is half of the resulting element separation

doPlot = 1; 
doHPBW = 1; 
doSave = 1;

%% Load bowtie antenna radiation pattern
[ant.theta,ant.RA,ant.RE,fig_bp] = antennaBP(ant.type,doPlot);

%% Calculate array geometry and array factor for different frequencies

% Scanning angles
thetaScA = ant.theta;
phiScA = 0;

for ii = 1:numel(dPhy)
    
    % Calculate virtual antenna locations
    switch ant.loc
        case 'store1'
            [ve,fig_ant] = antennaLoc(ant.loc,doPlot);
        case 'real'
            [ve,fig_ant] = antennaLoc(ant.loc,txrx,dPhy(ii),quadrant,off,doPlot);
        case 'virtual' 
            [ve,fig_ant] = antennaLoc(ant.loc,txrx,dPhy(ii),doPlot);
    end

    xPos = ve(:,1); yPos = ve(:,2);
    freq = freqs;
    w = ones(1,numel(xPos))/numel(xPos);
    
    W(ii,:) = arrayFactor(xPos,yPos,w,freq,c,thetaScA,phiScA,thetaStA); % Unweighted
    WWeight(ii,:) = W(ii,:) .* ant.RE; % Weight W with antenna beamform
end

%% Plot the array pattern for various frequencies 

[fig_w,ax_w] = plotBPPolar(db(WWeight(:,901:2701),'power'), thetaScA(:,901:2701), 50, dPhy); % Unweighted
set(gcf,'units','normalized','position',[0 0 1/2 1/2])

%% Calculate HPBW 

if doHPBW
    ant.hpbwRA(1) = hpbw(ant.theta,ant.RA); % Azimuth HPBW [deg]
    ant.hpbwRE(2) = hpbw(ant.theta,ant.RE); % Elevation HPBW [deg]
    W_hpbw = hpbw(thetaScA,W); % Array factor HPBW [deg]
    WWeight_hpbw = hpbw(thetaScA,WWeight); % Array factor HPBW [deg]
end

%% Saving fancies

if doSave
    % Create and cd to folder
    %fileLoc = '~/Google Drive/Academic/papers/paper3/figs/paper/';
    fileLoc = '~/Downloads/';
    cd(fileLoc)
    set(fig_w,'color','w')
    export_fig(fig_w,strcat(ant.type,'_AF_w_',num2str(thetaStA),'_',...
        num2str(txrx(1)),'-',num2str(txrx(2)),'_',num2str(dPhy*1000),'.png'),'-m2');
end

