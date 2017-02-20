%function beamPatternDphy(antType,txrx,dPhy)

%% Calculate and load array factor / beampattern
%
% TJ Young (adapted from arrayFactor scripts)
% 05 December 2016

%% Input parameters

% Antenna parameters
%ant.type = antType;
ant.type = 'pencil'; % 'bowtie' 'helix' 'pencil' 'isotropic' 'dipole'
thetaStA = 0; % Steering angle
freqs = 3e8; % Wave frequency [Hz]
c = 3e8; %/sqrt(3.18); % Wave speed [m/s] (Eq. 6.2 of Hubbard & Glasser 2005)

% Antenna locations
ant.loc = 'virtual'; % 'store1' 'other'
ant.wid = 0.00333;
txrx = [16 1]; % Default: [8 8]
quadrant = 4; % Default: 4
off = [2.49 -2.49]; % Default: [2.49 -2.49]
dPhy = 1; % Default: 0.83 % NOTE that dPhy is half of the resulting element separation

doPlot = 1; 
doHPBW = 1; 
doSave = 0;

%% Calculate array geometry and array factor for different frequencies
for ii = 1:length(ant.wid)
    
    % Load bowtie antenna radiation pattern
    [ant.theta,ant.RA(ii,:),ant.RE(ii,:),fig_bp(ii)] = antennaBP(ant.type,doPlot,ant.wid(ii));
    
    % Scanning angles
    thetaScA = ant.theta;
    phiScA = 0;
    
    % Calculate virtual antenna locations
    switch ant.loc
        case 'store1'
            [ve,fig_ant] = antennaLoc(ant.loc,doPlot);
        case 'real'
            [ve,fig_ant] = antennaLoc(ant.loc,txrx,dPhy,quadrant,off,doPlot);
        case 'virtual' 
            [ve,fig_ant] = antennaLoc(ant.loc,txrx,dPhy,doPlot);
    end

    xPos = ve(:,1); yPos = ve(:,2);
    freq = freqs;
    w = ones(1,numel(xPos))/numel(xPos);
    
    W(ii,:) = arrayFactor(xPos,yPos,w,freq,c,thetaScA,phiScA,thetaStA); % Unweighted
    WWeight(ii,:) = W(ii,:) .* ant.RE(ii,:); % Weight W with antenna beamform
end

%% Plot the array pattern for various frequencies 
WComb = [ant.RE; W(1,:); WWeight]; 
[fig_w,ax_w,pp] = plotBPRectangular(db(WComb(:,901:2701),'power'), thetaScA(:,901:2701), 50,[nan ant.wid]);%, txrx); % Unweighted
fig_w = gcf; 
set(fig_w,'units','normalized','position',[0 0 1 1/2])

set(pp(1),'color','k','lineStyle','--')
set(pp(2),'color', [51 153 255]/255,'lineStyle',':')
set(pp(3),'color', [0 76 153]/255,'lineWidth',2)

%% Calculate HPBW

if doHPBW
    for ii = 1:length(ant.wid)
    ant.hpbwRA(1) = hpbw(ant.theta,ant.RA(ii,:)); % Azimuth HPBW [deg]
    ant.hpbwRE(2) = hpbw(ant.theta,ant.RE(ii,:)); % Elevation HPBW [deg]
    %W_hpbw = hpbw(thetaScA,W); % Array factor HPBW [deg]
    WWeight_hpbw = hpbw(thetaScA,WWeight(ii,:)); % Array factor HPBW [deg]
    end
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

