%% Calculate and load array factor / beampattern
%
% TJ Young (adapted from arrayFactor scripts)
% 05 December 2016

ccc

%% Input parameters

% Antenna parameters
ant.type = 'isotropic'; % 'bowtie' 'helix' 'pencil' 'isotropic' 'dipole'
thetaStA = 0; % Steering angle
freqs = 3e8; % Wave frequency [Hz]
c = 3e8;%/sqrt(3.18); % Wave speed [m/s] (Eq. 6.2 of Hubbard & Glasser 2005)

% Antenna locations
ant.loc = 'virtual'; % 'store1' 'real' 'virtual'
txrx = [32 1]; % Default: [8 8]
quadrant = 4; % Default: 4
off = [2.49 -2.49]; % Default: [2.49 -2.49]
dPhy = 1 * c/3e8; %3 * c/3e8; % Default: 0.83
ampWeight = 'binomial'; % 'uniform' 'triangular' 'binomial'

doPlot = 1;
doHPBW = 1;
doSave = 0;

%% Load bowtie antenna radiation pattern
[ant.theta,ant.RA,ant.RE,fig_bp] = antennaBP(ant.type,doPlot);

%% Calculate amplitude weighting
lowpt = 0.25; highpt = 1;
if mod(txrx(1),2) == 0
    midpt = [txrx(1)/2 txrx(1)/2+1];
else
    midpt = txrx(1)/2 + 0.5;
end

switch ampWeight
    
    case 'uniform'
        w = ones(1,txrx(1)); % Weights
    case 'triangular'
        if mod(txrx(1),2) == 0
            w = interp1([1 midpt txrx(1)],[lowpt,highpt,highpt,lowpt],1:txrx(1));
        else
            w = interp1([1 midpt txrx(1)],[lowpt,highpt,lowpt],1:txrx(1));
        end
        w = w ./ sum(w);
    case 'binomial'
        w = binopdf(1:txrx(1),txrx(1),0.5);
        if mod(txrx(1),2) == 0
            w = [w(1:midpt(1)) w(midpt(1)) w(midpt(1)+1:end-1)];
        end
end

w = w ./ sum(w);

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
    
    xPos = ve(:,1); yPos = ve(:,2); % Virtual element positions
    freq = freqs;
    
    W(ii,:) = arrayFactor(xPos,yPos,w,freq,c,thetaScA,phiScA,thetaStA); % Unweighted
    WWeight(ii,:) = W(ii,:) .* ant.RE; % Weight W with antenna beamform
end

%% Plot the array pattern for various frequencies
fig_uw = plotBeamPattern(db(W(:,901:2701),'power'), freqs, thetaScA(:,901:2701)); % Unweighted
set(gcf,'units','normalized','position',[0 0 1/2 1])

fig_w = plotBeamPattern(db(WWeight(:,901:2701),'power'), freqs, thetaScA(:,901:2701)); % Weighted
set(gcf,'units','normalized','position',[1/2 0 1/2 1])

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
    fileLoc = strcat('~/Google Drive/Academic/papers/paper3/figs/wavelength/',...
        num2str(txrx(1)),'-',num2str(txrx(2)),'_',num2str(dPhy*1000));
    try
        cd(fileLoc);
    catch
        mkdir(fileLoc); cd(fileLoc);
    end
    set([fig_w fig_uw],'color','w')
    %export_fig(fig_uw,strcat(ant.type,'_AF_uw_',num2str(thetaStA),'.png'),'-m2');
    export_fig(fig_w,strcat(ant.type,'_AF_w_',num2str(thetaStA),'.png'),'-m2');
    if doPlot
        set([fig_bp fig_ant],'color','w')
        export_fig(fig_bp,strcat(ant.type,'_RP','.png'),'-m2')
        export_fig(fig_ant,strcat('arraySetup_',num2str(txrx(1)),'-',...
            num2str(txrx(2)),'_',num2str(dPhy*1000),'.png'),'-m2');
    end
end

