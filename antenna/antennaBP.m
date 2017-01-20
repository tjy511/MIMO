function varargout = antennaBP(antType,plotType,a)
% Calculate antenna array factor / beam pattern
%
% TJ Young
% 08 December 2016

theta =-180:0.1:180; % theta

if ~exist('a','var')
    a = 1;
end

% Calculate antenna beam pattern
switch antType
    
    % Isotropic reflector
    case 'isotropic'
    RA = ones(size(theta));
    RE = RA; 
    
    % Synthesised pencil beam
    case 'pencil'
        % Scanning angles
        RA = -a.*(theta).^2;
        RA = db2pow(RA);
        RE = RA;
        
    case 'dipole'
        % Currently similar to l = 0.1 * wavelength, perhaps too wide. 
        RA = ones(size(theta));
        RE1 = -1/(90^2) * theta(901:2700).^2 + 1;
        RE2 = fliplr(RE1); 
        RE = [RE1 0 RE2];
        RE = circshift(RE,[0 901]);
        
    % Bowtie antennas used on Store
    case 'bowtie'   
        % NOTE: Temporarily eyeballed radiation pattern from graph
        % , waiting on Lai Bun to send over data...
        %thetaRaw = 0:10:360; % theta
        thetaRaw = -180:10:180; % theta
        RARaw = [-9 -10 -14 -22 -18 -16 -15 -16 -20 ...
            -50 -25 -15 -11 -7.5 -5 -3 -1.5 -0.5 ...      
            0 -0.5 -1.5 -3 -5 -7.5 -11 -15 -25 ...
            -50 -20 -16 -15 -16 -18 -22 -14 -10 -9]; % Azimuth
        RERaw = [-9 -10 -12 -15 -14 -12.75 -11 -8.75 -6.75 ...
            -5 -3.25 -2 -1 -0.5 -0.25 -0.1 -0.025 -0.01 ...
            0 -0.01 -0.025 -0.1 -0.25 -0.5 -1 -2 -3.25 -5 ...
            -6.75 -8.75 -11 -12.75 -14 -15 -12 -10 -9]; % Elevation
        RARaw = db2pow(RARaw); RERaw = db2pow(RERaw);
        
        RA = interp1(thetaRaw,RARaw,theta,'pchip');
        RE = interp1(thetaRaw,RERaw,theta,'pchip');
    
    % Theoretical pattern for helical antennas    
    case 'helix'
        % Radiation pattern from Kraus (1947) / Kraus and Williamson (1948)
        % 12 degrees / 7-turn helix at 450 Mc
        thetaRaw = -180:10:180; % theta
        RARaw = [-100 -100 -100 -100 -100 -100 -100 -100 -100 ...
            -100 -50 -41 -36 -41 -50 -26 -14 -5 ...
            0 -5 -14 -26 -50 -41 -36 -41 -50 -100 ...
            -100 -100 -100 -100 -100 -100 -100 -100 -100]; 
        RERaw = RARaw;
        RARaw = db2pow(RARaw); RERaw = db2pow(RERaw);
        
        RA = interp1(thetaRaw,RARaw,theta,'pchip');
        RE = interp1(thetaRaw,RERaw,theta,'pchip');

    case 'helix2'
        % Helical antenna with a narrower beam than helix1
        thetaRaw = -180:370/41:180; % theta
        thetaRaw = [thetaRaw 180];
        RARaw = [-100 -100 -100 -100 -100 -100 -100 -100 -100 -100 -100 ...
            -100 -50 -41 -36 -41 -50 -26 -14 -5 ...
            0 -5 -14 -26 -50 -41 -36 -41 -50 -100 ...
            -100 -100 -100 -100 -100 -100 -100 -100 -100 -100 -100]; 
        RERaw = RARaw;
        RARaw = db2pow(RARaw); RERaw = db2pow(RERaw);
        
        RA = interp1(thetaRaw,RARaw,theta,'pchip');
        RE = interp1(thetaRaw,RERaw,theta,'pchip');

    case 'helix3'
        % Helical antenna with the narrowest beam
        thetaRaw = -180:370/45:180; % theta
        thetaRaw = [thetaRaw 180];
        RARaw = [-100 -100 -100 -100 -100 -100 -100 -100 -100 -100 -100 -100 -100 ...
            -100 -50 -41 -36 -41 -50 -26 -14 -5 ...
            0 -5 -14 -26 -50 -41 -36 -41 -50 -100 ...
            -100 -100 -100 -100 -100 -100 -100 -100 -100 -100 -100 -100 -100]; 
        RERaw = RARaw;
        RARaw = db2pow(RARaw); RERaw = db2pow(RERaw);
        
        RA = interp1(thetaRaw,RARaw,theta,'pchip');
        RE = interp1(thetaRaw,RERaw,theta,'pchip');

    case 'helix4'
        % Helical antenna with the narrowest beam
        thetaRaw = -180:360/32:180; % theta
        RARaw = [-100 -100 -100 -100 -100 -100 -100 ...
            -100 -50 -41 -36 -41 -50 -26 -14 -5 ...
            0 -5 -14 -26 -50 -41 -36 -41 -50 -100 ...
            -100 -100 -100 -100 -100 -100 -100]; 
        RERaw = RARaw;
        RARaw = db2pow(RARaw); RERaw = db2pow(RERaw);
        
        RA = interp1(thetaRaw,RARaw,theta,'pchip');
        RE = interp1(thetaRaw,RERaw,theta,'pchip');

    case 'helix5'
        % Helical antenna with the narrowest beam
        thetaRaw = -180:360/28:180; % theta
        RARaw = [-100 -100 -100 -100 -100 ...
            -100 -50 -41 -36 -41 -50 -26 -14 -5 ...
            0 -5 -14 -26 -50 -41 -36 -41 -50 -100 ...
            -100 -100 -100 -100 -100]; 
        RERaw = RARaw;
        RARaw = db2pow(RARaw); RERaw = db2pow(RERaw);
        
        RA = interp1(thetaRaw,RARaw,theta,'pchip');
        RE = interp1(thetaRaw,RERaw,theta,'pchip');
end

% Plotting fancies
if plotType == 0
    set(0,'DefaultFigureVisible','off')
end

fig = figure;

pax = polaraxes; subplot(2,1,1,pax);
polarplot(deg2rad(theta(db(RA,'power')>-50)),db(RA(db(RA,'power')>-50),'power'));
pax.RLim = [-50 0]; pax.ThetaZeroLocation = 'bottom';
title('Azimuth')

pax = polaraxes; subplot(2,1,2,pax);
polarplot(deg2rad(theta(db(RE,'power')>-50)),db(RE(db(RE,'power')>-50),'power'));
pax.RLim = [-50 0]; pax.ThetaZeroLocation = 'bottom';
title('Elevation')

% Export variables
if nargout == 3
    varargout{1} = theta;
    varargout{2} = RA;
    varargout{3} = RE;
elseif nargout == 4
    varargout{1} = theta;
    varargout{2} = RA;
    varargout{3} = RE;
    varargout{4} = fig;
end

set(0,'DefaultFigureVisible','on')