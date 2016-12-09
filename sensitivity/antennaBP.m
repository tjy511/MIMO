function [theta,RA,RE] = antennaBP(antType,plotType,a)
% Calculate antenna array factor / beam pattern
%
% TJ Young
% 08 December 2016

theta =-180:0.1:180; % theta

if ~exist('a','var')
    a = 0;
end

% Calculate antenna beam pattern
switch antType
    
    % Isotropic reflector
    case 'isotropic'
    RA = ones(1,length(theta));
    RE = RA; 
    
    % Synthesised pencil beam
    case 'pencil'
        % Scanning angles
        RA = -(theta-a).^2;
        RA = db2mag(RA);
        RE = RA;
        
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
        RARaw = db2mag(RARaw); RERaw = db2mag(RERaw);
        
        RA = interp1(thetaRaw,RARaw,theta,'pchip');
        RE = interp1(thetaRaw,RERaw,theta,'pchip');
end

% Plotting fancies
if plotType == 1
    figure
    
    pax = polaraxes; subplot(2,1,1,pax);
    polarplot(deg2rad(theta(db(RA)>-50)),db(RA(db(RA)>-50)));
    pax.RLim = [-50 0]; pax.ThetaZeroLocation = 'bottom';
    title('Azimuth')
    
    pax = polaraxes; subplot(2,1,2,pax);
    polarplot(deg2rad(theta(db(RE)>-50)),db(RE(db(RE)>-50)));
    pax.RLim = [-50 0]; pax.ThetaZeroLocation = 'bottom';
    title('Elevation')
end
