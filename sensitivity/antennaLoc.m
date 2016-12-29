function [varargout] = antennaLoc(varargin)
% antennaLoc(txrx,quadrant,off,dPhy,doPlot)
%
% Script to set location of antennas.
%
% Notes:
% - Tx runs on the x-axis; Rx runs on the y-axis
% - Tx / Rx is centred about the Cartesian origin
%
% Input variables
% - txrx: # Tx / # Rx [2-element array]
% - quadrant: Location of virtual antennas on the Cartesian system
% - off: Offset of Tx1 / Rx1 [2-element array]
% - dPhy: Physical antenna separation, centre-to-centre [m]
% - doPlot: Switch to plot antenna locations
%
% Output variables
% - tx: Location of Tx real positions
% - rx: Location of Rx real positions 
% - ve: Location of antenna virtual positions 
% 
% Syntax
% antennaLoc('store')
% antennaLoc(txrx,quadrant,off,dPhy)
% antennaLoc(txrx,quadrant,off,dPhy,doPlot)
%
% TJ Young
% 28 December 2016

%% Input variables
if nargin <= 2
    if strcmp(varargin{1},'store1') == 1;
        dPhy = 0.83; off = [2.49 -2.49];
        quadrant = 4; txrx = [8 8];
    end
    if nargin == 2
        doPlot = varargin{2};
    else
        doPlot = 1;
    end
elseif nargin >= 4
    txrx = varargin{1};
    quadrant = varargin{2};
    off = varargin{3};
    dPhy = varargin{4};
    doPlot = 1;
end

if nargin == 5
    doPlot = varargin{5};
end

xOff = off(1); % Offset from x-origin
yOff = off(2); % Offset from y-origin

%% Establish location of antenna real positions

% Establish quadrant
switch quadrant
    case 1
        dPhyx = dPhy;
        dPhyy = dPhy;
    case 2
        dPhyx = -dPhy;
        dPhyy = dPhy;
    case 3
        dPhyx = -dPhy;
        dPhyy = -dPhy;
    case 4
        dPhyx = dPhy;
        dPhyy = -dPhy;
end

% Tx antenna real positions [m]
for ii = 1:txrx(1) % x-coord
    if ii == 1
        tx(1,1) = xOff;
    else
        tx(ii,1) = tx(ii-1,1) + dPhyx;
    end
end
tx(:,2) = 0; % y-coord

% Rx antenna real positions [m]
rx(1:txrx(2),1) = 0; % x-coord
for ii = 1:txrx(2) % y-coord
    if ii == 1
        rx(1,2) = yOff;
    else
        rx(ii,2) = rx(ii-1,2) + dPhyy;
    end
end

%% Establish location of antenna virtual positions

counter = 0;
for ii = 1:txrx(1)
    for jj = 1:txrx(2)
        counter = counter + 1;
        ve(counter,1) = tx(ii,1) + (rx(jj,1)-tx(ii,1))/2;
        ve(counter,2) = rx(jj,2) + (tx(ii,2)-rx(jj,2))/2;
    end
end

%% Plot virtual elements

if doPlot==1
    fig = figure; hold on
    xlabel('x-position [m]')
    ylabel('y-position [m]')
    
    % Plot Tx and Rx antennas
    for k = 1:txrx(1)
        plot(tx(k,1),tx(k,2),'+')
        text(tx(k,1)-0.08,tx(k,2)+0.4,num2str(k))
    end
    
    for k = 1:txrx(2)
        plot(rx(k,1),rx(k,2),'o')
        text(rx(k,1)-0.4,rx(k,2)-0.04,num2str(k))
    end
    
    text(tx(1,1)-dPhy-0.08,tx(1,2)+0.4,'Tx')
    text(rx(1,1)-0.4,rx(1,2)+dPhy+0.04,'Rx')
    
    % Restrict axes
    xlim([floor(min([tx(:,1);rx(:,1)]))-1 ceil(max([tx(:,1);rx(:,1)]))+1])
    ylim([floor(min([tx(:,2);rx(:,2)]))-1 ceil(max([tx(:,2);rx(:,2)]))+1])
    
    % Plot virtual antennas
    for k = 1:size(ve,1)
        plot(ve(k,1),ve(k,2),'rx')
    end
    
end

%% Output variables
if nargout == 1
    varargout{1} = ve; 
elseif nargout == 2
    varargout{1} = ve;
    varargout{2} = fig;
elseif nargout == 3
    varargout{1} = tx;
    varargout{2} = rx; 
    varargout{3} = ve;
elseif nargout == 4
    varargout{1} = tx;
    varargout{2} = rx; 
    varargout{3} = ve;
    varargout{4} = fig;
end