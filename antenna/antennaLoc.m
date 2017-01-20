function [varargout] = antennaLoc(type,varargin)
% antennaLoc(txrx,quadrant,off,dPhy,doPlot)
%
% Script to set location of antennas.
%
% Notes:
% - Tx runs on the x-axis; Rx runs on the y-axis
% - Tx / Rx is centred about the Cartesian origin
%
% Input variables
% - type: Type of generation
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

% antennaLoc('store1')
% antennaLoc('real',txrx,dPhy,quadrant,off)
% antennaLoc('real',txrx,dPhy,quadrant,off,doPlot)
% antennaLoc('virtual',txrx,dPhy)
% antennaLoc('virtual',txrx,dPhy,doPlot)
%
% TJ Young
% 28 December 2016

%% Input variables

switch type
    case 'store1'
        dPhy = 0.83; off = (2350/sqrt(2) + 830/2)/1000;
        quadrant = 2; txrx = [8 8];
        if nargin == 1
            doPlot = varargin{1};
        else
            doPlot = 1;
        end
    case 'store2'
        dPhy = 0.83; off = (2700/sqrt(2) + 830/2)/1000;
        quadrant = 3; txrx = [8 8];
        if nargin == 1
            doPlot = varargin{1};
        else
            doPlot = 1;
        end
    case 'store3'
        dPhy = 0.83; off = (2480/sqrt(2) + 830/2)/1000;
        quadrant = 2; txrx = [8 8];
        if nargin == 1
            doPlot = varargin{1};
        else
            doPlot = 1;
        end
    case 'real'
        txrx = varargin{1};
        dPhy = varargin{2};
        quadrant = varargin{3};
        off = varargin{4};
        if nargin == 6
            doPlot = varargin{5};
        else
            doPlot = 1;
        end
        
    case 'virtual'
        txrx = varargin{1};
        dPhy = varargin{2};
        quadrant = 0; 
        off = [0 0]; 
        if nargin == 4
            doPlot = varargin{3};
        else
            doPlot = 1;
        end
end

%% Establish location of antenna real positions

% Establish quadrant
switch quadrant
    case 0
        dPhyx = 2*dPhy; dPhyy = 2*dPhy;
        xOff = off; yOff = off; 
    case 1
        dPhyx = dPhy; dPhyy = dPhy;
        xOff = off; yOff = off; 
    case 2
        dPhyx = -dPhy; dPhyy = dPhy;
        xOff = -off; yOff = off; 
    case 3
        dPhyx = -dPhy; dPhyy = -dPhy;
        xOff = -off; yOff = -off; 
    case 4
        dPhyx = dPhy; dPhyy = -dPhy;
        xOff = off; yOff = -off; 
end

% xOff = off(1); % Offset from x-origin
% yOff = off(2); % Offset from y-origin

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


if doPlot == 0
    set(0,'DefaultFigureVisible','off')
end

fig = figure; hold on, axis equal
xlabel('x-position [m]')
ylabel('y-position [m]')

% Plot Tx and Rx antennas
if strcmp(type,'virtual')
    % Restrict axes
    xlim([floor(min([ve(:,1);ve(:,1)]))-1 ceil(max([ve(:,1);ve(:,1)]))+1])
    ylim([floor(min([ve(:,2);ve(:,2)]))-1 ceil(max([ve(:,2);ve(:,2)]))+1])
else
    for k = 1:txrx(1)
        plot(tx(k,1),tx(k,2),'+')
        text(tx(k,1)-0.08,tx(k,2)+0.4,num2str(k))
    end
    
    for k = 1:txrx(2)
        plot(rx(k,1),rx(k,2),'o')
        text(rx(k,1)-0.4,rx(k,2)-0.04,num2str(k))
    end
    
    switch quadrant
        case 1
    text(tx(1,1)-dPhy-0.08,tx(1,2)+0.4,'Tx')
    text(rx(1,1)-0.4,rx(1,2)-dPhy-0.04,'Rx')
            case 2
    text(tx(1,1)+dPhy+0.08,tx(1,2)+0.4,'Tx')
    text(rx(1,1)-0.4,rx(1,2)-dPhy-0.04,'Rx')
            case 3
    text(tx(1,1)+dPhy+0.08,tx(1,2)+0.4,'Tx')
    text(rx(1,1)-0.4,rx(1,2)+dPhy+0.04,'Rx')
            case 4
    text(tx(1,1)-dPhy-0.08,tx(1,2)+0.4,'Tx')
    text(rx(1,1)-0.4,rx(1,2)+dPhy+0.04,'Rx')
    end
    
    % Restrict axes
    xlim([floor(min([tx(:,1);rx(:,1)]))-1 ceil(max([tx(:,1);rx(:,1)]))+1])
    ylim([floor(min([tx(:,2);rx(:,2)]))-1 ceil(max([tx(:,2);rx(:,2)]))+1])
end

% Plot virtual antennas
for k = 1:size(ve,1)
    plot(ve(k,1),ve(k,2),'rx')
end

set(0,'DefaultFigureVisible','on')

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