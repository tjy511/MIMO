% Tracing layers in MIMO in 3D
%
% TJ Young
% 10.11.2016

%% Load imagery file and associated parameters
close all
deployment = 3; % 1 2 3
switch deployment
    case 1
        fileIn = 'array2d_20140506-1813.mat';
        cfg.phiLim = [270-60 270+60]; % Limits for phi (degrees) from x-axis vector [225 315]
        cfg.pkthresh = -50; % dB threshold level for peaks
        cfg.pktolm = 10; % 2D tolerance for peaks (bins in z-direction)
        cfg.rThresh = 15; % Search range for pkselect (m in x-y direction); 
    case 2
        fileIn = 'array2d_20140726-1727.mat';
        cfg.phiLim = [270-45 270+60]; % Limits for phi (degrees) from x-axis vector
        cfg.pkthresh = -60; % dB threshold level for peaks
        cfg.pktolm = 15; % 2D tolerance for peaks (bins in z-direction)
        cfg.rThresh = 15; % Search range for pkselect (bins in x-y direction); 
    case 3
        fileIn = 'array2d_20150703-1221.mat';
        cfg.phiLim = [270-90 270+60]; % Limits for phi (degrees) from x-axis vector [225 315]
        cfg.pkthresh = -50; % dB threshold level for peaks
        cfg.pktolm = 10; % 2D tolerance for peaks (bins in z-direction)
        cfg.rThresh = 25; % Search range for pkselect (bins in x-y direction); 
end
load(fileIn,'xxPix','yyPix','imgPlane','dateStamp','R')

%% Config

% Parameters for filter
cfg.filter = 1; % Switch to filter/smooth profile
cfg.ftype = 'gaussian'; % Rotationally symmetric Gaussian lowpass filter
cfg.fparam = [7 1 2];

% Parameters for peak identification
cfg.pkselect = 1; % Use selected peaks from 2D processing
%cfg.pktolm = 10; % 2D tolerance for peaks (bins in z-direction)
%cfg.pkthresh = -50; % dB threshold level for peaks
cfg.pkprom = 0; % Filter by prominence threshold
%cfg.phiLim = [180-45 180+90]; % Limits for phi (degrees) from x-axis vector [225 315]
%cfg.rThresh = 8; % Search range for pkselect (bins in x-y direction); 

% Parameters for plotting
cfg.doPlot = 0; % Turn on intermediate plotting
cfg.doSave = 1; % Turn on figure exporting

% Activate config
if cfg.pkselect == 1
    intsx = load(strcat('intSelect_',fileIn(9:16),'x.mat'));
    intsy = load(strcat('intSelect_',fileIn(9:16),'y.mat'));
    int1.x = intsx.ints; int1.y = intsy.ints;
end
if cfg.doPlot == 0
    set(0,'DefaultFigureVisible','off')
end

% Parameters for depths
depths = [25:25:450];
%depths = int1.y(:,3);

%% Run through layers

clear pr pxy ppr int iIndx iIndy rIndx rIndy rInd
for cc = 1:numel(depths)
    
    %% Establish structure
    depth = depths(cc);
    xx = xxPix(depth,:);
    yy = yyPix(depth,:);
    zz = abs(imgPlane(:,:,depth));
    
    dx = mean(diff(xx)); dy = mean(diff(yy));
    pxy(cc,:) = xx; % Save 2D scale for plotting
    lyr.t = dateStamp; 
    
    %plotimgdepth_gland(xx,yy,db(zz));
    %title(['2D depth profile at Depth: ', num2str(depth) ,'m at Date/Time: ', datestr(dateStamp)])
    %view(45,30)
    
    %% Filter slice
    
    % Apply filters to remove noise
    if cfg.filter
        zz = pkConvol(zz,cfg.ftype,cfg.fparam);
    end
    ppr(cc,:,:) = zz; % Save power return for plotting
    
    %% Obtain peaks in 2D
    
    % Identify 2D maxima
    [max2.pwr,max2.idx,~,~]= extrema2(zz);
    
    % Set threshold
    max2.pwr = max2.pwr(db(max2.pwr) > cfg.pkthresh);
    max2.idx = max2.idx(db(max2.pwr) > cfg.pkthresh);
    
    % Calculate statistics
    max2.mean = mean(max2.pwr);
    max2.std = std(max2.pwr);
    
    % Convert indexing to [x,y]
    [smax.y,smax.x] = ind2sub(size(zz),max2.idx);
    max2.pwrOrig = max2.pwr; % Create unfiltered max2.pwr var
    
    try
        %% Convert to 3-D referencing
        
        % 3-dimensional location in Cartesian
        int2.x = xx(smax.x);
        int2.y = yy(smax.y);
        int2.z = -depth;
        
        % 3-dimensional location in Polar
        [int2.r,int2.theta,int2.phi] = cart2sph(int2.x,int2.y,int2.z); 
        
        %% Various filtering of peaks 
        
        % Limit search to specified angle
        if isempty(cfg.phiLim) == 0
            if cfg.phiLim(2) < cfg.phiLim(1)
                pInd = find(int2.phi > cfg.phiLim(1) | int2.phi < cfg.phiLim(2)); 
            else
                pInd = find(int2.phi > cfg.phiLim(1) & int2.phi < cfg.phiLim(2));
            end
            int2.x = int2.x(pInd); int2.y = int2.y(pInd);
            int2.r = int2.r(pInd); int2.theta = int2.theta(pInd); int2.phi = int2.phi(pInd);
            max2.idx = max2.idx(pInd); max2.pwr = max2.pwr(pInd);
        end
        
        %% Limit search to pre-selected peaks
        % Note: Perhaps may be useful to interpolate through identified layers?
        if cfg.pkselect == 1
            iIndx = find(abs(int1.x(:,3)-depth) < cfg.pktolm,1); % Find closest z within tolerance
            iIndy = find(abs(int1.y(:,3)-depth) < cfg.pktolm,1); % Find closest z within tolerance
            
            switch deployment % THIS IS TEMPORARY
                case {1,2}
                    rInd = find(abs(int2.y-int1.y(iIndy,1))<cfg.rThresh*dy & abs(int2.x-int1.y(iIndy,2))<cfg.rThresh*dx);
                case {4}
                    rInd = find(abs(int2.x-int1.x(iIndx,1))<cfg.rThresh*dy & abs(int2.x-int1.x(iIndx,2))<cfg.rThresh*dx);
                case {3}
                    try rIndx = find(abs(int2.x-int1.x(iIndx,1))<cfg.rThresh*dx); end; % Only use layers within range threshold
                    try rIndy = find(abs(int2.y-int1.y(iIndy,1))<cfg.rThresh*dy); end % Only use layers within range threshold
                    if exist('rIndx','var') && exist('rIndy','var')
                        rInd = intersect(rIndx,rIndy);
                    elseif exist('rIndx','var') && ~exist('rIndy','var')
                        rInd = rIndx;
                    elseif ~exist('rIndx','var') && exist('rIndy','var')
                        rInd = rIndy;
                    elseif ~exist('rIndx','var') && ~exist('rIndy','var')
                        rInd = [];
                    end
            end
            int2.x = int2.x(rInd); int2.y = int2.y(rInd);
            int2.r = int2.r(rInd); int2.theta = int2.theta(rInd); int2.phi = int2.phi(rInd); 
            max2.idx = max2.idx(rInd); max2.pwr = max2.pwr(rInd); 
        end
        
        %%  Pick absolute maximum
        [~,max2.mloc] = max(max2.pwr);
        [max2.my,max2.mx] = ind2sub(size(zz),max2.idx(max2.mloc));
        max2.mp = max2.pwr(max2.mloc); 
        if isempty(max2.mx) || isempty(max2.my)
            max2.mx = NaN; max2.my = NaN; max2.mp = NaN; 
        end
        
        % Check that identified max is prominent
        if cfg.pkprom == 1
            if max2.mp < max2.mean+1*max2.std
                max2.mx = NaN; max2.my = NaN; max2.mp = NaN;
            end
        end
        
        % Absolute maxima in Cartesian
        lyr.x(cc) = xx(max2.mx);
        lyr.y(cc) = yy(max2.my);
        lyr.z(cc) = int2.z;
        
        % Absolute maxima in Polar
        lyr.r(cc) = int2.r(max2.mloc);
        lyr.theta(cc) = int2.theta(max2.mloc);
        lyr.phi(cc) = int2.phi(max2.mloc);
        
        %% Plotting fancies
        plotimgdepth_gland(xx,yy,zz);
        plot3(0,0,1,'kp')
        plot3(xx(smax.x),yy(smax.y),max2.pwrOrig,'k.','markerSize',5)
        plot3(int2.x,int2.y,max2.pwr,'k.','markerSize',20)
        plot3(xx(max2.mx),yy(max2.my),max2.mp,'k*')
        title(['2D depth profile at Depth: ', num2str(depth) ,'m at Date/Time: ', datestr(dateStamp)])
        
    catch % If any criteria was not satisfied
        int2.x = NaN; int2.y = NaN; int2.z = -depth;
        int2.r = NaN; int2.theta = NaN; int2.phi = NaN;
        lyr.x(cc) = NaN; lyr.y(cc) = NaN; lyr.z(cc) = int2.z;
        lyr.r(cc) = NaN; lyr.theta(cc) = NaN; lyr.phi(cc) = NaN;
        plotimgdepth_gland(xx,yy,zz);
        title(['2D depth profile at Depth: ', num2str(depth) ,'m at Date/Time: ', datestr(dateStamp)])
    end
    
end

%% Manually remove odd layers

switch deployment
    case 1
        %exclude = 13;
    case 2
        %exclude = 6;
    case 3
        exclude = [5 11];
        lyr.x(exclude) = NaN; 
end

%% Visualise layer in 3 dimensions

plane.x0 = lyr.x;
plane.y0 = lyr.y;
plane.z0 = lyr.z;

set(0,'DefaultFigureVisible','on')

%fig = plotlayers_gland(plane.x0,plane.y0,plane.z0,pxy,ppr);
%fig = plotlayers_gland(plane.x0,plane.y0,plane.z0,pxy,lyr.theta);
fig = plotlayers_gland(plane.x0,plane.y0,plane.z0,repmat(pxy(18,:),[18 1]),lyr.theta);
title(['3D layer profile at at Date/Time: ', datestr(dateStamp)])
zLim = -325;
xlim([-150 150])
ylim([-150 150])
%xlim([zLim/2 -zLim/2])
%ylim([zLim/2 -zLim/2])
zlim([zLim 10])
view(-20,20);
box on

set(fig, 'Units','Normalized','Position', [0 0 1/3-0.044 1/3]);

%% Stats
nanmean(lyr.phi)
nanmedian(lyr.phi)

%% Export figures

if cfg.doSave
    % Create and cd to folder
    fileLoc = '~/Google Drive/Academic/papers/paper3/figs/processing/3d/';
    %fileLoc = '~/Downloads';
    try
        cd(fileLoc);
    catch
        mkdir(fileLoc); cd(fileLoc);
    end
    set(fig,'color','w')
    export_fig(fig,strcat(fileIn(9:16),'_3d.png'),'-m6');
end