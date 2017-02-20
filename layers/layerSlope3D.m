% Tracing layers in MIMO in 3D
%
% TJ Young
% 10.11.2016

%% Load imagery file and associated parameters
close all
deployment = 2; % 1 2 3
switch deployment
    case 1
        fileIn = 'array2d_20140506-1813.mat';
        cfg.phiLim = [90-45 90+60]; % Limits for phi (degrees) from x-axis vector [225 315]
        cfg.pkthresh = -50; % dB threshold level for peaks
        cfg.pktolm = 10; % 2D tolerance for peaks (bins in z-direction)
        cfg.rThresh = 10; % Search range for pkselect (m in x-y direction); 
    case 2
        fileIn = 'array2d_20140726-1727.mat';
        cfg.phiLim = [90 180]; % Limits for phi (degrees) from x-axis vector
        cfg.pkthresh = -60; % dB threshold level for peaks
        cfg.pktolm = 20; % 2D tolerance for peaks (bins in z-direction)
        cfg.rThresh = 15; % Search range for pkselect (bins in x-y direction); 
    case 3
        fileIn = 'array2d_20150703-1221.mat';
        cfg.phiLim = [180-45 180+90]; % Limits for phi (degrees) from x-axis vector [225 315]
        cfg.pkthresh = -50; % dB threshold level for peaks
        cfg.pktolm = 10; % 2D tolerance for peaks (bins in z-direction)
        cfg.rThresh = 8; % Search range for pkselect (bins in x-y direction); 
end
load(fileIn,'xxPix','yyPix','imgPlane','dateStamp','R')

% Parameters for depths
depths = [25:25:400];

%% Config

% Parameters for filter
cfg.filter = 1; % Switch to filter/smooth profile
cfg.ftype = 'gaussian'; % Rotationally symmetric Gaussian lowpass filter
cfg.fsize = 7; % Filter size
cfg.fsigma = 1; % Standard deviation threshold
cfg.fwindow = 4; % Window size in convolution

% Parameters for peak identification
cfg.pkselect = 1; % Use selected peaks from 2D processing
%cfg.pktolm = 10; % 2D tolerance for peaks (bins in z-direction)
%cfg.pkthresh = -50; % dB threshold level for peaks
cfg.pkprom = 0; % Filter by prominence threshold
%cfg.phiLim = [180-45 180+90]; % Limits for phi (degrees) from x-axis vector [225 315]
%cfg.rThresh = 8; % Search range for pkselect (bins in x-y direction); 

% Parameters for plotting
cfg.doPlot = 1; % Turn on intermediate plotting
cfg.doSave = 1; % Turn on figure exporting

% Activate config
if cfg.pkselect == 1
    intsx = load(strcat('intSelect_',fileIn(9:16),'x.mat'));
    intsy = load(strcat('intSelect_',fileIn(9:16),'y.mat'));
    intsx = intsx.ints; intsy = intsy.ints;
end
if cfg.doPlot == 0
    set(0,'DefaultFigureVisible','off')
end

%% Run through layers

clear pr pxy ppr int iIndx iIndy rIndx rIndy rInd
for cc = 1:numel(depths)
    
    %% Establish structure
    depth = depths(cc);
    xx = xxPix(depth,:);
    yy = yyPix(depth,:);
    zz = imgPlane(:,:,depth);
    
    dx = mean(diff(xx)); dy = mean(diff(yy));
    pxy(cc,:) = xx; % Save 2D scale for plotting
    lyr.t = dateStamp; 
    
    %plotimgdepth_gland(xx,yy,db(zz));
    %title(['2D depth profile at Depth: ', num2str(depth) ,'m at Date/Time: ', datestr(dateStamp)])
    %view(45,30)
    
    %% Filter slice
    
    % Apply filters to remove noise
    if cfg.filter
        filt = (fspecial(cfg.ftype,cfg.fsize,cfg.fsigma)); % Gaussian lowpass filter
        zz = medfilt2(zz,[cfg.fwindow,cfg.fwindow]); % Median filtering in 2 directions
        zz = conv2(zz,filt,'same'); % 2-D convolution with designed filter
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
  %%  
    try
        %% Convert to 3-D referencing
        
        % 3-dimensional location in Cartesian
        int.x = xx(smax.x);
        int.y = yy(smax.y);
        int.z = -depth;
        
        % 3-dimensional location in Polar
        [int.r,int.theta,int.phi] = cart2sph(int.x,int.y,int.z); 
        
        %% Various filtering of peaks 
        
        % Limit search to specified angle
        if isempty(cfg.phiLim) == 0
            pInd = find(int.phi > cfg.phiLim(1) & int.phi < cfg.phiLim(2));
            int.x = int.x(pInd); int.y = int.y(pInd);
            int.r = int.r(pInd); int.theta = int.theta(pInd); int.phi = int.phi(pInd);
            max2.idx = max2.idx(pInd); max2.pwr = max2.pwr(pInd);
        end
        
        %% Limit search to pre-selected peaks
        % Note: Perhaps may be useful to interpolate through identified layers?
        if cfg.pkselect == 1
            iIndx = find(abs(intsx(:,3)-depth) < cfg.pktolm,1); % Find closest z within tolerance
            iIndy = find(abs(intsy(:,3)-depth) < cfg.pktolm,1); % Find closest z within tolerance
            
            switch deployment % THIS IS TEMPORARY
                case 1
                    rInd = find(abs(int.y-intsy(iIndy,1))<cfg.rThresh*dy & abs(int.x-intsy(iIndy,2))<cfg.rThresh*dx);
                case {2,3}
                    try rIndx = find(abs(int.x-intsx(iIndx,1))<cfg.rThresh*dx); end; % Only use layers within range threshold
                    try rIndy = find(abs(int.y-intsy(iIndy,1))<cfg.rThresh*dy); end % Only use layers within range threshold
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
            int.x = int.x(rInd); int.y = int.y(rInd);
            int.r = int.r(rInd); int.theta = int.theta(rInd); int.phi = int.phi(rInd); 
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
        lyr.z(cc) = int.z;
        
        % Absolute maxima in Polar
        lyr.r(cc) = int.r(max2.mloc);
        lyr.theta(cc) = int.theta(max2.mloc);
        lyr.phi(cc) = int.phi(max2.mloc);
        
        %% Plotting fancies
        plotimgdepth_gland(xx,yy,zz);
        plot3(0,0,1,'kp')
        plot3(xx(smax.x),yy(smax.y),max2.pwrOrig,'k.','markerSize',5)
        plot3(int.x,int.y,max2.pwr,'k.','markerSize',20)
        plot3(xx(max2.mx),yy(max2.my),max2.mp,'k*')
        title(['2D depth profile at Depth: ', num2str(depth) ,'m at Date/Time: ', datestr(dateStamp)])
        
    catch % If any criteria was not satisfied
        int.x = NaN; int.y = NaN; int.z = -depth;
        int.r = NaN; int.theta = NaN; int.phi = NaN;
        lyr.x(cc) = NaN; lyr.y(cc) = NaN; lyr.z(cc) = int.z;
        lyr.r(cc) = NaN; lyr.theta(cc) = NaN; lyr.phi(cc) = NaN;
        plotimgdepth_gland(xx,yy,zz);
        title(['2D depth profile at Depth: ', num2str(depth) ,'m at Date/Time: ', datestr(dateStamp)])
    end
    
end

%% Manually remove odd layers

switch deployment
    case 1
        exclude = 13;
    case 2
        
    case 3
end
lyr.x(exclude) = NaN; 

%% Visualise layer in 3 dimensions

plane.x0 = lyr.x;
plane.y0 = lyr.y;
plane.z0 = lyr.z;

set(0,'DefaultFigureVisible','on')

%fig = plotlayers_gland(plane.x0,plane.y0,plane.z0,pxy,ppr);
fig = plotlayers_gland(plane.x0,plane.y0,plane.z0,pxy,lyr.theta);
title(['3D layer profile at at Date/Time: ', datestr(dateStamp)])
zLim = -325;
xlim([zLim/2 -zLim/2])
ylim([zLim/2 -zLim/2])
zlim([zLim 0])

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
    export_fig(fig,strcat(fileIn(9:16),'_3d.png'),'-m2');
end