% Tracing layers in MIMO in 3D
%
% TJ Young
% 10.11.2016

%% Load imagery file and associated parameters
load('array2d_20140506-1813_attended.mat','xxPix','yyPix','imgPlane','dateStamp','R')
close all

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
cfg.pktolm = 5; % 2D tolerance for peaks (bins)
cfg.pkthresh = -50; % dB threshold level for peaks
cfg.pkprom = 0; % Filter by prominence threshold
cfg.phiLim = [225 315]; % Limits for phi (degrees)
cfg.rThresh = 5; % Search range for pkselect; 

% Parameters for plotting
cfg.doPlot = 0; % Turn on intermediate plotting

% Activate config
if cfg.pkselect == 1
    load('intSelect.mat')
end
if cfg.doPlot == 0
    set(0,'DefaultFigureVisible','off')
end

%% Run through layers

clear pr pxy int
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
        filt = (fspecial(cfg.ftype,cfg.fsize,cfg.fsigma)); % Guassian lowpass filter
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
        
        % Limit search to pre-selected peaks
        % Note: Perhaps may be useful to interpolate through identified layers?
        if cfg.pkselect == 1
            iInd = find(abs(ints(:,3)-depth) < cfg.pktolm,1); % Find closest z within 5 m depth
            %rInd = find(abs(int.y-ints(iInd,1))<cfg.rThresh*dy & abs(int.x-ints(iInd,2))<cfg.rThresh*dx);
            rInd = find(abs(int.y-ints(iInd,1))<cfg.rThresh*dy); % Only use layers within range threshold
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
        plotimgdepth_gland(xx,yy,db(zz));
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
    end
    
end

%% Visualise layer in 3 dimensions

plane.x0 = lyr.x;
plane.y0 = lyr.y;
plane.z0 = lyr.z;

set(0,'DefaultFigureVisible','on')

fig = plotlayers_gland(plane.x0,plane.y0,plane.z0,pxy,ppr);
title(['3D layer profile at at Date/Time: ', datestr(dateStamp)])