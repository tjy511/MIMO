% Tracing layers in MIMO in 4D
%
% TJ Young
% 10.11.2016

%% Set up script

% cd to file
startup
cd(strcat(rwd,'data/process/mimo/unattended/'))

% List files to load
dirList = dir('*.mat');
for ii = 1:length(dirList)
    fileList{ii,1} = dirList(ii).name;
end

% Run initial script
layerSlope3D
lyr0 = lyr; clear lyr

% Establish depth vector
depths = [25:25:400];

%% Activate config

% Parameters for filter
cfg.filter = 1; % Switch to filter/smooth profile
cfg.ftype = 'gaussian'; % Rotationally symmetric Gaussian lowpass filter
cfg.fsize = 7; % Filter size
cfg.fsigma = 1; % Standard deviation threshold
cfg.fwindow = 4; % Window size in convolution

% Parameters for peak identification
cfg.pkselect = 1; % Use selected peaks from 2D processing
cfg.pkthresh = -55; % dB threshold level for peaks
cfg.pkprom = 0; % Filter by prominence threshold
cfg.phiLim = [225 315]; % Limits for phi (degrees)
cfg.rThresh = 5; % Search range for pkselect; 

% Parameters for display
cfg.verbose = 1; % Turn on display text
cfg.doPlot = 0; % Turn on intermediate plotting

cfg.doSave = 1; % Save images

% Activate config
if cfg.pkselect == 1
    load('intSelect.mat')
end
if cfg.doPlot == 0
    set(0,'DefaultFigureVisible','off')
end

%% Loop through files
for file = 1:numel(fileList)
    
    %% Set up file structure
    
    % Load file
    load(fileList{file},'xxPix','yyPix','imgPlane','dateStamp','R')
    disp(['Processing file: ',fileList{file}])
    lyr.t(file) = dateStamp;
    clear pr pxy int
    for cc = 1:numel(depths)
        
        % Assign parameters
        depth = depths(cc);
        xx = xxPix(depth,:);
        yy = yyPix(depth,:);
        zz = imgPlane(:,:,depth);
        
        dx = mean(diff(xx)); dy = mean(diff(yy));
        pxy(cc,:) = xx; % Save 2D scale for plotting
        
        %% Filter slice
        
        % Apply filters to remove noise
        if cfg.filter
            filt = (fspecial(cfg.ftype,cfg.fsize,cfg.fsigma)); % Gaussian lowpass filter
            zz = medfilt2(zz,[cfg.fwindow,cfg.fwindow]); % Median filtering in 2 directions
            zz = conv2(zz,filt,'same'); % 2-D convolution with designed filter
        end
        
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
        
        %% Convert to 3-D referencing
        try
            % 3-dimensional location in Cartesian
            int.x = xx(smax.x);
            int.y = yy(smax.y);
            int.z = -depth;
            
            % 3-dimensional location in Polar
            [int.r,int.theta,int.phi] = cart2sph(int.x,int.y,int.z);
            
            % Find nearest point within threshold
            trk.xythresh = sqrt((5*dx)^2+(5*dy)^2);
            trk.opt = sqrt((lyr0.x(cc)-int.x).^2 + (lyr0.y(cc)-int.y).^2);
            [trk.d,trk.i] = nanmin(trk.opt);
            if isnan(trk.d), trk.i = NaN; end
            
            % Check if layer fits threshold
            if trk.d > trk.xythresh
                warning(['layer: ',num2str(depth),'m fails threshold'])
            end
            
            % Reassign parameters from layer
            lyr.x(file,cc) = int.x(trk.i);
            lyr.y(file,cc) = int.y(trk.i);
            lyr.z(file,cc) = int.z;
            lyr.r(file,cc) = int.r(trk.i);
            lyr.theta(file,cc) = int.theta(trk.i);
            lyr.phi(file,cc) = int.phi(trk.i);
            lyr.dxy(file,cc) = trk.d;
            
            Disp(['Processed layer for depth:',num2str(depth),'m'])
            
            % Plotting fancies
            plotimgdepth_gland(xx,yy,db(zz));
            plot3(0,0,1,'kp')
            plot3(xx(smax.x),yy(smax.y),max2.pwrOrig,'k.','markerSize',5)
            title(['2D depth profile at Depth: ', num2str(depth) ,'m at Date/Time: ', datestr(dateStamp)])
            
        catch
            % Reassign parameters as NaN
            lyr.x(file,cc) = NaN;
            lyr.y(file,cc) =NaN;
            lyr.z(file,cc) = int.z;
            lyr.r(file,cc) = NaN;
            lyr.theta(file,cc) = NaN;
            lyr.phi(file,cc) = NaN;
            lyr.dxy(file,cc) = trk.d;
            
            %Disp(['Failed to obtain layer at depth: ',num2str(depth),'m'])
        end
    end
    
    % Renew reference file
    if file > 1
        lyr0.t = lyr.t(file-1);
        lyr0.x = lyr.x(file-1,:);
        lyr0.y = lyr.y(file-1,:);
        lyr0.z = lyr.z(file-1,:);
        lyr0.r = lyr.r(file-1,:);
        lyr0.theta = lyr.theta(file-1,:);
        lyr0.phi = lyr.phi(file-1,:);
    end
end

%% Plot figures

% Establish plane with normal vector and point
plane.x0 = lyr.x;
plane.y0 = lyr.y;
plane.z0 = lyr.z;
plane.a = 0-plane.x0;
plane.b = 0-plane.y0;
plane.c = 0-plane.z0;

set(0,'DefaultFigureVisible','on')
% Plotting fancies
for ff = 1:size(lyr.x,1)
    fig(ff) = figure; hold on, grid on
    plot3(0,0,0,'kp') % Location of radar unit
    quiver3(0,0,0,-R,'Color',[0.6 0.6 0.6],'lineWidth',1.5) % Reference depth vector
    
    clear ax
    for cc = 1:numel(depths)
        % Calculate plane
        plane.x = pxy(cc,:);
        [plane.X,plane.Y] = meshgrid(plane.x);
        plane.Z = (plane.a(ff,cc).*plane.X + plane.b(ff,cc).*plane.Y) ./ -plane.c(ff,cc) + plane.z0(ff,cc);
        
        % Plot layer components iteratively
        quiver3(0,0,0,lyr.x(ff,cc)*1.1,lyr.y(ff,cc)*1.1,lyr.z(ff,cc)*1.1,'k'); % Vector to identified peak
        plot3(lyr.x(ff,cc),lyr.y(ff,cc),lyr.z(ff,cc),'k*') % Location of identified peak
        ax(cc) = surf(plane.X,plane.Y,plane.Z);% Plot layer with uniform patch
        %patch(surf2patch(plane.X,plane.Y,plane.Z,db(zz)),'lineStyle','none'); % Plot layer with dB as patch
    end
    set(ax,'facealpha',0.5)
    shading flat
    legend = colorbar;
    legend.Label.String = 'Depth (m)';
    caxis([-depths(end) 0])
    zlim([-depths(end)-50 0])
    view(-15,15)
    xlabel('x-position (m)'); ylabel('y-position (m)'); zlabel('Depth (m)')
    title(['3D layer profile at at Date/Time: ', datestr(lyr.t(ff))])
    
    % Save image
    if cfg.doSave
        saveas(fig(ff),strcat(),'jpg')
    end
end
