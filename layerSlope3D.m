% Tracing layers in MIMO in 3D
%
% TJ Young
% 10.11.2016

%% 0. Load imagery file and associated parameters
load('array2d_20140506-1813_attended.mat')

% Parameters for depths
depths = [50:50:400];

% Parameters for filter
cfg.filter = 1; % Switch to filter/smooth profile
cfg.ftype = 'gaussian'; % Rotationally symmetric Gaussian lowpass filter
cfg.fsize = 7; % Filter size
cfg.fsigma = 1; % Standard deviation threshold
cfg.fwindow = 4; % Window size in convolution

% Parameters for peak identification
cfg.pkthresh = -50; % dB threshold level for peaks

% Parameters for plotting
cfg.doPlot = 0; % Turn on plotting
if cfg.doPlot == 0
    set(0,'DefaultFigureVisible','off')
end

%% Run through layers

clear pr int
for cc = 1:numel(depths)
    depth = depths(cc);
    
    xx = xxPix(depth,:);
    yy = yyPix(depth,:);
    zz = imgPlane(:,:,depth);
    pxy(cc,:) = xx; % Save 2D scale for plotting
    
    plotimgdepth_gland(xx,yy,db(zz));
    title(['2D depth profile at Depth: ', num2str(depth) ,'m at Date/Time: ', datestr(dateStamp)])
    %view(45,30)
    
    %% 1. Filter slice
    
    % Apply filters to remove noise
    if cfg.filter
        filt = (fspecial(cfg.ftype,cfg.fsize,cfg.fsigma)); % Guassian lowpass filter
        zz = medfilt2(zz,[cfg.fwindow,cfg.fwindow]); % Median filtering in 2 directions
        zz = conv2(zz,filt,'same'); % 2-D convolution with designed filter
    end
    ppr(cc,:,:) = zz; % Save power return for plotting
    
    %% 2. Obtain peaks in 2D
    
    % Identify 2D maxima
    [zmax,imax,~,~]= extrema2(zz);
    zmax = zmax(db(zmax) > cfg.pkthresh);
    imax = imax(db(zmax) > cfg.pkthresh);
    
    % Convert indexing to [x,y]
    [smax.y,smax.x] = ind2sub(size(zz),imax);
    
    % Pick absolute maximum
    [~,smax.mloc] = max(zmax);
    [smax.my,smax.mx] = ind2sub(size(zz),imax(smax.mloc));
    
    % Plotting fancies
    plotimgdepth_gland(xx,yy,db(zz));
    plot3(xx(smax.x),yy(smax.y),zmax,'k.')
    plot3(xx(smax.mx),yy(smax.my),zmax(smax.mloc),'k*')
    
    %% 3. Identify slope of layers
    
    % 3-dimensional location in Cartesian
    int.x = xx(smax.x);
    int.y = yy(smax.y);
    int.z = -depth;
    int.mx(cc) = xx(smax.mx);
    int.my(cc) = yy(smax.my);
    int.mz(cc) = int.z;
    
    % 3-dimensional location in Polar
    int.r = sqrt(int.x.^2+int.y.^2+int.z.^2);
    int.theta = acosd(abs(int.z)./int.r);
    int.phi = atand(int.y./int.x);
    
    % Correct abiguous values of phi (as period of tangent is pi)
    for ii = 1:numel(int.phi)
        if int.x(ii) < 0 && int.y(ii) > 0 % II quadrant
            int.phi(ii) = int.phi(ii) + 180;
        end
        if int.x(ii) < 0 && int.y(ii) < 0 % III quadrant
            int.phi(ii) = int.phi(ii) + 180;
        end
    end
    
    int.mr(cc) = int.r(smax.mloc);
    int.mtheta(cc) = int.theta(smax.mloc);
    int.mphi(cc) = int.phi(smax.mloc);
    
end

%% Visualise layer in 3 dimensions

% Establish plane with normal vector and point
plane.x0 = int.mx;
plane.y0 = int.my;
plane.z0 = int.mz;
plane.a = 0-plane.x0;
plane.b = 0-plane.y0;
plane.c = 0-plane.z0;
plane.d = 0;

% Plotting fancies
set(0,'DefaultFigureVisible','on')
figure, hold on
plot3(0,0,0,'kp') % Location of radar unit
quiver3(0,0,0,-R) % Reference depth vector

clear ax 
for cc = 1:numel(depths)
    % Calculate plane
    plane.x = pxy(cc,:);
    zz = squeeze(ppr(cc,:,:));
    [plane.X,plane.Y] = meshgrid(plane.x);
    plane.Z = (plane.a(cc).*plane.X + plane.b(cc).*plane.Y) ./ -plane.c(cc) + plane.z0(cc);
    
    % Plot layer components iteratively
    quiver3(0,0,0,int.mx(cc)*1.1,int.my(cc)*1.1,int.mz(cc)*1.1); % Vector to identified peak
    plot3(int.mx(cc),int.my(cc),int.mz(cc),'k*') % Location of identified peak
    ax(cc) = patch(surf2patch(plane.X,plane.Y,plane.Z,db(zz)),'lineStyle','none'); % Plot layer with dB as patch
end

legend = colorbar('Ticks',[-100 -80 -60 -40 -20]);
set(ax,'facealpha',0.5)
axis equal
colormap(jet)
caxis([-100 -20])
shading flat
zlim([-depths(end)-50 0])
view(-15,15)
xlabel('x-position (m)'); ylabel('y-position (m)'); zlabel('Depth (m)')