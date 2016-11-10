% Tracing layers in MIMO in 3D
%
% TJ Young
% 10.11.2016

%% 0. Load imagery file
load('array2d_20140506-1813_attended.mat')

depth = 200; 
xx = xxPix(depth,:);
yy = yyPix(depth,:); 
zz = imgPlane(:,:,depth);

plotimgdepth_gland(xx,yy,db(zz))
title(['2D depth profile at Depth: ', num2str(depth) ,'m at Date/Time: ', datestr(dateStamp)])
%view(45,30)

%% 1. Filter slice

cfg.filter = 1; % Switch to filter/smooth profile
cfg.ftype = 'gaussian'; % Rotationally symmetric Gaussian lowpass filter
cfg.fsize = 7; % Filter size
cfg.fsigma = 1; % Standard deviation threshold
cfg.fwindow = 4; % Window size in convolution

% Apply filters to remove noise
if cfg.filter
    filt = (fspecial(cfg.ftype,cfg.fsize,cfg.fsigma)); % Guassian lowpass filter 
    thres = (max([min(max(zz,[],1))  min(max(zz,[],2))])) ;
    zz = medfilt2(zz,[cfg.fwindow,cfg.fwindow]); % Median filtering in 2 directions
    zz = conv2(zz,filt,'same'); % 2-D convolution with designed filter
end

plotimgdepth_gland(xx,yy,db(zz))
%view(45,30)

%% 2. Obtain peaks in 2D

% Identify 2D maxima
thresh = -50;
[zmax,imax,~,~]= extrema2(zz);
zmax = zmax(db(zmax) > thresh);
imax = imax(db(zmax) > thresh);

% Convert indexing to [x,y]
[smax.y,smax.x] = ind2sub(size(zz),imax);

% Pick absolute maximum
[~,smax.mloc] = max(zmax);
[smax.my,smax.mx] = ind2sub(size(zz),imax(smax.mloc));

% Plotting fancies
plotimgdepth_gland(xx,yy,db(zz))
plot3(xx(smax.x),yy(smax.y),zmax,'k.')
plot3(xx(smax.mx),yy(smax.my),zmax(smax.mloc),'k*')

%% 3. Identify slope of layers

% 3-dimensional location in Cartesian
int.x = xx(smax.x);
int.y = yy(smax.y);
int.z = depth;
int.mx = xx(smax.mx);
int.my = yy(smax.my);
int.mz = int.z; 

% 3-dimensional location in Polar
int.r = sqrt(int.x.^2+int.y.^2+int.z.^2);
int.theta = acosd(int.z./int.r);
int.phi = atand(int.y./int.x);
int.mr = int.r(smax.mloc);
int.mtheta = int.theta(smax.mloc);
int.mphi = int.phi(smax.mloc);

% TO FIX: Tangent period is pi so int.phi in 3rd and 4th quadrant are
% ambiguous to 1st and 2nd quadrants
