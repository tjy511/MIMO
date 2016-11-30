% Tracing layers in MIMO in 2D
%
% TJ Young
% 24.10.2016

%% 0. Load imagery file
load('array2d_20140506-1813.mat')
cfg.slice = 'yy'; % 'xx' 'yy'

% Define 3-dimensional variables
xx = xxPix;%repmat([-50:49],641,1);
yy = repmat(Rs',1,100);
if cfg.slice == 'xx'
    zz = pp_slicex;
elseif cfg.slice == 'yy'
    zz = pp_slicey;
end

% Plot unaltered profile
fig1 = plotimgprofile_gland(xx,yy,zz);
%title(['Vertical 2D profile (y-direction) at Date/Time: ', datestr(dateStamp)])

%% 1. Filter profile

cfg.filter = 1; % Switch to filter/smooth profile
cfg.ftype = 'gaussian'; % Rotationally symmetric Gaussian lowpass filter
cfg.fsize = 7; % Filter size
cfg.fsigma = 1; % Standard deviation threshold
cfg.fwindow = 2; % Window size in convolution

% Apply filters to remove noise
if cfg.filter
    filt = (fspecial(cfg.ftype,cfg.fsize,cfg.fsigma)); % Guassian lowpass filter 
    thres = (max([min(max(zz,[],1))  min(max(zz,[],2))])) ;
    zz = medfilt2(zz,[cfg.fwindow,cfg.fwindow]); % Median filtering in 2 directions
    zz = conv2(zz,filt,'same'); % 2-D convolution with designed filter
end

%% 2. Obtain peaks in 2D

% Identify 2D maxima
thresh = -50;
[zmax,imax,~,~]= extrema2(zz);
zmax = zmax(db(zmax) > thresh);
imax = imax(db(zmax) > thresh);

% Convert indexing to [x,y]
[smax.x,smax.y] = ind2sub(size(zz),imax);
%%imax(:,1) = tmp.x; imax(:,2) = tmp.y;

% Plotting fancies
plotimgprofile_gland(xx,yy,zz);
for ii = 1:length(zmax)
    plot3(xx(smax.x(ii),smax.y(ii)),yy(smax.x(ii),smax.y(ii)),zmax(ii),'k.')
end

% Trace reflections through depth (user input)
[ibx,iby] = ginput; % Press return key to exit
close(gcf);
% Convert hand-drawn curve to real data
ibiy = interp1(iby,ibx,yy(:,1));
% Pick nearest points from drawn line
threshx = 3; % +/- bins in range from curve
for ii = 1:length(smax.x)
    int.select(ii) = abs(xx(smax.x(ii),smax.y(ii)) - ibiy(smax.x(ii))) < threshx;
end
plotimgprofile_gland(xx,yy,zz);
for ii = 1:length(int.select)
    if int.select(ii) == 1
        plot3(xx(smax.x(ii),smax.y(ii)),yy(smax.x(ii),smax.y(ii)),zmax(ii),'k.')
    end
end

%% 3. Identify slope of layers
% Slope angle phi is dictated by the amount deviated from nadir.

% 3-dimensional location in Cartesian

for ii = 1:length(int.select)
    tmp.pos(ii) = xx(smax.x(ii),smax.y(ii)); 
    tmp.dep(ii) = yy(smax.x(ii),smax.y(ii));
end

if cfg.slice == 'xx';
    int.x = zeros(size(int.select)); % Location along x-axis
    int.y = tmp.pos; % Location along y-axis
elseif cfg.slice == 'yy';
    int.x = tmp.pos; % Location along x-axis
    int.y = zeros(size(int.select)); % Location along y-axis
end
int.z = tmp.dep; % Location through depth

% 3-dimensional location in Polar
[int.r,int.theta,int.phi] = cart2sph(int.x,int.y,int.z);

%% Export picked slopes

% Subset picked slope
intS.x = int.x(int.select); 
intS.y = int.y(int.select);
intS.z = int.z(int.select);
intS.r = int.r(int.select);
intS.theta = int.theta(int.select);
intS.phi = int.phi(int.select);

ints = [intS.x' intS.y' intS.z'];
ints = sortrows(ints,3);

startup
cd(strcat(rwd,'/results/mimo/'));
save('intSelect','ints');