% Tracing layers in MIMO
%
% TJ Young
% 10.10.2016

%% 0. Load imagery file
load('array2d_20140506-1813_attended.mat')

% Define 3-dimensional variables
xx = xxPix;%repmat([-50:49],641,1);
yy = repmat(Rs',1,100);
zz = pp_slicey;

% Plot unaltered profile
plotimgprofile_gland(xx,yy,zz);
title(['Vertical 2D profile (y-direction) at Date/Time: ', datestr(dateStamp)])

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

%% 2. Identify layers

cfg.thresh = -50; % Power threshold for peaks
cfg.vecLength = 3; % Minimum vector length of layers

% Pre-allocate arrays
int.pks = nan(size(xx)); int.locs = int.pks;

% Find all peaks (1D --> 2D)
for ii = 1:size(xx,2)
    vec = zz(:,ii); % Depth vector of power
    [pks,locs] = findpeaks(vec,'minPeakHeight',db2mag(cfg.thresh)); % Identify peaks
    for jj = 1:length(pks)
        int.pks(jj,ii) = pks(jj);
        int.locs(jj,ii) = locs(jj);
    end
end

% Re-size array to original dimensions
xxPk = nan(size(xx)); yyPk = xxPk; zzPk = xxPk;
for ii = 1:size(xx,2)
    for jj = 1:size(xx,1);
        if ismember(jj,int.locs(:,ii))
            xxPk(jj,ii) = xx(jj,ii);
            yyPk(jj,ii) = yy(jj,ii);
            zzPk(jj,ii) = zz(jj,ii);
        end
    end
end

% % Plotting fancies
% plotimgprofile_gland(xx,yy,zz);
% plot3(xxPk,yyPk,zzPk,'k.')

% Identify connected components (layers)
BW = ~isnan(zzPk);
CC = bwconncomp(BW); 
layers.idx = CC.PixelIdxList;

% To do: Remove vectors shorter than n pixels
for ii = 1:CC.NumObjects % Replaces short vectors with NaN
    if size(layers.idx{ii},1) < cfg.vecLength
        layers.idx{ii} = NaN; 
    end
end
fh = @(x) all(isnan(x(:)));
layers.idx(cellfun(fh,layers.idx)) = []; % Removes NaN values from cell array
CC.NumObjects = length(layers.idx); % Re-calculate array length

% Plotting fancies
plotimgprofile_gland(xx,yy,zz);
for ii = 1:CC.NumObjects
    layer = layers.idx{ii};
    x = xxPk(layer); y = yyPk(layer); z = zzPk(layer);
    plot3(x,y,z,'k')
end

%% 3. Identify slope of layers

for ii = 1:CC.NumObjects
    layer = layers.idx{ii};
    
    % Fit regression
    x = xxPk(layer); y = yyPk(layer); 
    %breg = robustfit(x,y);
    breg = regress(y,[ones(length(x),1) x]);
    
    % Re-assign output to polar coordinates
    layers.r(ii) = sqrt(range(x)^2+range(y)^2); % Rho
    layers.t(ii) = atand(range(y)/range(x)); % Theta
end

%% SCRATCH

% Identify 2D maxima
thresh = -50;
[zmax,imax,~,~]= extrema2(zz);
zmax = zmax(db(zmax) > thresh);
imax = imax(db(zmax) > thresh);

% Convert indexing to [x,y]
[imax.x,imax.y] = ind2sub(size(zz),imax);
%imax(:,1) = tmp.x; imax(:,2) = tmp.y;

pt = 90;

% Plotting fancies
plotimgprofile_gland(xx,yy,zz);
for ii = 1:length(zmax)
plot3(xx(imax.x(ii),imax.y(ii)),yy(imax.x(ii),imax.y(ii)),zmax(ii),'k.')
end
plot3(xx(imax.x(pt),imax.y(pt)),yy(imax.x(pt),imax.y(pt)),zmax(pt),'k*')

