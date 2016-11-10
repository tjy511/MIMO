% Tracing layers in MIMO in 2D
%
% TJ Young
% 24.10.2016

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
threshx = 5; % +/- bins in range from curve
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
if zz == pp_slicex
    int.x = 0; % Location along x-axis
    int.y = xx(smax.x); % Location along y-axis
elseif zz == pp_slicey
    int.x = xx(smax.x); % Location along x-axis
    int.y = 0; % Location along y-axis
end
int.z = yy(smax.y); % Location through depth

% 3-dimensional location in Polar
int.r = sqrt(int.x.^2+int.y.^2+int.z.^2);
int.theta = acosd(int.z./int.r); 
int.phi = atand(int.y./int.x); 

%% SCRATCH

% % Pick out strongest layers by depth bin
% clear idx
% for ii = 1:size(xx,1)
%     idx.depth = find(smax.x == ii); % Search by depth bin
%     if size(idx.depth,1) == 0
%         smaxd.x(ii) = NaN;
%         smaxd.y(ii) = NaN; 
%         smaxd.z(ii) = NaN; 
%     else
%         [smaxd.z(ii),idx.power] = max(zz(idx.depth));
%         smaxd.x(ii) = smax.x(idx.depth(idx.power));
%         smaxd.y(ii) = smax.y(idx.depth(idx.power));
%     end
% end
% % Remove NaN from indices
% smaxd.x(isnan(smaxd.x)) = [];
% smaxd.y(isnan(smaxd.y)) = [];
% smaxd.z(isnan(smaxd.z)) = [];

% for ii = 1:length(smaxd.x)
%     plot3(xx(smaxd.x(ii),smaxd.y(ii)),yy(smaxd.x(ii),smaxd.y(ii)),smaxd.z(ii),'k.')
% end
%pt = 90; plot3(xx(smax.x(pt),smax.y(pt)),yy(smax.x(pt),smax.y(pt)),zmax(pt),'k*')

% %% 2. Identify layers
% 
% cfg.thresh = -50; % Power threshold for peaks
% cfg.vecLength = 3; % Minimum vector length of layers
% 
% % Pre-allocate arrays
% int.pks = nan(size(xx)); int.locs = int.pks;
% 
% % Find all peaks (1D --> 2D)
% for ii = 1:size(xx,2)
%     vec = zz(:,ii); % Depth vector of power
%     [pks,locs] = findpeaks(vec,'minPeakHeight',db2mag(cfg.thresh)); % Identify peaks
%     for jj = 1:length(pks)
%         int.pks(jj,ii) = pks(jj);
%         int.locs(jj,ii) = locs(jj);
%     end
% end
% 
% % Re-size array to original dimensions
% xxPk = nan(size(xx)); yyPk = xxPk; zzPk = xxPk;
% for ii = 1:size(xx,2)
%     for jj = 1:size(xx,1);
%         if ismember(jj,int.locs(:,ii))
%             xxPk(jj,ii) = xx(jj,ii);
%             yyPk(jj,ii) = yy(jj,ii);
%             zzPk(jj,ii) = zz(jj,ii);
%         end
%     end
% end
% 
% % % Plotting fancies
% % plotimgprofile_gland(xx,yy,zz);
% % plot3(xxPk,yyPk,zzPk,'k.')
% 
% % Identify connected components (layers)
% BW = ~isnan(zzPk);
% CC = bwconncomp(BW); 
% layers.idx = CC.PixelIdxList;
% 
% % To do: Remove vectors shorter than n pixels
% for ii = 1:CC.NumObjects % Replaces short vectors with NaN
%     if size(layers.idx{ii},1) < cfg.vecLength
%         layers.idx{ii} = NaN; 
%     end
% end
% fh = @(x) all(isnan(x(:)));
% layers.idx(cellfun(fh,layers.idx)) = []; % Removes NaN values from cell array
% CC.NumObjects = length(layers.idx); % Re-calculate array length
% 
% % Plotting fancies
% plotimgprofile_gland(xx,yy,zz);
% for ii = 1:CC.NumObjects
%     layer = layers.idx{ii};
%     x = xxPk(layer); y = yyPk(layer); z = zzPk(layer);
%     plot3(x,y,z,'k')
% end
% 
% for ii = 1:CC.NumObjects
%     layer = layers.idx{ii};
%     
%     % Fit regression
%     x = xxPk(layer); y = yyPk(layer); 
%     %breg = robustfit(x,y);
%     breg = regress(y,[ones(length(x),1) x]);
%     
%     % Re-assign output to polar coordinates
%     layers.r(ii) = sqrt(range(x)^2+range(y)^2); % Rho
%     layers.t(ii) = atand(range(y)/range(x)); % Theta
% end
% % May also want to identify distance away from nadir...