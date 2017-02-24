% Tracing layers in MIMO in 2D
%
% TJ Young
% 24.10.2016

%% Config parameters

% File and dimensions
deployment = 2; % 1 2 3 
cfg.slice = 'xx'; % 'xx' 'yy'

% Filter
cfg.filter = 1; % Switch to filter/smooth profile
cfg.ftype = 'gaussian'; % Rotationally symmetric Gaussian lowpass filter
cfg.fparam = [7 1 2]; % Filter size

% Peak identification
threshx = 3; % +/- bins in range from curve

% Export
doExport = 1; % Export selected points to .mat file
doSave = 1;

%% 0. Load imagery file and associated parameters
close all
switch deployment
    case 1
        fileIn = 'array2d_20140506-1813.mat';
        thresh = -50;
    case 2
        fileIn = 'array2d_20140726-1727.mat';
        thresh = -60;
    case 3
        fileIn = 'array2d_20150703-1221.mat';
        thresh = -60;
end
load(fileIn,'xxPix','yyPix','pp_slicex','pp_slicey','Rs','dateStamp','R')

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
if cfg.filter
    zz = pkConvol(zz,cfg.ftype,cfg.fparam);
end

%% 2. Obtain peaks in 2D

% Identify 2D maxima
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
for ii = 1:length(smax.x)
    int.select(ii) = abs(xx(smax.x(ii),smax.y(ii)) - ibiy(smax.x(ii))) < threshx;
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
int.p = zmax'; % Power of internal layers

% 3-dimensional location in Polar
[int.r,int.theta,int.phi] = cart2sph(int.x,int.y,int.z);

% Subset picked slope
intS.x = int.x(int.select);
intS.y = int.y(int.select);
intS.z = int.z(int.select);
intS.r = int.r(int.select);
intS.theta = int.theta(int.select);
intS.phi = int.phi(int.select);
intS.pow = int.p(int.select); 

ints = [intS.x' intS.y' intS.z' intS.r' intS.theta' intS.phi' intS.pow'];
ints = sortrows(ints,3);

% Remove duplicate values at same depths
uDepths = unique(ints(:,3)); 
for ii = 1:length(uDepths)
    ii = uDepths(ii);
    disp(num2str(ii))
    multIdx = find(ismember(ints(:,3),ii)); % Find multiple values of internal layer
    if numel(multIdx) > 1 % If multiple values exist
        disp(['   ',num2str(ii)])
        [~,idx] = max(ints(multIdx,7));
        multIdx(idx) = [];
        ints(multIdx,:) = [];
    end
end

%% Final product
fig2 = plotimgprofile_gland(xx,yy,zz);
for ii = 1:length(int.select)
    if int.select(ii) == 1
        if cfg.slice == 'xx'
            plot3(ints(:,2),ints(:,3),db(ints(:,7),'voltage')-10,'k.', 'markerSize',8)
        elseif cfg.slice == 'yy'
            plot3(ints(:,1),ints(:,3),db(ints(:,7),'voltage')-10,'k.', 'markerSize',8)
        end
    end
end

%% Export picked slopes

if doExport
    startup
    fileOut = strcat('intSelect_',fileIn(9:16),cfg.slice(1));
    cd(strcat(rwd,'/results/mimo/'));
end
save(fileOut,'ints');

%% Export figures

if doSave
    % Create and cd to folder
    fileLoc = '~/Google Drive/Academic/papers/paper3/figs/2d/';
    %fileLoc = '~/Downloads';
    try
        cd(fileLoc);
    catch
        mkdir(fileLoc); cd(fileLoc);
    end
    set([fig1 fig2],'color','w')
    export_fig(fig1,strcat(fileIn(9:16),'_slice',cfg.slice(1),'_orig.png'),'-m2');
    export_fig(fig2,strcat(fileIn(9:16),'_slice',cfg.slice(1),'_layers_select.png'),'-m2');
end