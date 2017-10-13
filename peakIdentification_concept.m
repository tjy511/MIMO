% Tracing layers in MIMO in 2D
%
% TJ Young
% 24.10.2016

%% Config parameters

% File and dimensions
deployment = 1; % 1 2 3 
%cfg.slice = 'yy'; % 'xx' 'yy'

% Filter
cfg.filter = 1; % Switch to filter/smooth profile
cfg.ftype = 'gaussian'; % Rotationally symmetric Gaussian lowpass filter
cfg.fparam = [7 1 2]; % Filter size

% Peak identification
threshx = 3; % +/- bins in range from curve
doPeak = 0; % Re-identify peaks

% Export
doExport = 0; % Export selected points to .mat file
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
load(fileIn,'xxPix','yyPix','imgPlane','pp_slicex','pp_slicey','Rs','dateStamp','R')

% Vertical sections
xxvx = xxPix;
xxvy = zeros(size(xxPix)); 
yyvx = xxvy;
yyvy = xxvx;
zzv = -repmat(Rs',1,100);
ccvx = pp_slicex;
ccvy = pp_slicey;

% Horizontal sections
depth = 250; depth = depth+9;
xxh = xxPix(depth,:);
yyh = yyPix(depth,:);
zzh = repmat(-depth,[size(xxPix,2) size(yyPix,2)]);
cch = abs(imgPlane(:,:,depth));

% Filter profile
if cfg.filter
    ccvx = pkConvol(ccvx,cfg.ftype,cfg.fparam);
    ccvy = pkConvol(ccvy,cfg.ftype,cfg.fparam);
    cch = pkConvol(cch,cfg.ftype,cfg.fparam);
end

%% Plot profile

fig = figure; 
hold on, box on, grid on

hx = surf(xxvx,yyvx,zzv,db(ccvx,'voltage'),'EdgeColor','none');
hy = surf(yyvx,yyvy,zzv,db(ccvy,'voltage'),'EdgeColor','none');
hh = surf(xxh,yyh,zzh,db(cch,'voltage'),'EdgeColor','none');

plot3([0 0],[0 0],[min(zzv(:,1)) max(zzv(:,1))],'k:','lineWidth',1);
plot3([min(xxvx(depth,:)) max(xxvx(depth,:))],[0 0],[-depth -depth],'k:','lineWidth',1);
plot3([0 0],[min(xxvx(depth,:)) max(xxvx(depth,:))],[-depth -depth],'k:','lineWidth',1);

scatter3(0,0,0,50,'k','filled');

ext = 50;
xlim([min(xxvx(depth+ext,:)) max(xxvx(depth+ext,:))]);
ylim([min(xxvx(depth+ext,:)) max(xxvx(depth+ext,:))]);
zlim([zzv(depth+ext,1) 1]);
%xlim([min(xxvx(end,:)) max(xxvx(end,:))]);
%ylim([min(xxvx(end,:)) max(xxvx(end,:))]);
%zlim([min(zzv(:,1)) max(zzv(:,1))]);

view(45,30)
colormap(jet)
caxis([-100 -20])
legend = colorbar('Ticks',[-100 -80 -60 -40 -20]);
legend.Label.String = 'dB (V_{RMS})';

xlabel('Cross-range [m] (X)')
ylabel('Along-range [m] (Y)')
zlabel('Depth [m] (Z)')

%% Export figures

if doSave
    % Create and cd to folder
    fileLoc = '~/Google Drive/Academic/papers/paper3/figs/for_paper/fieldmethods/';
    %fileLoc = '~/Downloads';
    try
        cd(fileLoc);
    catch
        mkdir(fileLoc); cd(fileLoc);
    end
    set(fig,'color','w')
    export_fig(fig,'peak_identification_concept.png','-m6');
end