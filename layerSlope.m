% Tracing layers in MIMO
%
% TJ Young
% 10.10.2016

%% Load imagery file
load('array2d_20140506-1813_attended.mat')

%% Set depth parameters
depth = 600;
depths = 590:620;

%% Plot 2D imagery

% X-axis
figure, hold on
surf(xxPix,Rs,db(pp_slicex),'EdgeColor','none')
view(0,-90)
colormap(jet)
caxis([-100 -20])
legend = colorbar('Ticks',[-100 -80 -60 -40 -20]);
legend.Label.String = 'dB (Vrms)';
xlabel('Range (m)')
ylabel('Depth (m)')
title(['Vertical 2D profile (x-direction) at Date/Time: ', datestr(dateStamp)])

% Y-axis
figure,hold on
surf(yyPix,Rs,db(pp_slicey),'EdgeColor','none')
view(0,-90)
colormap(jet)
caxis([-100 -20])
legend = colorbar('Ticks',[-100 -80 -60 -40 -20]);
legend.Label.String = 'dB (Vrms)';
xlabel('Range (m)')
ylabel('Depth (m)')
title(['Vertical 2D profile (y-direction) at Date/Time: ', datestr(dateStamp)])

%% Identify strong peaks in 2D

% Define variables
xx = repmat([-50:49],641,1);%xxPix;
yy = repmat(Rs',1,100);
zz = pp_slicex;

res = 2; % 1 2

% Apply filters to remove noise
filt = (fspecial('gaussian', 7,1));
thres = (max([min(max(zz,[],1))  min(max(zz,[],2))])) ;
zz = medfilt2(zz,[3,3]);
zz = conv2(single(zz),filt,'same'); 

% Identify 2D maxima
thresh = -50;
[zmax,imax,~,~]= extrema2(zz);
zmax = zmax(db(zmax) > thresh);
imax = imax(db(zmax) > thresh);

% Plotting fancies
figure, hold on
surf(xx,yy,db(zz),'EdgeColor','none')
view(180,90)
colormap(jet)
caxis([-100 -20])

plot3(xx(imax),yy(imax),zmax,'k.')

% To do: Obtain ellipses for layers, determine location of maximum axis,
% determine angle from level of this axis. 