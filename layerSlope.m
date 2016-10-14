% Tracing layers in MIMO
%
% TJ Young
% 10.10.2016

%% Load imagery file
load('array2d_20140506-1813_attended.mat')

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
xx = xxPix;%repmat([-50:49],641,1);
yy = repmat(Rs',1,100);
zz = pp_slicey;

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

% Convert indexing to [x,y]
[imax.x,imax.y] = ind2sub(size(zz),imax);
%imax(:,1) = tmp.x; imax(:,2) = tmp.y;

pt = 90;

% Plotting fancies
figure, hold on
surf(xx,yy,db(zz),'EdgeColor','none')
view(180,90)
colormap(jet)
caxis([-100 -20])
for ii = 1:length(zmax)
plot3(xx(imax.x(ii),imax.y(ii)),yy(imax.x(ii),imax.y(ii)),zmax(ii),'k.')
end
plot3(xx(imax.x(pt),imax.y(pt)),yy(imax.x(pt),imax.y(pt)),zmax(pt),'k*')
axis equal

%%
% To do: Obtain ellipses for layers, determine location of maximum axis,
% determine angle from level of this axis. 

% From these 2D maxima, obtain ellipse clouds

% Extract 2D space
sr.x = 10; sr.y = 20;  % Search range (+/- pixels)
wk.r = imax.x(pt)-sr.x:imax.x(pt)+sr.x;
wk.c = imax.y(pt)-sr.y:imax.y(pt)+sr.y;

wk.x = xx(wk.r,wk.c);
wk.y = yy(wk.r,wk.c);
wk.z = zz(wk.r,wk.c);

figure
surf(wk.x,wk.y,db(wk.z),'edgeColor','none')
view(180,90)
colormap(jet)
caxis([-100 -20])
axis equal

% New idea: run peak identification function in 1D vertically, then stitch
% together to make 2D space. 